#pragma once
#include "libload.h"
#include "options.h"
#include <array>
#include <cstdint>
#include <deque>
#include <string_view>
void trim(string& str,const std::string& whitespace = " \t\r\n");
bool isGZfile(const string fi);

// Input streams read in large blocks; small reads are slow on network and cluster file systems.
inline constexpr size_t __inputBufferSize = size_t(1) << 20;


// Rank names, each stored once at an address that never changes; ID 0 is the unknown name. One pool
// serves every reference database of a run, so assignments keep their names after a RefTax is gone.
// intern() is not thread-safe; name() is, while no thread interns.
class NamePool
{
public:
	NamePool() { names.push_back(&__unkwnTax); }
	NamePool(const NamePool&) = delete;            // the index refers into this object
	NamePool& operator=(const NamePool&) = delete;
	uint32_t intern(std::string_view name);
	const string& name(uint32_t id) const { return *names[id]; }
	size_t size() const { return names.size(); }
private:
	std::deque<string> store;
	unordered_map<std::string_view, uint32_t> index;
	vector<const string*> names;                   // by ID
};


struct TaxObj
{
	TaxObj(std::string_view, int, bool nativeSLV, bool doNotCheckTax, NamePool&);
	TaxObj(const TaxObj* t);
	TaxObj(int d) : SavedTaxs(), Subj(""), hitDB(""), perID(0.f), repID(false),
		hasHitDB(false), speciesUncertain(false), depth(0) { SavedTaxs.reserve(d); }
	//functions
	string getWriteString(const vector<double>&);
	void limitByIdentity(double identity, const vector<double>& thresholds) {
		int eligible = 0;
		while (eligible < depth && eligible < (int)thresholds.size() && thresholds[eligible] <= identity) {
			++eligible;
		}
		depth = eligible;
	}
	void setRepID(bool x) { repID = x; }
	//get tax at depth x
	const string& get(int x) const {
		if (x < 0 || x >= depth || x >= (int)SavedTaxs.size()) { return __unkwnTax; }
		return *SavedTaxs[x];
	}
	// name must outlive this object: a NamePool's name, or __unkwnTax.
	void set(int x, const string* name) {
       if (x < 0) { return; }
		if (x >= (int)SavedTaxs.size()) { SavedTaxs.resize(x + 1, &__unkwnTax); }
		SavedTaxs[x] = name;
		if (x >= depth) { depth = x + 1; }
	}
	//check if other tax is better and copies if so these vals over itself
	bool evalAcpyTax(const TaxObj* oth);
	void copyOver(const TaxObj* oth);
	void setHitDB(const string& x) { hitDB = x; hasHitDB = true; }
	const string& getHitDB() const { return hitDB; }
	bool reportsHitDB() const { return hasHitDB; }

	//variables
	vector<const string*> SavedTaxs; // rank names, owned by a NamePool
	string Subj;
	string hitDB;
	double perID;
	bool repID;
	bool hasHitDB;
	bool speciesUncertain;
	int depth;//number of taxonomic levels assigned/available; not the configured maximum

};


// Reference taxonomy. Rank names are stored once, in a NamePool that must outlive it; a reference's
// lineage is a run of name IDs (0 = unknown), and subjects are found through an open-addressing hash table.
class RefTax
{
public:
	struct Lineage {
		uint64_t subjectOffset; // into subjectChars
		uint32_t subjectLength;
		uint32_t firstName;     // into nameIds; depth entries follow
		uint16_t depth;
		bool speciesUncertain;
	};

	static constexpr size_t defaultBlockSize = size_t(4) << 20;
	// threads parse the file in parallel; blockSize is the amount of text handed to a thread at once. With useIndex,
	// the taxonomy is read from its index (indexPath) if that was made from the same file with the same settings, and
	// the index is written otherwise.
	RefTax(const string&, NamePool&, int tdep, bool, bool, int threads = 1, size_t blockSize = defaultBlockSize,
		bool useIndex = false);
	RefTax(const RefTax&) = delete;            // lineages and the subject index are large
	RefTax& operator=(const RefTax&) = delete;
	static string indexPath(const string& taxonomyFile) { return taxonomyFile + ".idx"; }
	void stats();
	// The reference's lineage, or nullptr for an unknown subject.
	const Lineage* find(std::string_view subject) const;
	uint32_t nameId(const Lineage& lineage, int rank) const {
		return rank < lineage.depth ? nameIds[lineage.firstName + rank] : 0;
	}
	const string& name(uint32_t id) const { return pool.name(id); }
	const string& rankName(const Lineage& lineage, int rank) const { return name(nameId(lineage, rank)); }
	std::string_view subject(const Lineage& lineage) const {
		return std::string_view(subjectChars).substr(lineage.subjectOffset, lineage.subjectLength);
	}
	size_t size() const { return uniqueSubjects; }
	size_t distinctNames() const { return pool.size() - 1; } // in the pool, so of every taxonomy loaded into it
	void setTaxLvls(vector<string> x) { tlevels = x; }
private:
	struct ParsedBlock;
	struct IndexHeader;
	void readText(const string& inF, int tdep, bool checktaxStr, int threads, size_t blockSize);
	static bool indexKey(const string& inF, int tdep, bool checktaxStr, IndexHeader& key);
	bool readIndex(const string& path, const IndexHeader& key);
	void writeIndex(const string& path, const IndexHeader& key) const;
	void mergeBlock(ParsedBlock& block);
	void buildIndex();
	bool betterLineage(const Lineage& a, const Lineage& b) const;

	string TaxFile;
	vector<string> tlevels;
	NamePool& pool;
	vector<uint32_t> nameIds;
	string subjectChars;
	vector<Lineage> lineages;                          // every record, duplicates included
	vector<uint64_t> slots;                            // hash tag << 32 | (lineage index + 1)
	size_t uniqueSubjects;
	size_t duplicateSubjects;
	size_t malformedLines, firstMalformedLine, extraColumnLines, firstExtraColumnLine;
};

// A reference lineage as seen through one hit: ranks beyond the hit's eligible depth, and a suppressed
// uncertain species, read as unknown. It refers to the reference taxonomy instead of copying it.
struct HitLineage
{
	HitLineage(const RefTax& reference, const RefTax::Lineage& refLineage) :
		ref(&reference), lineage(&refLineage), depth(refLineage.depth), speciesMasked(false) {}
	uint32_t nameId(int rank) const {
		if (rank >= depth || (speciesMasked && rank == 6)) { return 0; }
		return ref->nameId(*lineage, rank);
	}
	const string& get(int rank) const { return ref->name(nameId(rank)); }
	bool speciesUncertain() const { return lineage->speciesUncertain; }
	void limitByIdentity(double identity, const vector<double>& thresholds) {
		int eligible = 0;
		while (eligible < depth && eligible < (int)thresholds.size() && thresholds[eligible] <= identity) {
			++eligible;
		}
		depth = eligible;
	}

	const RefTax* ref;
	const RefTax::Lineage* lineage;
	int depth;
	bool speciesMasked;
};

// Boundaries of the tab-separated fields of one mapping line; a 13th field holds the rest of the line.
struct TabFields
{
	static constexpr size_t maxFields = 13;
	explicit TabFields(const string& line);
	std::array<std::pair<size_t, size_t>, maxFields> span;
	size_t count;
};

struct BlastRes
{
    BlastRes();
	BlastRes(const string&,int);
   bool parseFromLine(const string&, int);
	bool parseFields(const string&, const TabFields&, int);
	// Negative when this hit ranks above the other: identity, then query coverage, then alignment length.
	int compareRank(const BlastRes&) const;
	bool betterThan(const BlastRes&) const;
	static bool isColumnHeader(const string&);
	static int supportedColumnCount(const string&);
	static int supportedColumnCount(const TabFields&);
	static bool extractQueryToken(const string&, string&);

	string Query, Sbj;
	int alLen;
	double perID, eval,score;
	double Qcoverage;
	int queryCovered; // covered query bases, min(query span, alignment length)
	int queryLength;
	bool queryCoverageKnown; // false for legacy 12-column BLAST rows without qlen
	bool fail;
};

class BlastReader
{
public:
	BlastReader(const string&, const string&);
 ~BlastReader();
	// Replaces the batch with the next query's alignments; false when the input is exhausted.
	bool nextBatch(vector<BlastRes>& batch);
  vector<BlastRes> getResBatch();
	size_t queriesRead() const { return queryBlocks; }
	// Non-contiguous query blocks are detected among this many most recently completed queries.
	static const size_t maxTrackedQueries = 100000;
private:
	bool readRecord(BlastRes& result);
	void rememberCompletedQuery(const string&);
	void skipRecord(const char* problem);
	void finishFile();

	string fileName;
	vector<char> inputBuffer;
	bool processedBatch;
	    bool hasLastBlast;
	BlastRes lastBlast;
	BlastRes current;
	istream* blast;
	bool allRead;
	bool finished;
	bool legacyNoticeShown;
	int inptFmt;
	int detectedColumns;
	size_t queryBlocks;
	size_t lineNumber;
	size_t validRecords;
	size_t skippedRecords;
	size_t firstSkippedLine;
	  string lineBuffer;
	unordered_set<std::string_view> recentQueries; // views into recentOrder
	vector<string> recentOrder; // ring buffer of recently completed queries; never reallocated once filled
	size_t nextEviction;
};
