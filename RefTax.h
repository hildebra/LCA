#pragma once
#include "libload.h"
#include "options.h"
void trim(string& str,const std::string& whitespace = " \t\r\n");
bool isGZfile(const string fi);


struct TaxObj
{
	TaxObj(const string&, int, bool nativeSLV, bool doNotCheckTax);
	TaxObj(const TaxObj* t);
	TaxObj(int d) : SavedTaxs(), Subj(""), hitDB(""), perID(0.f), repID(false),
		hasHitDB(false), speciesUncertain(false), depth(0) { SavedTaxs.reserve(d); }
	//functions
	string getWriteString(const vector<double>&);
	void copy_vals(const TaxObj* t) {
		SavedTaxs = t->SavedTaxs;
		depth = t->depth;
		speciesUncertain = t->speciesUncertain;
	}
	void setRepID(bool x) { repID = x; }
	void makeSpeciesUnknown() {
		const int speciesRank = 6;
		if (speciesUncertain && depth > speciesRank) {
			if ((int)SavedTaxs.size() <= speciesRank) { SavedTaxs.resize(speciesRank + 1, __unkwnTax); }
			SavedTaxs[speciesRank] = __unkwnTax;
		}
	}
	//get tax at depth x
	const string& get(int x) const {
		if (x < 0 || x >= depth || x >= (int)SavedTaxs.size()) { return __unkwnTax; }
		return SavedTaxs[x];
	}
	void set(int x, string v) { 
       if (x < 0) { return; }
		if (x >= (int)SavedTaxs.size()) { SavedTaxs.resize(x + 1, __unkwnTax); }
		SavedTaxs[x] = v;
		if (x >= depth) { depth = x + 1; }
	}
	//check if other tax is better and copies if so these vals over itself
	bool evalAcpyTax(const TaxObj* oth);
	void copyOver(const TaxObj* oth);
	void setHitDB(const string& x) { hitDB = x; hasHitDB = true; }
	const string& getHitDB() const { return hitDB; }
	bool reportsHitDB() const { return hasHitDB; }

	//int dept() { return depth; }
	//variables
	vector<string> SavedTaxs;
	string Subj;
	string hitDB;
	float perID;
	bool repID;
	bool hasHitDB;
	bool speciesUncertain;
	int depth;//number of taxonomic levels assigned/available; not the configured maximum

};


class RefTax
{
public:
	RefTax(const string&,int tdep,bool,bool);
	~RefTax();
	void stats();
	int depth() const { return (int) tlevels.size(); }
	unordered_map <string, TaxObj*>::const_iterator find(const string& s) const {
		return Tlink.find(s);
	}
	unordered_map <string, TaxObj*>::const_iterator end() const {
		return Tlink.end();
	}
	void setTaxLvls(vector<string> x) { tlevels = x; }
	//const string getLevels();
	//vector<string> taxLevls() { return tlevels; }
private:
	string TaxFile;
	unordered_map <string, TaxObj*> Tlink;
	vector<string> tlevels;
};

struct BlastRes
{
    BlastRes();
	BlastRes(const string&,int);
   bool parseFromLine(const string&, int);
	static bool isColumnHeader(const string&);
	static int supportedColumnCount(const string&);
	static bool extractQueryToken(const string&, string&);
	bool isSameQuery(const string &q) const {	if (q == Query) { return true; } return false;	}
	
	string Query, Sbj;
	int alLen;
	double perID, eval,score;
	float Qcoverage;
	bool queryCoverageKnown; // false for legacy 12-column BLAST rows without qlen
	bool fail;
};

class BlastReader
{
public:
	BlastReader(const string&, const string&);
 ~BlastReader();
  vector<BlastRes> getResBatch();
private:

	bool processedBatch;
	    bool hasLastBlast;
	BlastRes lastBlast;
	istream* blast;
	bool allRead;
	bool seenData;
	bool legacyNoticeShown;
	int inptFmt;
	int detectedColumns;
	int blastCnter;
	size_t lineNumber;
	  string lineBuffer;
	unordered_map<string, size_t> foundSbjs;
	unordered_set<string> completedQueries;
	vector<BlastRes> batchBuffer;
};
