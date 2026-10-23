#include "RefTax.h"
#include "LCAlib.h"
#include "Pipeline.h"
#include <cerrno>
#include <cctype>
#include <charconv>
#include <climits>
#include <filesystem>
#include <limits>

namespace {
// Taxonomy naming rules, shared with the library.
using lca::iequals_ascii;
using lca::LineageViews;
using lca::parseLineage;

inline std::string_view field_view(const string& line, const std::pair<size_t, size_t>& span) {
	return std::string_view(line).substr(span.first, span.second - span.first);
}

inline std::string_view trim_number(std::string_view value) {
	size_t begin = 0, end = value.size();
	while (begin < end && std::isspace(static_cast<unsigned char>(value[begin]))) { ++begin; }
	while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1]))) { --end; }
	return value.substr(begin, end - begin);
}

// Like strtod/strtol on a whole field: surrounding whitespace and a leading '+' are accepted.
inline bool parse_double(std::string_view text, double& out) {
	text = trim_number(text);
	if (text.empty()) { return false; }
	if (text.front() == '+') {
		text.remove_prefix(1);
		if (text.empty() || text.front() == '+' || text.front() == '-') { return false; }
	}
#if defined(__cpp_lib_to_chars)
	const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), out);
	return error == std::errc() && end == text.data() + text.size() && std::isfinite(out);
#else
	const string copy(text);
	char* parseEnd = nullptr;
	errno = 0;
	out = std::strtod(copy.c_str(), &parseEnd);
	return parseEnd == copy.c_str() + copy.size() && errno != ERANGE && std::isfinite(out);
#endif
}

inline bool parse_int(std::string_view text, int& out) {
	text = trim_number(text);
	if (text.empty()) { return false; }
	if (text.front() == '+') {
		text.remove_prefix(1);
		if (text.empty() || text.front() == '+' || text.front() == '-') { return false; }
	}
	const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), out);
	return error == std::errc() && end == text.data() + text.size();
}

inline int known_tax_count(const TaxObj& tax) {
	int count = 0;
	for (int rank = 0; rank < tax.depth; ++rank) {
		if (tax.get(rank) != __unkwnTax) { ++count; }
	}
	return count;
}

// Lexicographic order of the reported ranks; names hidden beyond the depth do not count.
inline int compare_visible_lineage(const TaxObj& a, const TaxObj& b) {
	const int ranks = std::max(a.depth, b.depth);
	for (int rank = 0; rank < ranks; ++rank) {
		const int order = a.get(rank).compare(b.get(rank));
		if (order != 0) { return order; }
	}
	return 0;
}

// A block of complete taxonomy lines, and the line number of its first line.
struct TaxBlock {
	string text;
	size_t firstLine = 0;
};
}

//generic functions
void trim(string& str,
	const std::string& whitespace)
{
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos) {
		str.clear();
		return;
	}
	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	str = str.substr(strBegin, strRange);
}

bool isGZfile(const string fi) {
	if (fi.size() < 3) { return false; }
	string suffix = fi.substr(fi.size() - 3);
	std::transform(suffix.begin(), suffix.end(), suffix.begin(), [](unsigned char c) {
		return static_cast<char>(std::tolower(c));
	});
	return suffix == ".gz";
}



TaxObj::TaxObj(const TaxObj* t): SavedTaxs(t->SavedTaxs), Subj(t->Subj), hitDB(t->hitDB),
	perID(t->perID), repID(t->repID), hasHitDB(t->hasHitDB),
	speciesUncertain(t->speciesUncertain), depth(t->depth) {
}

string TaxObj::getWriteString(const vector<double>& ids) {
	int outDepth = (int)ids.size();
	if (outDepth <= 0) {
		outDepth = depth;
	}
	if (outDepth <= 0) {
		return "";
	}
	string ret;
	ret.reserve(128);
	for (int i = 0; i < outDepth; i++) {
		const string& cur = get(i);
		const bool unknown = i >= depth || cur == __unkwnTax ||
			(i < (int)ids.size() && perID < ids[i]);
		if (i > 0) { ret += __defaultTaxSep; }
		ret += unknown ? __unkwnTaxWR : cur;
	}
	if (hasHitDB) {
		ret += __defaultTaxSep;
		ret += hitDB.empty() ? __unkwnTaxWR : hitDB;
	}
	if (repID) {
		ret += __defaultTaxSep;
		ret += to_string(perID);
	}
	return ret;
}



TaxObj::TaxObj(std::string_view X,int d, bool /*nativeSLV*/, bool doNotCheckTax):
	SavedTaxs(), Subj(""), hitDB(""), perID(0.f), repID(false), hasHitDB(false),
	speciesUncertain(false), depth(0) {
	thread_local LineageViews lineage;
	parseLineage(X, d, doNotCheckTax, lineage);
	SavedTaxs.reserve(d);
	SavedTaxs.resize(lineage.depth);
	for (int rank = 0; rank < lineage.depth; ++rank) {
		SavedTaxs[rank].assign(lineage.ranks[rank].data(), lineage.ranks[rank].size());
	}
	depth = lineage.depth;
	speciesUncertain = lineage.speciesUncertain;
}

bool TaxObj::evalAcpyTax(const TaxObj* oth) {
	const int knownHere = known_tax_count(*this);
	const int knownOther = known_tax_count(*oth);
	const bool sameQuality = knownOther == knownHere && oth->depth == depth && oth->perID == perID;
	const int lineageOrder = sameQuality ? compare_visible_lineage(*oth, *this) : 0;
	const bool otherBetter = knownOther > knownHere ||
		(knownOther == knownHere && oth->depth > depth) ||
		(knownOther == knownHere && oth->depth == depth && oth->perID > perID) ||
		(sameQuality && oth->speciesUncertain != speciesUncertain && !oth->speciesUncertain) ||
		(sameQuality && oth->speciesUncertain == speciesUncertain && lineageOrder < 0) ||
		(sameQuality && oth->speciesUncertain == speciesUncertain && lineageOrder == 0 && oth->hitDB < hitDB);
	if (!otherBetter) { return false; }
	copyOver(oth);
	return true;
}

void TaxObj::copyOver(const TaxObj* oth) {
	SavedTaxs = oth->SavedTaxs;
	depth = oth->depth;
	perID = oth->perID;
	speciesUncertain = oth->speciesUncertain;
	hitDB = oth->hitDB;
	hasHitDB = oth->hasHitDB;
}


//*******************************************************
//        RefTax
//*******************************************************

// One parsed block: names are interned per block first, so the serial merge looks up each
// distinct name of a block once instead of every rank of every line.
struct RefTax::ParsedBlock {
	struct Record {
		std::string_view subject;
		uint32_t firstName; // into nameIds
		uint16_t depth;
		bool speciesUncertain;
	};
	string text;                          // owns the lines that the views refer to
	vector<std::string_view> names;       // block-local names; local ID = index + 1, 0 = unknown
	vector<uint32_t> nameIds;
	vector<Record> records;
	size_t malformed = 0, firstMalformed = 0, extraColumns = 0, firstExtraColumn = 0;
};

RefTax::RefTax(const string& inF, int tdep, bool /*nativeSLV*/, bool checktaxStr, int threads, size_t blockSize):
	TaxFile(inF), tlevels(tdep, ""), uniqueSubjects(0), duplicateSubjects(0),
	malformedLines(0), firstMalformedLine(0), extraColumnLines(0), firstExtraColumnLine(0)
{
	cout << "Loading tax DB.." << inF<<endl;
	names.push_back(&__unkwnTax);
	if (blockSize == 0) { blockSize = 1; }

	// Blocks are read in one call each, so the stream's own buffer is bypassed.
	ifstream in(inF.c_str(), ios::in | ios::binary);
	if (!in) { throw FatalError{11, "Cant open file " + inF + "\n"}; }

	string carry; // incomplete last line of the previous block
	size_t nextLine = 1;
	bool atEnd = false;
	auto produce = [&](TaxBlock& block) -> bool {
		block.text = std::move(carry);
		carry.clear();
		while (!atEnd) {
			const size_t old = block.text.size();
			block.text.resize(old + blockSize);
			in.read(&block.text[old], static_cast<std::streamsize>(blockSize));
			const size_t got = static_cast<size_t>(in.gcount());
			block.text.resize(old + got);
			if (got < blockSize) {
				if (in.bad()) { throw FatalError{13, "Failed while reading taxonomy file " + inF + "\n"}; }
				atEnd = true;
				break;
			}
			const size_t lastNewline = block.text.rfind('\n');
			if (lastNewline != string::npos) {
				carry.assign(block.text, lastNewline + 1, string::npos);
				block.text.resize(lastNewline + 1);
				break;
			}
		}
		if (block.text.empty()) { return false; }
		block.firstLine = nextLine;
		nextLine += static_cast<size_t>(std::count(block.text.begin(), block.text.end(), '\n'));
		return true;
	};

	const bool doNotCheckTax = !checktaxStr;
	auto work = [tdep, doNotCheckTax](TaxBlock& block, ParsedBlock& out) {
		out.text = std::move(block.text);
		thread_local LineageViews lineage;
		unordered_map<std::string_view, uint32_t> localIds;
		localIds.reserve(4096);
		const std::string_view text(out.text);
		size_t lineNumber = block.firstLine;
		for (size_t pos = 0; pos < text.size();) {
			size_t end = text.find('\n', pos);
			if (end == std::string_view::npos) { end = text.size(); }
			std::string_view line = text.substr(pos, end - pos);
			pos = end + 1;
			const size_t number = lineNumber++;
			if (!line.empty() && line.back() == '\r') { line.remove_suffix(1); }
			if (line.empty() || line[0] == '#') { continue; }
			const size_t tab = line.find('\t');
			if (tab == std::string_view::npos || tab == 0) {
				if (out.malformed++ == 0) { out.firstMalformed = number; }
				continue;
			}
			// Only the second column holds the lineage; later columns would otherwise leak into rank names.
			size_t lineageEnd = line.find('\t', tab + 1);
			if (lineageEnd == std::string_view::npos) {
				lineageEnd = line.size();
			} else if (line.find_first_not_of(" \t", lineageEnd) != std::string_view::npos) {
				if (out.extraColumns++ == 0) { out.firstExtraColumn = number; }
			}
			parseLineage(line.substr(tab + 1, lineageEnd - tab - 1), tdep, doNotCheckTax, lineage);
			ParsedBlock::Record record;
			record.subject = line.substr(0, tab);
			record.firstName = static_cast<uint32_t>(out.nameIds.size());
			record.depth = static_cast<uint16_t>(lineage.depth);
			record.speciesUncertain = lineage.speciesUncertain;
			for (int rank = 0; rank < lineage.depth; ++rank) {
				const std::string_view name = lineage.ranks[rank];
				uint32_t id = 0;
				if (!name.empty()) {
					const auto entry = localIds.try_emplace(name, static_cast<uint32_t>(out.names.size() + 1));
					if (entry.second) { out.names.push_back(name); }
					id = entry.first->second;
				}
				out.nameIds.push_back(id);
			}
			out.records.push_back(record);
		}
	};

	std::error_code sizeError;
	const uintmax_t fileSize = std::filesystem::file_size(inF, sizeError);
	bool reserved = sizeError || fileSize == 0;
	runOrdered<TaxBlock, ParsedBlock>(threads, produce, work, [&](ParsedBlock& block) {
		// Size the arrays once from the first block, instead of doubling them while they grow.
		if (!reserved && !block.text.empty()) {
			const double scale = 1.05 * static_cast<double>(fileSize) / static_cast<double>(block.text.size());
			lineages.reserve(static_cast<size_t>(scale * static_cast<double>(block.records.size())));
			nameIds.reserve(static_cast<size_t>(scale * static_cast<double>(block.nameIds.size())));
			size_t subjectBytes = 0;
			for (const auto& record : block.records) { subjectBytes += record.subject.size(); }
			subjectChars.reserve(static_cast<size_t>(scale * static_cast<double>(subjectBytes)));
			reserved = true;
		}
		mergeBlock(block);
	});
	buildIndex();
	// Release a poor size estimate; a close one is kept to avoid a second copy.
	if (lineages.capacity() > lineages.size() + lineages.size() / 4) { lineages.shrink_to_fit(); }
	if (nameIds.capacity() > nameIds.size() + nameIds.size() / 4) { nameIds.shrink_to_fit(); }
	if (subjectChars.capacity() > subjectChars.size() + subjectChars.size() / 4) { subjectChars.shrink_to_fit(); }

	if (malformedLines > 0) {
		if (uniqueSubjects == 0) {
			throw FatalError{12, "No valid taxonomy records in " + inF + "; line " + to_string(firstMalformedLine) +
				" lacks an identifier, a tab, and a lineage.\n"};
		}
		cerr << "Warning: skipped " << malformedLines << " taxonomy line(s) in " << inF
			<< " without an identifier, a tab, and a lineage (first at line " << firstMalformedLine << ").\n";
	}
	if (extraColumnLines > 0) {
		cerr << "Warning: ignored columns after the lineage on " << extraColumnLines << " taxonomy line(s) in "
			<< inF << " (first at line " << firstExtraColumnLine << ").\n";
	}
	cout << duplicateSubjects << " of " << duplicateSubjects + uniqueSubjects << " are duplicate entries\n";
	this->stats();
}

void RefTax::mergeBlock(ParsedBlock& block) {
	vector<uint32_t> globalIds(block.names.size() + 1, 0);
	for (size_t i = 0; i < block.names.size(); ++i) {
		auto found = nameIndex.find(block.names[i]);
		if (found == nameIndex.end()) {
			nameStore.emplace_back(block.names[i]);
			found = nameIndex.emplace(std::string_view(nameStore.back()), static_cast<uint32_t>(names.size())).first;
			names.push_back(&nameStore.back());
		}
		globalIds[i + 1] = found->second;
	}
	for (const auto& record : block.records) {
		if (nameIds.size() + record.depth > std::numeric_limits<uint32_t>::max() ||
			lineages.size() + 1 >= std::numeric_limits<uint32_t>::max()) {
			throw FatalError{14, "Taxonomy file " + TaxFile + " has too many entries.\n"};
		}
		Lineage lineage;
		lineage.subjectOffset = subjectChars.size();
		lineage.subjectLength = static_cast<uint32_t>(record.subject.size());
		subjectChars.append(record.subject);
		lineage.firstName = static_cast<uint32_t>(nameIds.size());
		lineage.depth = record.depth;
		lineage.speciesUncertain = record.speciesUncertain;
		for (uint16_t rank = 0; rank < record.depth; ++rank) {
			nameIds.push_back(globalIds[block.nameIds[record.firstName + rank]]);
		}
		lineages.push_back(lineage);
	}
	if (block.malformed > 0) {
		if (malformedLines == 0) { firstMalformedLine = block.firstMalformed; }
		malformedLines += block.malformed;
	}
	if (block.extraColumns > 0) {
		if (extraColumnLines == 0) { firstExtraColumnLine = block.firstExtraColumn; }
		extraColumnLines += block.extraColumns;
	}
}

// Duplicate identifiers keep the same lineage regardless of record order: more known ranks, then
// more ranks, then a certain species, then the lexicographically smaller lineage.
bool RefTax::betterLineage(const Lineage& a, const Lineage& b) const {
	auto known = [this](const Lineage& lineage) {
		int count = 0;
		for (int rank = 0; rank < lineage.depth; ++rank) { count += nameId(lineage, rank) != 0; }
		return count;
	};
	const int knownA = known(a), knownB = known(b);
	if (knownA != knownB) { return knownA > knownB; }
	if (a.depth != b.depth) { return a.depth > b.depth; }
	if (a.speciesUncertain != b.speciesUncertain) { return !a.speciesUncertain; }
	for (int rank = 0; rank < a.depth; ++rank) {
		const int order = rankName(a, rank).compare(rankName(b, rank));
		if (order != 0) { return order < 0; }
	}
	return false;
}

void RefTax::buildIndex() {
	size_t capacity = 16;
	while (capacity < lineages.size() * 2) { capacity <<= 1; }
	slots.assign(capacity, 0);
	const size_t mask = capacity - 1;
	for (size_t i = 0; i < lineages.size(); ++i) {
		const uint64_t hash = std::hash<std::string_view>{}(subject(lineages[i]));
		const uint64_t tag = hash >> 32;
		for (size_t pos = hash & mask;; pos = (pos + 1) & mask) {
			uint64_t& slot = slots[pos];
			if (slot == 0) {
				slot = (tag << 32) | (i + 1);
				uniqueSubjects++;
				break;
			}
			if ((slot >> 32) == tag) {
				const size_t other = (slot & 0xffffffffu) - 1;
				if (subject(lineages[other]) == subject(lineages[i])) {
					duplicateSubjects++;
					if (betterLineage(lineages[i], lineages[other])) { slot = (tag << 32) | (i + 1); }
					break;
				}
			}
		}
	}
}

const RefTax::Lineage* RefTax::find(std::string_view s) const {
	if (slots.empty()) { return nullptr; }
	const uint64_t hash = std::hash<std::string_view>{}(s);
	const uint64_t tag = hash >> 32;
	const size_t mask = slots.size() - 1;
	for (size_t pos = hash & mask;; pos = (pos + 1) & mask) {
		const uint64_t slot = slots[pos];
		if (slot == 0) { return nullptr; }
		if ((slot >> 32) == tag) {
			const Lineage& lineage = lineages[(slot & 0xffffffffu) - 1];
			if (subject(lineage) == s) { return &lineage; }
		}
	}
}

void RefTax::stats() {
	int cnt = 0;
	vector<int> hist(tlevels.size() + 1, 0);
	int maxD = 0;
	for (const uint64_t slot : slots) {
		if (slot == 0) { continue; }
		const Lineage& lineage = lineages[(slot & 0xffffffffu) - 1];
		const int dep = lineage.depth;
		if (dep >= (int)hist.size()) {
			cerr << "Tax depth " << dep << " of object " << subject(lineage) << " is outside the configured range\n";
			continue;
		}
		if (maxD < dep) { maxD = dep; }
		hist[dep]++;
		cnt++;
	}
	cout << "TaxDB " << TaxFile << " contained " << cnt << " entries, depth distribution is:\n";
	for (int i = 0; i < maxD+1; i++) {
		cout << i << ":" << hist[i] << " ";
	}
	cout << endl;
}


//*******************************************************
//        BlastRes
//*******************************************************

TabFields::TabFields(const string& line) : span(), count(0) {
	const char* data = line.data();
	const size_t size = line.size();
	size_t begin = 0;
	while (true) {
		const void* tab = count + 1 < maxFields ? std::memchr(data + begin, '\t', size - begin) : nullptr;
		const size_t end = tab != nullptr ? static_cast<size_t>(static_cast<const char*>(tab) - data) : size;
		span[count++] = {begin, end};
		if (tab == nullptr) { break; }
		begin = end + 1;
	}
}

int BlastRes::compareRank(const BlastRes& other) const {
	if (perID != other.perID) { return perID > other.perID ? -1 : 1; }
	if (queryCoverageKnown != other.queryCoverageKnown) { return queryCoverageKnown ? -1 : 1; }
	if (queryCoverageKnown && Qcoverage != other.Qcoverage) { return Qcoverage > other.Qcoverage ? -1 : 1; }
	if (alLen != other.alLen) { return alLen > other.alLen ? -1 : 1; }
	return 0;
}

bool BlastRes::betterThan(const BlastRes& other) const {
	const int rank = compareRank(other);
	return rank != 0 ? rank < 0 : Sbj < other.Sbj;
}

BlastRes::BlastRes() :
	Query(""), Sbj(""), alLen(0), perID(0.f), eval(-1.f), score(0.f),
	Qcoverage(0.0), queryCovered(0), queryLength(0), queryCoverageKnown(false), fail(true) {
}

BlastRes::BlastRes(const string& line, int inptFmt):
		Query(""), Sbj(""), alLen(0), perID(0.f), eval(-1.f), score(0.f),
		Qcoverage(0.0), queryCovered(0), queryLength(0), queryCoverageKnown(false), fail(true) {
	parseFromLine(line, inptFmt);
}

bool BlastRes::extractQueryToken(const string& line, string& query) {
	const TabFields fields(line);
	query.assign(field_view(line, fields.span[0]));
	return query.find_first_not_of(" \t\r\n\v\f") != string::npos;
}

int BlastRes::supportedColumnCount(const TabFields& fields) {
	return fields.count == 11 || fields.count == 12 ? static_cast<int>(fields.count) : 0;
}

int BlastRes::supportedColumnCount(const string& line) {
	return supportedColumnCount(TabFields(line));
}

bool BlastRes::isColumnHeader(const string& line) {
	static const char* common[] = {
		"qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
		"qstart", "qend", "sstart", "send"
	};
	const TabFields fields(line);
	if (fields.count != 11 && fields.count != 12) { return false; }
	auto label = [&](size_t i) { return trim_number(field_view(line, fields.span[i])); };
	for (size_t i = 0; i < 10; ++i) {
		if (!iequals_ascii(label(i), common[i])) { return false; }
	}
	if (fields.count == 11) {
		return iequals_ascii(label(10), "qlen");
	}
	return iequals_ascii(label(10), "evalue") && iequals_ascii(label(11), "bitscore");
}

bool BlastRes::parseFromLine(const string& line, int inptFmt) {
	return parseFields(line, TabFields(line), inptFmt);
}

bool BlastRes::parseFields(const string& line, const TabFields& fields, int inptFmt) {
	fail = true;
	Query.clear();
	Sbj.clear();
	alLen = 0;
	perID = 0.0;
	eval = -1.0;
	score = 0.0;
	Qcoverage = 0.0;
	queryCovered = 0;
	queryLength = 0;
	queryCoverageKnown = false;

	if (inptFmt != 0 || line.empty()) { return false; }
	if (fields.count != 11 && fields.count != 12) { return false; }
	auto field = [&](size_t i) { return field_view(line, fields.span[i]); };
	Query.assign(field(0));
	Sbj.assign(field(1));
	if (Query.find_first_not_of(" \r\n\v\f") == string::npos ||
		Sbj.find_first_not_of(" \r\n\v\f") == string::npos) { return false; }
	if (!parse_double(field(2), perID)) { return false; }
	if (!parse_int(field(3), alLen)) { return false; }
	int mismatches = 0, gaps = 0, qstart = 0, qstop = 0;
	int sstart = 0, sstop = 0;
	if (!parse_int(field(4), mismatches) || !parse_int(field(5), gaps) || !parse_int(field(6), qstart) ||
		!parse_int(field(7), qstop) || !parse_int(field(8), sstart) || !parse_int(field(9), sstop)) {
		return false;
	}
	if (perID < 0.0 || perID > 100.0 || alLen <= 0 || mismatches < 0 || gaps < 0 ||
		qstart <= 0 || qstop <= 0 || sstart <= 0 || sstop <= 0) {
		return false;
	}

	if (fields.count == 11) {
		int qlen = 0;
		if (!parse_int(field(10), qlen) || qlen <= 0 || qstart > qlen || qstop > qlen) {
			return false;
		}
		const long long querySpan = std::llabs(static_cast<long long>(qstop) - qstart) + 1;
		queryCovered = static_cast<int>(std::min<long long>(querySpan, alLen));
		queryLength = qlen;
		Qcoverage = std::min(1.0, static_cast<double>(queryCovered) / qlen);
		queryCoverageKnown = true;
	} else {
		if (!parse_double(field(10), eval) || !parse_double(field(11), score) || eval < 0.0 || score < 0.0) {
			return false;
		}
	}
	fail = false;
	return true;
}


//*******************************************************
//        BlastReader
//*******************************************************


BlastReader::BlastReader(const string& inf, const string& inFmt): fileName(inf), inputBuffer(), processedBatch(false),
	hasLastBlast(false), blast(NULL), allRead(false), finished(false), legacyNoticeShown(false),
	inptFmt(-1), detectedColumns(0), queryBlocks(0), lineNumber(0), validRecords(0), skippedRecords(0),
	firstSkippedLine(0), lineBuffer(), recentQueries(), recentOrder(), nextEviction(0) {
	if (inFmt != "bl8") {
		throw FatalError{26, "Unsupported input format '" + inFmt + "'. Only the custom bl8 format is supported.\n"};
	}
	inptFmt = 0;
	if (isGZfile(inf)) {
#ifdef _gzipread
		blast = new igzstream(inf.c_str(), ios::in);
#else
		throw FatalError{50, "gzip not supported in your LCA build\n"};
#endif
	} else {
		inputBuffer.resize(__inputBufferSize);
		ifstream* file = new ifstream();
		file->rdbuf()->pubsetbuf(inputBuffer.data(), inputBuffer.size());
		file->open(inf.c_str(), ios::in | ios::binary);
		blast = file;
	}
	if (!*blast) {
		delete blast;
		blast = NULL;
		throw FatalError{23, "Blast input file " + inf + " could not be opened. exiting..\n"};
	}
	lineBuffer.reserve(512);
}

BlastReader::~BlastReader() {
	if (blast != NULL) {
		delete blast;
		blast = NULL;
	}
}

void BlastReader::rememberCompletedQuery(const string& query) {
	if (recentOrder.size() < maxTrackedQueries) {
		// The ring is reserved once, so the views in recentQueries stay valid.
		if (recentOrder.capacity() < maxTrackedQueries) { recentOrder.reserve(maxTrackedQueries); }
		recentOrder.push_back(query);
		recentQueries.insert(recentOrder.back());
		return;
	}
	string& slot = recentOrder[nextEviction];
	recentQueries.erase(slot);
	slot = query;
	recentQueries.insert(slot);
	nextEviction = (nextEviction + 1) % maxTrackedQueries;
}

void BlastReader::skipRecord(const char* problem) {
	static const size_t maxListedRecords = 5;
	if (skippedRecords++ == 0) { firstSkippedLine = lineNumber; }
	if (skippedRecords <= maxListedRecords) {
		cerr << "Warning: skipping m8 record at line " << lineNumber << " of " << fileName << ": " << problem << ".\n";
	}
}

void BlastReader::finishFile() {
	if (finished) { return; }
	finished = true;
	if (blast->bad() || (blast->fail() && !blast->eof())) {
		throw FatalError{28, "Failed while reading blast input after line " + to_string(lineNumber) + "\n"};
	}
	if (skippedRecords == 0) { return; }
	if (validRecords == 0) {
		throw FatalError{25, "No valid m8 records in " + fileName + "; first malformed record at line " +
			to_string(firstSkippedLine) + ". Expected tab-separated fields: either 11 columns ending in qlen, or the "
			"legacy 12-column BLAST layout ending in evalue and bitscore.\n"};
	}
	cerr << "Warning: skipped " << skippedRecords << " malformed m8 record(s) in " << fileName
		<< " (first at line " << firstSkippedLine << ").\n";
}

bool BlastReader::readRecord(BlastRes& result) {
	while (getline(*blast, lineBuffer, '\n')) {
		lineNumber++;
		if (!lineBuffer.empty() && lineBuffer.back() == '\r') { lineBuffer.pop_back(); }
		const size_t contentStart = lineBuffer.find_first_not_of(" \t\r");
		if (contentStart == string::npos || lineBuffer[contentStart] == '#') { continue; }
		const TabFields fields(lineBuffer);
		const int columns = BlastRes::supportedColumnCount(fields);
		const char* problem = nullptr;
		if (columns == 0) {
			problem = "expected 11 or 12 tab-separated columns";
		} else if (detectedColumns != 0 && columns != detectedColumns) {
			problem = detectedColumns == 11 ? "the file uses the 11-column layout, but this row has 12 columns"
				: "the file uses the 12-column layout, but this row has 11 columns";
		} else if (!result.parseFields(lineBuffer, fields, inptFmt)) {
			problem = "a field is empty, not a number, or out of range";
		}
		if (problem == nullptr) {
			if (detectedColumns == 0) {
				detectedColumns = columns;
				if (columns == 12 && !legacyNoticeShown) {
					cerr << "Notice: legacy 12-column BLAST input has no qlen; query-coverage filtering "
						<< "is skipped for these records.\n";
					legacyNoticeShown = true;
				}
			}
			validRecords++;
			return true;
		}
		// Column headers may repeat, e.g. in concatenated files.
		if (BlastRes::isColumnHeader(lineBuffer)) { continue; }
		skipRecord(problem);
	}
	finishFile();
	return false;
}

vector<BlastRes> BlastReader::getResBatch() {
	vector<BlastRes> batch;
	nextBatch(batch);
	return batch;
}

bool BlastReader::nextBatch(vector<BlastRes>& batch) {
	batch.clear();
	if (!processedBatch) {
		hasLastBlast = readRecord(lastBlast);
		processedBatch = true;
		if (!hasLastBlast) {
			allRead = true;
			return false;
		}
	}

	if (!hasLastBlast || allRead) {
		return false;
	}

	batch.push_back(std::move(lastBlast));
	queryBlocks++;

	while (readRecord(current)) {
		if (current.Query != batch.front().Query) {
			rememberCompletedQuery(batch.front().Query);
			if (recentQueries.find(current.Query) != recentQueries.end()) {
				throw FatalError{27, "Query '" + current.Query + "' occurs in multiple non-contiguous blocks at line " +
					to_string(lineNumber) + ". Sort/group the custom m8 input by query.\n"};
			}
			std::swap(lastBlast, current);
			hasLastBlast = true;
			return true;
		}

		// Eligibility depends on command-line filters, so deduplicate later.
		batch.push_back(current);
	}

	hasLastBlast = false;
	allRead = true;
	return true;
}
