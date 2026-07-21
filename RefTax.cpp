#include "RefTax.h"
#include <cerrno>
#include <cctype>
#include <climits>
#include <limits>

namespace {
inline bool iequals_ascii_range(const string& src, size_t begin, size_t end, const char* txt) {
	const size_t len = end - begin;
	size_t i = 0;
	for (; txt[i] != '\0'; ++i) {
		if (i >= len) { return false; }
		if (std::tolower(static_cast<unsigned char>(src[begin + i])) != std::tolower(static_cast<unsigned char>(txt[i]))) {
			return false;
		}
	}
	return i == len;
}

inline bool istarts_with_ascii_range(const string& src, size_t begin, size_t end, const char* txt) {
	const size_t len = end - begin;
	for (size_t i = 0; txt[i] != '\0'; ++i) {
		if (i >= len) { return false; }
		if (std::tolower(static_cast<unsigned char>(src[begin + i])) != std::tolower(static_cast<unsigned char>(txt[i]))) {
			return false;
		}
	}
	return true;
}

inline bool parse_double_range(const string& src, size_t begin, size_t end, double& out) {
	if (begin >= end) { return false; }
	const char* first = src.data() + begin;
	const char* last = src.data() + end;
	char* parseEnd = nullptr;
	errno = 0;
	out = std::strtod(first, &parseEnd);
	return parseEnd == last && errno != ERANGE && std::isfinite(out);
}

inline bool parse_int_range(const string& src, size_t begin, size_t end, int& out) {
	if (begin >= end) { return false; }
	const char* first = src.data() + begin;
	const char* last = src.data() + end;
	char* parseEnd = nullptr;
	errno = 0;
	long parsed = std::strtol(first, &parseEnd, 10);
	if (parseEnd != last || errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX) { return false; }
	out = static_cast<int>(parsed);
	return true;
}

inline vector<pair<size_t, size_t> > whitespace_fields(const string& line) {
	vector<pair<size_t, size_t> > fields;
	fields.reserve(12);
	size_t pos = 0;
	while (pos < line.size()) {
		while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos]))) { ++pos; }
		if (pos == line.size()) { break; }
		const size_t begin = pos;
		while (pos < line.size() && !std::isspace(static_cast<unsigned char>(line[pos]))) { ++pos; }
		fields.emplace_back(begin, pos);
	}
	return fields;
}

inline int tax_rank_from_prefix(char rank) {
	switch (std::tolower(static_cast<unsigned char>(rank))) {
	case 'd':
	case 'k': return 0;
	case 'p': return 1;
	case 'c': return 2;
	case 'o': return 3;
	case 'f': return 4;
	case 'g': return 5;
	case 's': return 6;
	case 't': return 7;
	default: return -1;
	}
}

inline int known_tax_count(const TaxObj& tax) {
	return static_cast<int>(std::count_if(tax.SavedTaxs.begin(), tax.SavedTaxs.end(),
		[](const string& value) { return value != __unkwnTax; }));
}

inline bool contains_species_sp_marker(const string& value) {
	for (size_t pos = 0; pos + 2 < value.size(); ++pos) {
		if (std::tolower(static_cast<unsigned char>(value[pos])) != 's' ||
			std::tolower(static_cast<unsigned char>(value[pos + 1])) != 'p' || value[pos + 2] != '.') {
			continue;
		}
		const bool leftBoundary = pos == 0 || std::isspace(static_cast<unsigned char>(value[pos - 1]));
		const bool rightBoundary = pos + 3 == value.size() || std::isspace(static_cast<unsigned char>(value[pos + 3]));
		if (leftBoundary && rightBoundary) { return true; }
	}
	return false;
}

inline bool better_blast_hit(const BlastRes& lhs, const BlastRes& rhs) {
	if (lhs.perID != rhs.perID) { return lhs.perID > rhs.perID; }
	if (lhs.alLen != rhs.alLen) { return lhs.alLen > rhs.alLen; }
	if (lhs.queryCoverageKnown != rhs.queryCoverageKnown) { return lhs.queryCoverageKnown; }
	if (lhs.queryCoverageKnown && lhs.Qcoverage != rhs.Qcoverage) { return lhs.Qcoverage > rhs.Qcoverage; }
	return lhs.Sbj < rhs.Sbj;
}
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
	for (int i = 0; i < outDepth; i++) {
		const string& cur = get(i);
		const bool unknown = i >= depth || cur == __unkwnTax ||
			(i < (int)ids.size() && perID < ids[i]);
		if (i > 0) { ret += __defaultTaxSep; }
		ret += unknown ? __unkwnTaxWR : cur;
	}
	if (hasHitDB) {
		ret += __defaultTaxSep + (hitDB.empty() ? __unkwnTaxWR : hitDB);
	}
	if (repID) {
		ret += __defaultTaxSep + to_string(perID);
	}
	return ret;
}



TaxObj::TaxObj(const string& X,int d, bool /*nativeSLV*/, bool doNotCheckTax):
	SavedTaxs(), Subj(""), hitDB(""), perID(0.f), repID(false), hasHitDB(false),
	speciesUncertain(false), depth(0) {
	SavedTaxs.reserve(d);
	int sequentialRank = 0;
	size_t tokenBegin = 0;
	while (tokenBegin <= X.size()) {
		const size_t separator = X.find(';', tokenBegin);
		const size_t tokenEnd = separator == string::npos ? X.size() : separator;
		string token = X.substr(tokenBegin, tokenEnd - tokenBegin);
		trim(token, " \t\r\n");

		int rank = sequentialRank;
		if (token.size() >= 3 && token[1] == '_' && token[2] == '_') {
			const int prefixedRank = tax_rank_from_prefix(token[0]);
			if (prefixedRank >= 0) {
				rank = prefixedRank;
				token.erase(0, 3);
				trim(token, " \t\r\n");
			}
		}
		sequentialRank = std::max(sequentialRank + 1, rank + 1);

		if (rank >= 0 && rank < d) {
			const bool taxKnown = !token.empty() &&
				!iequals_ascii_range(token, 0, token.size(), "unclassified") &&
				!iequals_ascii_range(token, 0, token.size(), "uncultured bacterium") &&
				!iequals_ascii_range(token, 0, token.size(), "uncultured") && token != "?";

			if (rank == 6 && !token.empty()) {
				const bool currentSpeciesUncertain =
					istarts_with_ascii_range(token, 0, token.size(), "uncultured") ||
					istarts_with_ascii_range(token, 0, token.size(), "unclassified") ||
					contains_species_sp_marker(token);
				speciesUncertain = currentSpeciesUncertain;
				if (taxKnown && !currentSpeciesUncertain) {
					const size_t firstSpace = token.find(' ');
					if (firstSpace != string::npos) {
						size_t truncateAt = token.find(' ', firstSpace + 1);
						if (istarts_with_ascii_range(token, 0, token.size(), "candidatus ") && truncateAt != string::npos) {
							truncateAt = token.find(' ', truncateAt + 1);
						}
						if (truncateAt != string::npos) { token.resize(truncateAt); }
					}
				}
			}

			if (doNotCheckTax || taxKnown) {
				if ((int)SavedTaxs.size() <= rank) { SavedTaxs.resize(rank + 1, __unkwnTax); }
				SavedTaxs[rank] = token;
				depth = std::max(depth, rank + 1);
			}
		}

		if (separator == string::npos) { break; }
		tokenBegin = separator + 1;
	}
	if ((int)SavedTaxs.size() > depth) {
		SavedTaxs.resize(depth);
	}
}

bool TaxObj::evalAcpyTax(const TaxObj* oth) {
	const int knownHere = known_tax_count(*this);
	const int knownOther = known_tax_count(*oth);
	const bool sameQuality = knownOther == knownHere && oth->depth == depth && oth->perID == perID;
	const bool otherBetter = knownOther > knownHere ||
		(knownOther == knownHere && oth->depth > depth) ||
		(knownOther == knownHere && oth->depth == depth && oth->perID > perID) ||
		(sameQuality && oth->speciesUncertain != speciesUncertain && !oth->speciesUncertain) ||
		(sameQuality && oth->speciesUncertain == speciesUncertain && oth->SavedTaxs < SavedTaxs) ||
		(sameQuality && oth->speciesUncertain == speciesUncertain && oth->SavedTaxs == SavedTaxs &&
			oth->hitDB < hitDB);
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



RefTax::RefTax(const string& inF, int tdep, bool nativeSLV,bool checktaxStr):TaxFile(inF),
tlevels(tdep,"")
{
	//ini constants
	//#my @taxLvls = ("domain", "phylum", "class", "order", "family", "genus");
	//cout << "DEBUG\n\n"; return;
	cout << "Loading tax DB.." << inF<<endl;

	string line;
	ifstream in(inF.c_str());
	if (!in) { cerr << "Cant open file " << inF << endl; std::exit(11); }
	size_t TaxDbl(0), TaxSingl(0), taxLineNumber(0);
	while (getline(in, line, '\n')) {
		taxLineNumber++;
		if (!line.empty() && line.back() == '\r') { line.pop_back(); }
		if (line.empty() || line[0] == '#') { continue; }
		size_t dlmt = line.find("\t");
		if (dlmt == std::string::npos) {
			cerr << "Malformed taxonomy record at line " << taxLineNumber
				<< ": expected an identifier, a tab, and taxonomy.\n";
			exit(12);
		}
		string ID = line.substr(0, dlmt);
		if (ID.empty()) {
			cerr << "Malformed taxonomy record at line " << taxLineNumber
				<< ": the identifier is empty.\n";
			exit(12);
		}
		//string ttax = line.substr(dlmt+1);
		TaxObj* t = new TaxObj(line.substr(dlmt + 1), tdep, nativeSLV, !checktaxStr);
		auto fnd = Tlink.find(ID);
		if (fnd == Tlink.end()){//all good
			Tlink[ID] = t; TaxSingl++;
		} else {//not good: tax is double annotated
			// Select the same annotation regardless of duplicate record order.
			fnd->second->evalAcpyTax(t);
           delete t;
			TaxDbl++;
		}
	}
	//cerr << "C1\n";
	cout << TaxDbl << " of " << TaxDbl + TaxSingl << " are duplicate entries\n";
	//cerr << "C2\n";
	this->stats();
	//cerr << "C3\n";
}

RefTax::~RefTax()
{
	for (auto it = Tlink.begin(); it != Tlink.end(); ++it){
		//std::cout << " " << it->first << ":" << it->second;
		delete it->second;
	}
}
void RefTax::stats() {
	int cnt = 0;
	vector<int> hist(tlevels.size() + 1, 0);
	int maxD = 0;
	for (auto it = Tlink.begin(); it != Tlink.end(); ++it) {
		int dep = it->second->depth;
		if (dep < 0 || dep >= (int)hist.size()) {
			cerr << "Tax depth " << dep << " of object " << it->first << " is outside the configured range\n";
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

BlastRes::BlastRes() :
	Query(""), Sbj(""), alLen(0), perID(0.f), eval(-1.f), score(0.f),
	Qcoverage(0.f), queryCoverageKnown(false), fail(true) {
}

BlastRes::BlastRes(const string& line, int inptFmt):
		Query(""), Sbj(""), alLen(0), perID(0.f), eval(-1.f), score(0.f),
		Qcoverage(0.f), queryCoverageKnown(false), fail(true) {
	parseFromLine(line, inptFmt);
}

bool BlastRes::extractQueryToken(const string& line, string& query) {
	const vector<pair<size_t, size_t> > fields = whitespace_fields(line);
	if (fields.empty()) { return false; }
	query.assign(line, fields[0].first, fields[0].second - fields[0].first);
	return true;
}

int BlastRes::supportedColumnCount(const string& line) {
	int count = 0;
	size_t pos = 0;
	while (pos < line.size()) {
		while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos]))) { ++pos; }
		if (pos == line.size()) { break; }
		++count;
		if (count > 12) { return 0; }
		while (pos < line.size() && !std::isspace(static_cast<unsigned char>(line[pos]))) { ++pos; }
	}
	return count == 11 || count == 12 ? count : 0;
}

bool BlastRes::isColumnHeader(const string& line) {
	static const char* common[] = {
		"qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
		"qstart", "qend", "sstart", "send"
	};
	const vector<pair<size_t, size_t> > fields = whitespace_fields(line);
	if (fields.size() != 11 && fields.size() != 12) { return false; }
	for (size_t i = 0; i < 10; ++i) {
		if (!iequals_ascii_range(line, fields[i].first, fields[i].second, common[i])) { return false; }
	}
	if (fields.size() == 11) {
		return iequals_ascii_range(line, fields[10].first, fields[10].second, "qlen");
	}
	return iequals_ascii_range(line, fields[10].first, fields[10].second, "evalue") &&
		iequals_ascii_range(line, fields[11].first, fields[11].second, "bitscore");
}

bool BlastRes::parseFromLine(const string& line, int inptFmt) {
	fail = true;
	Query.clear();
	Sbj.clear();
	alLen = 0;
	perID = 0.0;
	eval = -1.0;
	score = 0.0;
	Qcoverage = 0.f;
	queryCoverageKnown = false;

	if (inptFmt != 0 || line.empty()) { return false; }

	const vector<pair<size_t, size_t> > fields = whitespace_fields(line);
	if (fields.size() != 11 && fields.size() != 12) { return false; }
	Query.assign(line, fields[0].first, fields[0].second - fields[0].first);
	Sbj.assign(line, fields[1].first, fields[1].second - fields[1].first);
	if (Query.empty() || Sbj.empty()) { return false; }
	if (!parse_double_range(line, fields[2].first, fields[2].second, perID)) { return false; }
	if (!parse_int_range(line, fields[3].first, fields[3].second, alLen)) { return false; }
	int mismatches = 0, gaps = 0, qstart = 0, qstop = 0;
	int sstart = 0, sstop = 0;
	if (!parse_int_range(line, fields[4].first, fields[4].second, mismatches)) { return false; }
	if (!parse_int_range(line, fields[5].first, fields[5].second, gaps)) { return false; }
	if (!parse_int_range(line, fields[6].first, fields[6].second, qstart)) { return false; }
	if (!parse_int_range(line, fields[7].first, fields[7].second, qstop)) { return false; }
	if (!parse_int_range(line, fields[8].first, fields[8].second, sstart)) { return false; }
	if (!parse_int_range(line, fields[9].first, fields[9].second, sstop)) { return false; }
	if (perID < 0.0 || perID > 100.0 || alLen <= 0 || mismatches < 0 || gaps < 0 ||
		qstart < 0 || qstop < 0 || sstart < 0 || sstop < 0) {
		return false;
	}

	if (fields.size() == 11) {
		int qlen = 0;
		if (!parse_int_range(line, fields[10].first, fields[10].second, qlen) ||
			qlen <= 0 || qstart > qlen || qstop > qlen) {
			return false;
		}
		const long long querySpan = std::llabs(static_cast<long long>(qstop) - qstart) + 1;
		const long long coveredQueryBases = std::min<long long>(querySpan, alLen);
		Qcoverage = static_cast<float>(std::min(1.0, static_cast<double>(coveredQueryBases) / qlen));
		queryCoverageKnown = true;
	} else {
		if (!parse_double_range(line, fields[10].first, fields[10].second, eval) ||
			!parse_double_range(line, fields[11].first, fields[11].second, score) ||
			eval < 0.0 || score < 0.0) {
			return false;
		}
	}
	fail = false;
	return true;
}


//*******************************************************
//        BlastReader
//*******************************************************


BlastReader::BlastReader(const string& inf, const string& inFmt): processedBatch(false),
	hasLastBlast(false), blast(NULL), allRead(false), seenData(false), legacyNoticeShown(false),
	inptFmt(-1), detectedColumns(0), blastCnter(0), lineNumber(0),
	lineBuffer(), foundSbjs(), completedQueries(), batchBuffer() {
#ifdef DEBUG
	cerr << "ini blast file\n";
#endif // DEBUG

	if (isGZfile(inf)) {
#ifdef _gzipread
		blast = new igzstream(inf.c_str(), ios::in);
#else
		cerr << "gzip not supported in your LCA build\n"; exit(50);
#endif
	} else { blast = new ifstream(inf.c_str(), ios::in); }
	if (!*blast) {
		cerr << "Blast input file " << inf << " could not be opened. exiting..\n";
		exit(23);
	}

	if (inFmt == "bl8") {
		inptFmt = 0;
	} else {
		cerr << "Unsupported input format '" << inFmt << "'. Only the custom bl8 format is supported.\n";
		exit(26);
	}
	lineBuffer.reserve(512);
	foundSbjs.reserve(256);
	completedQueries.reserve(1024);
	batchBuffer.reserve(256);

}

BlastReader::~BlastReader() {
	if (blast != NULL) {
		delete blast;
		blast = NULL;
	}
}

vector<BlastRes> BlastReader::getResBatch() {
	batchBuffer.clear();
	foundSbjs.clear();
	blastCnter++;
	auto readRecord = [&](BlastRes& result) -> bool {
		auto registerColumns = [&](int columns) {
			if (detectedColumns != 0 && detectedColumns != columns) {
				cerr << "Inconsistent m8 layout at line " << lineNumber << ": file started with "
					<< detectedColumns << " columns but this row has " << columns << ".\n";
				exit(25);
			}
			detectedColumns = columns;
			if (columns == 12 && !legacyNoticeShown) {
				cerr << "Notice: legacy 12-column BLAST input has no qlen; query-coverage filtering "
					<< "is skipped for these records.\n";
				legacyNoticeShown = true;
			}
		};
		while (getline(*blast, lineBuffer, '\n')) {
			lineNumber++;
			if (!lineBuffer.empty() && lineBuffer.back() == '\r') { lineBuffer.pop_back(); }
			const size_t contentStart = lineBuffer.find_first_not_of(" \t\r");
			if (contentStart == string::npos || lineBuffer[contentStart] == '#') { continue; }
			const int columns = BlastRes::supportedColumnCount(lineBuffer);
			if (!seenData && BlastRes::isColumnHeader(lineBuffer)) {
				registerColumns(columns);
				continue;
			}
			if (columns != 0) { registerColumns(columns); }
			if (!result.parseFromLine(lineBuffer, inptFmt)) {
				cerr << "Malformed m8 record at line " << lineNumber
					<< ". Expected either 11 columns ending in qlen, or the legacy 12-column "
					<< "BLAST layout ending in evalue and bitscore.\n";
				exit(25);
			}
			seenData = true;
			return true;
		}
		return false;
	};

	if (!processedBatch) {
		hasLastBlast = readRecord(lastBlast);
		processedBatch = true;
		if (!hasLastBlast) {
			allRead = true;
			return batchBuffer;
		}
	}

	if (!hasLastBlast || allRead || lastBlast.fail) {
		return batchBuffer;
	}

	batchBuffer.push_back(lastBlast);
	foundSbjs[lastBlast.Sbj] = 0;
	const string cmpQu = lastBlast.Query;

	BlastRes cur;
	while (readRecord(cur)) {
		if (cur.Query != cmpQu) {
			completedQueries.insert(cmpQu);
			if (completedQueries.find(cur.Query) != completedQueries.end()) {
				cerr << "Query '" << cur.Query << "' occurs in multiple non-contiguous blocks at line "
					<< lineNumber << ". Sort/group the custom m8 input by query.\n";
				exit(27);
			}
			lastBlast = cur;
			hasLastBlast = true;
			return batchBuffer;
		}

		auto existing = foundSbjs.find(cur.Sbj);
		if (existing == foundSbjs.end()) {
			foundSbjs[cur.Sbj] = batchBuffer.size();
			batchBuffer.push_back(cur);
		} else if (better_blast_hit(cur, batchBuffer[existing->second])) {
			batchBuffer[existing->second] = cur;
		}
	}

	completedQueries.insert(cmpQu);
	hasLastBlast = false;
	allRead = true;
	return batchBuffer;
}
