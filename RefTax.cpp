#include "RefTax.h"
#include <cctype>

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
	const char* first = src.data() + begin;
	const char* last = src.data() + end;
   char* parseEnd = nullptr;
	out = std::strtod(first, &parseEnd);
	return parseEnd == last;
}

inline bool parse_int_range(const string& src, size_t begin, size_t end, int& out) {
	const char* first = src.data() + begin;
	const char* last = src.data() + end;
   char* parseEnd = nullptr;
	long parsed = std::strtol(first, &parseEnd, 10);
	if (parseEnd != last) { return false; }
	out = static_cast<int>(parsed);
	return true;
}

inline bool query_matches_line(const string& line, const string& query) {
	const size_t sep = line.find('\t');
	if (sep == string::npos || sep != query.size()) { return false; }
	return line.compare(0, sep, query) == 0;
}
}

//generic functions
void trim(string& str,
	const std::string& whitespace)
{
	auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		strBegin = 0;

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	str = str.substr(strBegin, strRange);
}

bool isGZfile(const string fi) {
	string subst = fi.substr(fi.length() - 3);
	if (subst == ".gz") {
		return true;
	}
	return false;
}



TaxObj::TaxObj(TaxObj* t):SavedTaxs(t->SavedTaxs), Subj(t->Subj), perID(t->perID),
speciesUncertain(t->speciesUncertain),depth(t->depth) {
	;
}

string TaxObj::getWriteString(const vector<double>& ids) {
   int outDepth = (int)ids.size();
	if (outDepth <= 0) {
		outDepth = depth;
	}
	if (outDepth <= 0) {
		return "";
	}
	string ret(get(0));
	for (int i = 1; i < outDepth; i++) {
		string& cur = get(i);
		if (i >= depth || cur == __unkwnTax || (i < 7 && i < (int)ids.size() && perID < ids[i])) {
			ret += __defaultTaxSep + __unkwnTaxWR;
		} else {
          ret += __defaultTaxSep + cur;
		}
	}
	if (repID) {
		ret += __defaultTaxSep + to_string(perID);
	}
	return ret;
}



TaxObj::TaxObj(const string& X,int d, bool nativeSLV, bool doNotCheckTax):
   SavedTaxs(), Subj(""), perID(0.f), speciesUncertain(false),depth(0){
	SavedTaxs.reserve(d);
 size_t fnd = nativeSLV ? 0 : X.find("__", 0);
	if (!nativeSLV) {
		fnd = (fnd == string::npos) ? string::npos : fnd + 2;
	}

	const size_t strL = X.length();
	int cnt(0);
	while (fnd != string::npos && cnt < d) {
		size_t f2 = X.find(";", fnd);
		size_t end = (f2 == string::npos) ? strL : f2;

		size_t begin = fnd;
		while (begin < end && (X[begin] == ' ' || X[begin] == '\t')) { ++begin; }
		while (end > begin && (X[end - 1] == ' ' || X[end - 1] == '\t')) { --end; }
        const size_t tokenLen = end - begin;

        bool taxKnown = tokenLen > 0 && !iequals_ascii_range(X, begin, end, "unclassified") && !iequals_ascii_range(X, begin, end, "uncultured bacterium")
			&& !iequals_ascii_range(X, begin, end, "uncultured") && !(tokenLen == 1 && X[begin] == '?');

		if (taxKnown && cnt == 6) {
           size_t pos = X.find(' ', begin);
			if (pos != string::npos && pos < end) {
				pos = X.find(' ', pos + 1);
				if (pos == string::npos || pos >= end) {
					pos = string::npos;
				}
			}
            size_t pos2 = X.find("sp.", begin);
			if (pos2 == string::npos || pos2 >= end) {
				pos2 = string::npos;
			}
			if (istarts_with_ascii_range(X, begin, end, "uncultured") || istarts_with_ascii_range(X, begin, end, "unclassified")) {
				speciesUncertain = true;
           } else if (pos != string::npos && pos2 == pos - 3) {
				const size_t cmpLen = pos2 - begin - 1;
              if (cnt > 0 && cnt - 1 < (int)SavedTaxs.size() && SavedTaxs[cnt - 1].size() == cmpLen && SavedTaxs[cnt - 1].compare(0, cmpLen, X, begin, cmpLen) == 0) {
				speciesUncertain = true;
              }
			} else if (pos != string::npos) {
				end = pos;
			}
		}

		if (doNotCheckTax || taxKnown) {
			depth = cnt + 1;
         if ((int)SavedTaxs.size() <= cnt) {
				SavedTaxs.resize(cnt + 1, __unkwnTax);
			}
          SavedTaxs[cnt].assign(X, begin, end - begin);
		}

		cnt++;
		if (!nativeSLV) {
			if (f2 == string::npos) { break; }
			size_t next = X.find("__", f2);
			fnd = (next == string::npos) ? string::npos : next + 2;
		}
		else {
			if (f2 == string::npos || f2 + 1 >= strL) { break; }
			fnd = f2 + 1;
		}
	}
   if ((int)SavedTaxs.size() > depth) {
		SavedTaxs.resize(depth);
	}
}

bool TaxObj::evalAcpyTax(TaxObj* oth) {
	if (perID != 0.f && oth->perID != 0.f) {
		if (perID > (oth->perID *0.99) ) {
			return false;
		}
		else if ( (oth->perID * 0.985) > (perID ) ) { // real advantage
			if (oth->depth >= depth) { // also make sure that this is not a hit to "?"
				copyOver(oth);
				return true;
			}
		}
	}
	if (oth->depth > depth) {//replace tax
		copyOver(oth);
		return true;
	}
	return false;
}

void TaxObj::copyOver(TaxObj* oth) {
	SavedTaxs = oth->SavedTaxs;
	depth = oth->depth;
	perID = oth->perID;
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
	int TaxDbl(0), TaxSingl(0);
	while (getline(in, line, '\n')) {
		size_t dlmt = line.find("\t");
		if (dlmt == std::string::npos) {
			cerr << "Line " << line << " does not contain \\t character.. ignoring\n";
			continue;
		}
		string ID = line.substr(0, dlmt);
		//string ttax = line.substr(dlmt+1);
		TaxObj* t = new TaxObj(line.substr(dlmt + 1),7, nativeSLV, !checktaxStr);
		auto fnd = Tlink.find(ID);
		if (fnd == Tlink.end()){//all good
			Tlink[ID] = t; TaxSingl++;
		} else {//not good: tax is double annotated
			//cerr << "Tax ID: " << ID << " is double used!\n";
			//exit(930);
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
	int cnt = 0; vector<int> hist(10, 0);
	int maxD = 0;
	for (auto it = Tlink.begin(); it != Tlink.end(); ++it) {
		int dep = it->second->depth ; if (maxD < dep) { maxD = dep; }
		if (dep > 10) { cerr << "Tax depth " << dep << " of object " << it->first << " is too great!\n"; }
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
	Qcoverage(0.f), fail(true) {
}

BlastRes::BlastRes(const string& line, int inptFmt):
		Query(""), Sbj(""), alLen(0), perID(0.f), eval(-1.f), score(0.f),
		Qcoverage(0.f), fail(true) {
	parseFromLine(line, inptFmt);
}

bool BlastRes::extractQueryToken(const string& line, string& query) {
	if (line.empty()) { return false; }
	const size_t sep = line.find('\t');
	if (sep == string::npos) { return false; }
	query.assign(line, 0, sep);
	return true;
}

bool BlastRes::parseFromLine(const string& line, int inptFmt) {
	(void)inptFmt;
	fail = true;
	Query.clear();
	Sbj.clear();
	alLen = 0;
	perID = 0.0;
	Qcoverage = 0.f;

 if (line.empty()) { return false; }

	size_t lineEnd = line.size();
	if (lineEnd > 0 && line[lineEnd - 1] == '\r') {
		lineEnd--;
	}
	if (lineEnd == 0) { return false; }

	size_t cursor = 0;
	size_t begin = 0;
	size_t end = 0;
	auto nextField = [&](size_t& outBegin, size_t& outEnd) -> bool {
     if (cursor > lineEnd) { return false; }
		outBegin = cursor;
		outEnd = line.find('\t', cursor);
       if (outEnd == string::npos || outEnd > lineEnd) {
			outEnd = lineEnd;
			cursor = lineEnd + 1;
		}
		else {
			cursor = outEnd + 1;
		}
		return outBegin <= outEnd;
	};

	if (!nextField(begin, end)) { return false; }
	Query.assign(line, begin, end - begin);
	if (!nextField(begin, end)) { return false; }
	Sbj.assign(line, begin, end - begin);
	if (!nextField(begin, end) || !parse_double_range(line, begin, end, perID)) { return false; }
	if (!nextField(begin, end) || !parse_int_range(line, begin, end, alLen)) { return false; }
	//mismatch
	if (!nextField(begin, end)) { return false; }
	//inserts
	if (!nextField(begin, end)) { return false; }
	//qstart
	if (!nextField(begin, end)) { return false; }
	//qstop
	if (!nextField(begin, end)) { return false; }
	//sstart
	if (!nextField(begin, end)) { return false; }
	//stop
	if (!nextField(begin, end)) { return false; }
	//qlen
	if (!nextField(begin, end)) { return false; }
	double ql = 0.0;
	if (!parse_double_range(line, begin, end, ql)) { return false; }

	Qcoverage = 1.f; //could be eval in old m8 format..
	if (ql > 1) {
		Qcoverage = static_cast<float>(static_cast<double>(alLen) / ql);
	}
	fail = false;
	return true;
}


//*******************************************************
//        BlastReader
//*******************************************************


BlastReader::BlastReader(const string& inf, const string& inFmt):openedGZ(false), processedBatch(false),
hasLastBlast(false), blast(NULL), allRead(false), inptFmt(-1), blastCnter(0), lineBuffer(), foundSbjs(), batchBuffer(){
#ifdef DEBUG
	cerr << "ini blast file\n";
#endif // DEBUG

	if (isGZfile(inf)) {
		openedGZ = true;
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
	} else if (inFmt == "uc") {
		inptFmt = 1;
	}
	lineBuffer.reserve(512);
	foundSbjs.reserve(256);
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

	if (!processedBatch) {
		while (getline(*blast, lineBuffer, '\n')) {
			if (lineBuffer.empty()) { continue; }
			if (lastBlast.parseFromLine(lineBuffer, inptFmt)) {
				hasLastBlast = true;
				break;
			}
		}
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
	foundSbjs.insert(lastBlast.Sbj);
	const string cmpQu = lastBlast.Query;

	while (getline(*blast, lineBuffer, '\n')) {
		if (lineBuffer.empty()) { continue; }
		string lineQuery;
		if (!BlastRes::extractQueryToken(lineBuffer, lineQuery)) { continue; }

     if (lineQuery != cmpQu) {
			hasLastBlast = lastBlast.parseFromLine(lineBuffer, inptFmt);
			if (!hasLastBlast) { continue; }
			return batchBuffer;
		}

      batchBuffer.emplace_back();
		BlastRes& cur = batchBuffer.back();
		if (!cur.parseFromLine(lineBuffer, inptFmt)) {
			batchBuffer.pop_back();
			continue;
		}
		if (!foundSbjs.insert(cur.Sbj).second) {
			batchBuffer.pop_back();
		}
	}

	hasLastBlast = false;
	allRead = true;
	return batchBuffer;
}
