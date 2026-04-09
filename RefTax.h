#pragma once
#include "libload.h"
#include "options.h"
void trim(string& str,const std::string& whitespace = " \t");
bool isGZfile(const string fi);


struct TaxObj
{
	TaxObj(const string&, int, bool nativeSLV, bool doNotCheckTax);
	TaxObj(TaxObj*t);
  TaxObj(int d) :SavedTaxs(), Subj(""), perID(0.f), speciesUncertain(false),depth(d) {}
	//functions
	string getWriteString(const vector<double>&);
	void copy_vals(TaxObj*t) { SavedTaxs = t->SavedTaxs; depth = t->depth; }
	void setRepID(bool x) { repID = x; }
   void makeSpeciesUnknown() {
		if (speciesUncertain && depth > 0) {
			if ((int)SavedTaxs.size() < depth) { SavedTaxs.resize(depth, __unkwnTax); }
			SavedTaxs[depth - 1] = __unkwnTax;
		}
	}
	//get tax at depth x
  string& get(int x) {
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
	bool evalAcpyTax(TaxObj* oth);
	inline void copyOver(TaxObj* oth);
   void addHitDB(string x) { SavedTaxs.push_back(x); depth = static_cast<int>(SavedTaxs.size()); }

	//int dept() { return depth; }
	//variables
	vector<string> SavedTaxs;
	string Subj;
	float perID;
	bool repID;
	bool speciesUncertain;
	int depth;//saves explicitly the depth, taking ? etc into account

};


class RefTax
{
public:
	RefTax(const string&,int tdep,bool,bool);
	~RefTax();
	void stats();
	int depth() { return (int) tlevels.size(); }
	unordered_map <string, TaxObj*>::const_iterator find(const string& s) {
		return Tlink.find(s);
	}
	unordered_map <string, TaxObj*>::const_iterator end() {
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
	static bool extractQueryToken(const string&, string&);
	bool isSameQuery(const string &q) const {	if (q == Query) { return true; } return false;	}
	
	string Query, Sbj;
	int alLen;
	double perID, eval,score;
	float Qcoverage;
	bool fail;
};

class BlastReader
{
public:
	BlastReader(const string&, const string&);
 ~BlastReader();
  vector<BlastRes> getResBatch();
private:

	bool openedGZ,processedBatch;
    bool hasLastBlast;
	BlastRes lastBlast;
	istream* blast;
	bool allRead;
	int inptFmt;
	int blastCnter;
  string lineBuffer;
	unordered_set<string> foundSbjs;
	vector<BlastRes> batchBuffer;
};
