#pragma once
#include "libload.h"
#include "RefTax.h"
inline const string __taxSepMat = ";";
inline const string __MatSep = "\t";

class Matrix
{
public:
	Matrix(int depth, vector<string> taxs, bool reportRead);
	~Matrix();
	void add(TaxObj*);
	vector<string> outputPaths(const string&) const;
	bool writeAllLevels(const string&);
//vars
	vector< string > colIDs; //rows = features, cols = tax levels
	vector<unordered_map<string, size_t> > rowIDs;
	bool readReport;
};

