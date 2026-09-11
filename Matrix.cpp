#include "Matrix.h"



Matrix::Matrix(int depth, vector<string> taxs, bool reportRead):mat(depth), 
		colIDs(taxs), rowIDs(depth), readReport(reportRead)
{
	if (readReport){//means that not only species level is reported, but also singular hit (shini table)
		mat.resize(depth + 1);
		rowIDs.resize(depth + 1); colIDs.resize(depth + 1, "Hit2DB");
	}
}
Matrix::~Matrix()
{
}
void Matrix::add(TaxObj* t) {
	string addStr("");// t->get(0));
	for (size_t DL = 0; DL < rowIDs.size(); DL++) {
		string curT = readReport && DL + 1 == rowIDs.size()
			? t->getHitDB() : t->get(static_cast<int>(DL));
		if (curT == __unkwnTax) {curT = __unkwnTaxWR;}
		if (DL == 0) {
			addStr += curT;
		} else {
			addStr += __taxSepMat + curT;
		}
		auto fnd = rowIDs[DL].find(addStr);
		if (fnd == rowIDs[DL].end()) {
			rowIDs[DL][addStr] = 1;
		} else {
			rowIDs[DL][addStr]++;
		}
	}
}
vector<string> Matrix::outputPaths(const string& outF) const {
	vector<string> paths;
	for (const auto& level : colIDs) { paths.push_back(outF + "_" + level); }
	return paths;
}

bool Matrix::writeAllLevels(const string& outF) {
	const auto paths = outputPaths(outF);
	for (size_t DL = 0; DL < rowIDs.size(); DL++) {
		const string& outF1 = paths[DL];
		ofstream of(outF1.c_str());
		if (!of) {
			cerr << "Could not create matrix output file " << outF1 << endl;
			return false;
		}
		vector<pair<string, size_t> > rows(rowIDs[DL].begin(), rowIDs[DL].end());
		std::sort(rows.begin(), rows.end(), [](const pair<string, size_t>& a, const pair<string, size_t>& b) {
			return a.first < b.first;
		});
		for (const auto& row : rows) {
			of << row.first << __MatSep << row.second << '\n';
		}
		of.close();
		if (!of) {
			cerr << "Failed while writing matrix output file " << outF1 << endl;
			return false;
		}
	}
	return true;
}
