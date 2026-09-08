// LCA.cpp 

#include "libload.h"
#include "RefTax.h"
#include "LCAimpl.h"
#include "Matrix.h"
//0.24: fixed bug of not reading "k__?; p__?; c__?; .." strings
//0.26: 28.3.26: fixed various smallish bugs, including wrongly reported %id in some cases, and some parallelization issues. 
//0.27: 4.9.26: performance improvements
//0.28: 21.7.26: deterministic hit filtering; corrected LCA/depth, hit-read output,
//taxonomy/rank and CRLF parsing, query-span coverage, legacy 12-column BLAST input,
//configurable tdep, strict input,
//single-thread execution, multi-database reporting, and output error handling.
//0.29: 8.9.26: fixed whitespace parsing bug
const char* LCA_ver = "0.29";

void helpMsg() {
	cout << "LCA requires at least 3 arguments (-i, -r, -o)\n For more help and options, use \"./LCA -h\"\n";
#ifdef _gzipread
	cout << "Compiled with gzip support\n";
#else
	cout << "No gzip support compiled\n";
#endif
}
void welcomeMsg() {
	cout << "Least common ancestor (LCA) assignments ver " << LCA_ver << endl;
}
int main(int argc, char* argv[])
{
	clock_t tStart = clock();
	options OPT(argc, argv, __default_depth);
	if (OPT.version) {
		cout << LCA_ver << endl;
		return 0;
	}
	welcomeMsg();

	const size_t refDbCount = OPT.refDBs.size();
	const bool highLvl = OPT.calcHighMats;
	Matrix mat(OPT.taxDepth, OPT.Taxlvls, OPT.hitRD);
	unordered_map<string, TaxObj*> assign;
	unordered_set<string> inputQueries;

	const string IDname = OPT.isReads ? "Reads" : "OTU";
	ofstream O(OPT.outF.c_str());
	if (!O) {
		cerr << "Could not create output file " << OPT.outF << endl;
		return 30;
	}
	O << IDname << "\t" << OPT.TaxLvl2string();
	if (OPT.hitRD) { O << "\tHit2DB"; }
	if (OPT.reportID) { O << "\t%ID"; }
	O << '\n';

	ofstream HITPAT;
	const bool checkHitPat = !OPT.repHitPattern.empty();
	if (checkHitPat) {
		HITPAT.open(OPT.repHitPattern.c_str());
		if (!HITPAT) {
			cerr << "Could not create hit-pattern file " << OPT.repHitPattern << endl;
			return 31;
		}
		HITPAT << "OTU/ASV\tID\tTaxDepth\n";
	}

	size_t duplicateQueries = 0;
	size_t reassignedQueries = 0;
	size_t taxWritten = 0;
	const bool multiDBuse = refDbCount > 1;
	for (size_t xi = 0; xi < refDbCount; xi++) {
		RefTax RT(OPT.refDBs[xi], OPT.taxDepth, OPT.nativeSlVdb, OPT.checkTaxoUnkw);
		RT.setTaxLvls(OPT.Taxlvls);
		BlastReader BR(OPT.blFiles[xi], OPT.input_format);

		while (true) {
			vector<BlastRes> hits = BR.getResBatch();
			if (hits.empty()) { break; }
			inputQueries.insert(hits.front().Query);
			TaxObj* result = LCA(hits, &RT, &OPT);
			if (result == NULL) { continue; }

			if (!multiDBuse) {
				O << result->Subj << '\t' << result->getWriteString(OPT.idThr) << '\n';
				if (checkHitPat) {
					HITPAT << result->Subj << '\t' << result->perID << '\t' << result->depth << '\n';
				}
				taxWritten++;
				if (highLvl) { mat.add(result); }
				delete result;
				continue;
			}

			auto previous = assign.find(result->Subj);
			if (previous == assign.end()) {
				assign[result->Subj] = result;
			} else {
				duplicateQueries++;
				if (previous->second->evalAcpyTax(result)) { reassignedQueries++; }
				delete result;
			}
		}
		cout << "Done Blast File reading\n";
	}

	vector<string> orderedQueries;
	orderedQueries.reserve(assign.size());
	for (const auto& entry : assign) { orderedQueries.push_back(entry.first); }
	std::sort(orderedQueries.begin(), orderedQueries.end());
	for (const auto& query : orderedQueries) {
		TaxObj* result = assign[query];
		O << result->Subj << '\t' << result->getWriteString(OPT.idThr) << '\n';
		if (checkHitPat) {
			HITPAT << result->Subj << '\t' << result->perID << '\t' << result->depth << '\n';
		}
		taxWritten++;
		if (highLvl) { mat.add(result); }
		delete result;
		assign[query] = NULL;
	}

	O.flush();
	if (!O) {
		cerr << "Failed while writing output file " << OPT.outF << endl;
		return 32;
	}
	if (checkHitPat) {
		HITPAT.flush();
		if (!HITPAT) {
			cerr << "Failed while writing hit-pattern file " << OPT.repHitPattern << endl;
			return 33;
		}
	}
	cout << "Wrote " << taxWritten << "/" << inputQueries.size() << " LCA tax assignments\n";

	if (highLvl && !mat.writeAllLevels(OPT.outF)) {
		return 34;
	}

	if (duplicateQueries > 0) {
		cout << "Found " << duplicateQueries << " duplicate query assignments across " << refDbCount
			<< " reference databases; reassigned " << reassignedQueries << " of these." << endl;
	}
	printf("LCA finished. Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	return 0;
}
