// LCA.cpp

#include "libload.h"
#include "RefTax.h"
#include "LCAimpl.h"
#include "Matrix.h"
#include "Pipeline.h"
#include <chrono>
#include <filesystem>
//0.24: fixed bug of not reading "k__?; p__?; c__?; .." strings
//0.26: 28.3.26: fixed various smallish bugs, including wrongly reported %id in some cases, and some parallelization issues.
//0.27: 4.9.26: performance improvements
//0.28: 21.7.26: deterministic hit filtering; corrected LCA/depth, hit-read output,
//taxonomy/rank and CRLF parsing, query-span coverage, legacy 12-column BLAST input,
//configurable tdep, strict input,
//single-thread execution, multi-database reporting, and output error handling.
//Identity cutoffs became inclusive (a hit at exactly 97% now reaches species; before, identity had
//to exceed the cutoff), and -id accepts decimals (before, values were truncated to integers).
//0.29: 8.9.26: fixed whitespace parsing bug
//0.30: 11.9.26: preserve spaces in tab-separated IDs; fix hit selection/filtering,
//identity/depth consistency, and input/output error handling.
//0.31: 8.10.26: highest identity anchors the primary filter, with the relative filter on query
//coverage instead of alignment length; matching names below a missing rank are used by the LCA;
//tied -reportBestHit hits resolve by consensus; malformed mapping/taxonomy lines are skipped with
//warnings; taxonomy columns after the lineage are ignored; more placeholder species names are
//recognized; query-block tracking is limited to 100,000 queries.
//0.32: 9.10.26: faster I/O and parsing, compact taxonomy storage (about 4x less memory), and -t N
//threads for taxonomy loading and assignment; output is identical to 0.31 for any thread count.
//0.33: 23.10.26: LCAlib.h, the assignment of one query as a header-only library (lca::Assign);
//the sources moved to src/;
//output is identical to 0.32.
const char* LCA_ver = "0.33";

namespace {
bool samePath(const string& first, const string& second) {
	std::error_code error;
	if (std::filesystem::equivalent(first, second, error)) { return true; }
	const auto canonicalFirst = std::filesystem::weakly_canonical(first, error);
	if (error) { return first == second; }
	const auto canonicalSecond = std::filesystem::weakly_canonical(second, error);
	return error ? first == second : canonicalFirst == canonicalSecond;
}

bool validateOutputPaths(const options& opt, const Matrix& matrix) {
	vector<string> inputs = opt.blFiles;
	inputs.insert(inputs.end(), opt.refDBs.begin(), opt.refDBs.end());
	vector<string> outputs = {opt.outF};
	if (!opt.repHitPattern.empty()) { outputs.push_back(opt.repHitPattern); }
	if (opt.calcHighMats) {
		const auto paths = matrix.outputPaths(opt.outF);
		outputs.insert(outputs.end(), paths.begin(), paths.end());
	}
	for (size_t i = 0; i < outputs.size(); ++i) {
		for (const auto& input : inputs) {
			if (samePath(outputs[i], input)) {
				cerr << "Output path " << outputs[i] << " overlaps input file " << input << endl;
				return false;
			}
		}
		for (size_t j = 0; j < i; ++j) {
			if (samePath(outputs[i], outputs[j])) {
				cerr << "Output paths overlap: " << outputs[i] << " and " << outputs[j] << endl;
				return false;
			}
		}
	}
	return true;
}

// Queries handed to a worker at once, and the text and assignments computed for them. The query
// vectors travel back with the result and are reused, so they keep their capacity.
struct QueryChunk {
	vector<vector<BlastRes> > queries;
	size_t count = 0;
};
struct ChunkResult {
	string text;
	string hitPattern;
	vector<TaxObj*> results;
	size_t written = 0;
	vector<vector<BlastRes> > queries;
};
// Threads need large chunks to amortize synchronization; one thread is faster with small ones,
// which keep freshly parsed hits in cache.
size_t chunkMaxQueries(int threads) { return threads > 1 ? 2048 : 16; }
size_t chunkMaxHits(int threads) { return threads > 1 ? 16384 : 256; }

int runLCA(options& OPT) {
	const auto start = std::chrono::steady_clock::now();
	const size_t refDbCount = OPT.refDBs.size();
	const bool highLvl = OPT.calcHighMats;
	Matrix mat(OPT.taxDepth, OPT.Taxlvls, OPT.hitRD);
	if (!validateOutputPaths(OPT, mat)) { return 29; }
	unordered_map<string, TaxObj*> assign;
	size_t queriesRead = 0;

	const string IDname = OPT.isReads ? "Reads" : "OTU";
	vector<char> outputBuffer(__inputBufferSize);
	ofstream O;
	O.rdbuf()->pubsetbuf(outputBuffer.data(), outputBuffer.size());
	O.open(OPT.outF.c_str());
	if (!O) {
		cerr << "Could not create output file " << OPT.outF << endl;
		return 30;
	}
	O << IDname << "\t" << OPT.TaxLvl2string();
	if (OPT.hitRD) { O << "\tHit2DB"; }
	if (OPT.reportID) { O << "\t%ID"; }
	O << '\n';

	vector<char> hitPatternBuffer;
	ofstream HITPAT;
	const bool checkHitPat = !OPT.repHitPattern.empty();
	if (checkHitPat) {
		hitPatternBuffer.resize(__inputBufferSize);
		HITPAT.rdbuf()->pubsetbuf(hitPatternBuffer.data(), hitPatternBuffer.size());
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
		RefTax RT(OPT.refDBs[xi], OPT.taxDepth, OPT.nativeSlVdb, OPT.checkTaxoUnkw, OPT.numThr);
		RT.setTaxLvls(OPT.Taxlvls);
		BlastReader BR(OPT.blFiles[xi], OPT.input_format);

		// The reader groups queries on this thread; workers assign whole chunks, and results are
		// written in input order, so the output does not depend on the number of threads.
		vector<vector<vector<BlastRes> > > recycled; // used on this thread only
		const size_t maxQueries = chunkMaxQueries(OPT.numThr);
		const size_t maxHits = chunkMaxHits(OPT.numThr);
		auto produce = [&BR, &recycled, maxQueries, maxHits](QueryChunk& chunk) {
			if (!recycled.empty()) {
				chunk.queries = std::move(recycled.back());
				recycled.pop_back();
			}
			size_t hits = 0;
			while (chunk.count < maxQueries && hits < maxHits) {
				if (chunk.count == chunk.queries.size()) { chunk.queries.emplace_back(); }
				if (!BR.nextBatch(chunk.queries[chunk.count])) { break; }
				hits += chunk.queries[chunk.count].size();
				chunk.count++;
			}
			return chunk.count > 0;
		};
		auto work = [&](QueryChunk& chunk, ChunkResult& out) {
			for (size_t i = 0; i < chunk.count; ++i) {
				TaxObj* result = LCA(chunk.queries[i], &RT, &OPT);
				if (result == NULL) { continue; }
				if (!multiDBuse) {
					out.text += result->Subj;
					out.text += '\t';
					out.text += result->getWriteString(OPT.idThr);
					out.text += '\n';
					if (checkHitPat) {
						out.hitPattern += result->Subj + '\t' + to_string(result->perID) + '\t' +
							to_string(result->depth) + '\n';
					}
					out.written++;
					if (!highLvl) {
						delete result;
						continue;
					}
				}
				out.results.push_back(result);
			}
			out.queries = std::move(chunk.queries);
		};
		auto consume = [&](ChunkResult& out) {
			O << out.text;
			if (checkHitPat) { HITPAT << out.hitPattern; }
			taxWritten += out.written;
			for (TaxObj* result : out.results) {
				if (!multiDBuse) {
					mat.add(result);
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
			out.results.clear();
			recycled.push_back(std::move(out.queries));
		};
		runOrdered<QueryChunk, ChunkResult>(OPT.numThr, produce, work, consume);
		queriesRead += BR.queriesRead();
		cout << "Done Blast File reading (" << BR.queriesRead() << " queries)\n";
	}

	vector<string> orderedQueries;
	orderedQueries.reserve(assign.size());
	for (const auto& entry : assign) { orderedQueries.push_back(entry.first); }
	std::sort(orderedQueries.begin(), orderedQueries.end());
	for (const auto& query : orderedQueries) {
		TaxObj* result = assign[query];
		O << result->Subj << '\t' << result->getWriteString(OPT.idThr) << '\n';
		if (checkHitPat) {
			HITPAT << result->Subj << '\t' << to_string(result->perID) << '\t' << result->depth << '\n';
		}
		taxWritten++;
		if (highLvl) { mat.add(result); }
		delete result;
		assign[query] = NULL;
	}

	O.close();
	if (!O) {
		cerr << "Failed while writing output file " << OPT.outF << endl;
		return 32;
	}
	if (checkHitPat) {
		HITPAT.close();
		if (!HITPAT) {
			cerr << "Failed while writing hit-pattern file " << OPT.repHitPattern << endl;
			return 33;
		}
	}
	if (multiDBuse) {
		cout << "Wrote " << taxWritten << " LCA tax assignments from " << queriesRead << " query blocks in "
			<< refDbCount << " reference databases\n";
	} else {
		cout << "Wrote " << taxWritten << "/" << queriesRead << " LCA tax assignments\n";
	}

	if (highLvl && !mat.writeAllLevels(OPT.outF)) {
		return 34;
	}

	if (duplicateQueries > 0) {
		cout << "Found " << duplicateQueries << " duplicate query assignments across " << refDbCount
			<< " reference databases; reassigned " << reassignedQueries << " of these." << endl;
	}
	printf("LCA finished. Time taken: %.2fs\n",
		std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
	return 0;
}
}

void welcomeMsg() {
	cout << "Least common ancestor (LCA) assignments ver " << LCA_ver << endl;
}
int main(int argc, char* argv[])
{
	options OPT(argc, argv, __default_depth);
	if (OPT.version) {
		cout << LCA_ver << endl;
		return 0;
	}
	welcomeMsg();
	try {
		return runLCA(OPT);
	} catch (const FatalError& error) {
		cerr << error.message << std::flush;
		return error.code;
	}
}
