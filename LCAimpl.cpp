#include "LCAimpl.h"
#include <algorithm>
#include <cmath>

namespace {
void deduplicateBlastHits(vector<BlastRes>& BR) {
	unordered_map<string, size_t> subjects;
	size_t retained = 0;
	for (size_t pos = 0; pos < BR.size(); ++pos) {
		const BlastRes& hit = BR[pos];
		auto found = subjects.find(hit.Sbj);
		if (found == subjects.end()) {
			subjects[hit.Sbj] = retained;
			if (retained != pos) { BR[retained] = hit; }
			++retained;
		} else if (hit.betterThan(BR[found->second])) {
			BR[found->second] = hit;
		}
	}
	BR.resize(retained);
}
}

TaxObj* LCA(vector<BlastRes>& BR, RefTax* RT, options* opt) {
	//hard coded -> replace later for dynamic allocation

	if (BR.size() == 0) { return NULL; } // empty, just return empty string
	
	
	//primary blast filter
	double bestPerID(0.f);
   filterBlastPrimary(BR, opt, bestPerID);
	if (BR.size() == 0) { return NULL; } // empty, just return empty string

	double avgPerID(0.0);
   vector<TaxObj*> allTax = BlastToTax(BR, RT, opt, avgPerID);
	
	//routine that performs actual LCA matching etc
	//string ret((*BR.begin())->Query + "\t");
	
#ifdef LCAdebg
   cout << BR.front().Query << endl;
#endif
	bool singularHit(false);
	
	TaxObj* ret = LCAcore(allTax,singularHit,opt->LCAfract,opt->taxDepth);
	ret->Subj = BR.front().Query;
	ret->perID = avgPerID;
	ret->limitByIdentity(avgPerID, opt->idThr);

	if (opt->hitRD) {
		auto bestIt = std::min_element(BR.begin(), BR.end(), [](const BlastRes& a, const BlastRes& b) {
			return a.betterThan(b);
		});
		const bool reportSubject = bestIt != BR.end() && (singularHit || opt->reportBestHit);
		ret->setHitDB(reportSubject ? bestIt->Sbj : __unkwnTax);
	}
	ret->setRepID(opt->reportID);
	
/*	if (allTax.size() > 1) {
		int ii = 0;
	}*/
	
	//finalize 
	//ret += "\n";

	//cleanup
  for (auto* x : allTax) {
		delete x;
	}
	return ret;
}

TaxObj* LCAcore(const vector<TaxObj*>& TO, bool &hitRd, double LCAfrac, int tdepth) {
	TaxObj* ret = new TaxObj(tdepth);
	hitRd = false;
	const int maxSetCount = static_cast<int>(TO.size());
	if (maxSetCount == 1) {
		ret->copy_vals(TO.front());
		hitRd = true;
		return ret;
	}
	if (maxSetCount == 0) { return ret; }

	vector<char> active(maxSetCount, 1);

	for (int DL = 0; DL < tdepth; DL++) {
		unordered_map<string, int> cntOptions;
		cntOptions.reserve(maxSetCount + 1);
		int activeCount = 0;
		int unknownCount = 0;
		for (int pos = 0; pos < maxSetCount; ++pos) {
			if (!active[pos]) { continue; }
			activeCount++;
			const string& curTax = TO[pos]->get(DL);
			if (curTax == __unkwnTax) { unknownCount++; }
			else { cntOptions[curTax]++; }
#ifdef LCAdebg
			cout << curTax << " ";
#endif // LCAdebg
		}
		const int knownCount = activeCount - unknownCount;
		if (knownCount <= 0) { break; }

		int bestCount = 0;
		string consensus;
		for (const auto& option : cntOptions) {
			if (option.second > bestCount ||
				(option.second == bestCount && (consensus.empty() || option.first < consensus))) {
				bestCount = option.second;
				consensus = option.first;
			}
		}

		if (bestCount < static_cast<double>(knownCount) * LCAfrac) { break; }
		ret->set(DL, consensus);
#ifdef LCAdebg
		cout << " : " << consensus;
#endif

#ifdef LCAdebg
		cout << endl;
#endif
		// Only descendants of the accepted parent can vote at the next level.
		for (int pos = 0; pos < maxSetCount; ++pos) {
			if (active[pos] && TO[pos]->get(DL) != consensus) { active[pos] = 0; }
		}
	}
	return ret;
}

double filterBlastPrimary(vector<BlastRes>& BR, options* opt, double& bestID) {
	bestID = 0.0;
	if (BR.empty()) { return 0.0; }
	const float minCov = opt->minCover;
	auto passesCoverage = [minCov](const BlastRes& hit) {
		return !hit.queryCoverageKnown || hit.Qcoverage >= minCov;
	};
	const int minAliLen = static_cast<int>(std::ceil(opt->minAliLen));

	if (opt->BLfilter) {
		BR.erase(std::remove_if(BR.begin(), BR.end(), [&](const BlastRes& hit) {
			return hit.alLen < minAliLen || !passesCoverage(hit);
		}), BR.end());
	}
	if (BR.empty()) { return 0.0; }

	if (opt->reportBestHit) {
		const auto best = std::min_element(BR.begin(), BR.end(), [](const BlastRes& a, const BlastRes& b) {
			return a.betterThan(b);
		});
		if (best != BR.begin()) { BR.front() = *best; }
		BR.resize(1);
	}
	if (!opt->BLfilter || opt->reportBestHit) {
		const auto best = std::min_element(BR.begin(), BR.end(), [](const BlastRes& a, const BlastRes& b) {
			return a.betterThan(b);
		});
		bestID = best->perID;
		const double score = bestID * static_cast<double>(best->alLen);
		deduplicateBlastHits(BR);
		return score;
	}

	const BlastRes* highestIdentity = &*std::min_element(BR.begin(), BR.end(), [](const BlastRes& a, const BlastRes& b) {
		return a.betterThan(b);
	});

	// Order-independent version of the former running compromise: prefer the
	// highest-identity hit unless a hit within 10% identity is at least 20% longer.
	const BlastRes* anchor = highestIdentity;
	const double substantialLength = static_cast<double>(highestIdentity->alLen) * 1.2;
	for (const auto& hit : BR) {
		if (static_cast<double>(hit.alLen) < substantialLength || hit.perID < highestIdentity->perID * 0.9) { continue; }
		if (anchor == highestIdentity || hit.betterThan(*anchor)) { anchor = &hit; }
	}
	bestID = anchor->perID;
	int maxL = anchor->alLen;

	//filter parameters
	double lengthToleranceF(0.85);
	double tolerance(1.5);
	if (bestID >= 100) { tolerance = 0.05f; }
	else if (bestID >= 99.5) { tolerance = 0.15f; }
	else if (bestID >= 99) { tolerance = 0.25f; }
	else if (bestID >= 98) { tolerance = 0.75f; }
	else if (bestID >= 97) { tolerance = 1.0f; }
	
	maxL = static_cast<int>(std::ceil(static_cast<double>(maxL) * lengthToleranceF));
	if (maxL < minAliLen) {
		maxL = minAliLen;
	}

	size_t writePos = 0;
	for (size_t readPos = 0; readPos < BR.size(); ++readPos) {
		BlastRes& cur = BR[readPos];
		if ((cur.perID + tolerance) < bestID ||
			cur.alLen < maxL ||
			!passesCoverage(cur)) {
			continue;
		}
		if (writePos != readPos) {
			BR[writePos] = cur;
		}
		++writePos;
	}
	BR.resize(writePos);
	// Only surviving alignments may represent a subject in the LCA vote.
	deduplicateBlastHits(BR);
	//cerr << "\n";

	//filter done
	return bestID*(double)maxL;
}


vector<TaxObj*> BlastToTax(const vector<BlastRes>& BR, RefTax* RT, options* opt, double& consPerID) {
	vector<TaxObj*> ret;
	ret.reserve(BR.size());
	vector<double> &thr = opt->idThr;


	consPerID = 0.0f;
	if (BR.empty()) { return ret; }
	long double identitySum = 0.0L;
	bool anySpeciesCertain(false);
	for (auto it = BR.begin(); it != BR.end(); it++) {
		double curID(it->perID);
		identitySum += curID;
		auto fnd = RT->find(it->Sbj);
		if (fnd != RT->end()) {
			TaxObj* F = new TaxObj(fnd->second);
			F->limitByIdentity(curID, thr);
			ret.emplace_back(F);
			if (F->get(6) != __unkwnTax && !F->speciesUncertain) { anySpeciesCertain = true; }
		}	else {
            cerr << "Could not find tax for Subject " << it->Sbj << endl;
			exit(74);
		}
	}
	consPerID = static_cast<double>(identitySum / BR.size());

	//remove uncertain species, in case of enough good hits
	if (anySpeciesCertain && ret.size() > 1) {
		for (auto* t : ret) {
			if (t->depth > 6 && t->speciesUncertain) {
				t->makeSpeciesUnknown();
			}
		}
	}


	return ret;
}
