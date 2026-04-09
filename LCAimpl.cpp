#include "LCAimpl.h"
#include <algorithm>
#include <cmath>

TaxObj* LCA(vector<BlastRes>& BR, RefTax* RT, options* opt) {
	//hard coded -> replace later for dynamic allocation

	if (BR.size() == 0) { return NULL; } // empty, just return empty string
	
	
	//primary blast filter
	double bestPerID(0.f);
   filterBlastPrimary(BR, opt, bestPerID);
	if (BR.size() == 0) { return NULL; } // empty, just return empty string

	float avgPerID(0.f);
   vector<TaxObj*> allTax = BlastToTax(BR, RT, opt, avgPerID);
	
	//routine that performs actual LCA matching etc
	//string ret((*BR.begin())->Query + "\t");
	
#ifdef LCAdebg
   cout << BR.front().Query << endl;
#endif
	bool singularHit(false);
	
	TaxObj* ret = LCAcore(allTax,singularHit,opt->LCAfract,opt->taxDepth) ;
	ret->setRepID(opt->reportID);
   ret->Subj = BR.front().Query;
	ret->perID = avgPerID;

	if (opt->hitRD ){
        auto bestIt = std::max_element(BR.begin(), BR.end(), [](const BlastRes& a, const BlastRes& b) {
			if (a.perID != b.perID) { return a.perID < b.perID; }
			return a.alLen < b.alLen;
		});
		if (bestIt != BR.end()) {
			ret->addHitDB(bestIt->Sbj);
			ret->perID = (float)bestIt->perID;
		} else {
			ret->addHitDB(__unkwnTax);
			ret->perID = (float)bestPerID;
		}
		/*if (opt->reportBestHit && (*BR.begin())->perID >= opt->idThr.back()) {
			ret->addHitDB((*BR.begin())->Sbj);
		} else if (singularHit) {
			ret->addHitDB((*BR.begin())->Sbj);
		} else {
			ret->addHitDB(__unkwnTax);
		}*/
	}
	
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
	
	int FullSetCnt((int)TO.size());
	int MaxSetCnt(FullSetCnt);
	if ( FullSetCnt == 1 ){//early return
       ret->copy_vals(TO.front());
		hitRd = true;
		return ret;
	}

   vector<char> DO(MaxSetCnt, 1);

	for (int DL(0); DL < tdepth; DL++) {
		unordered_map<string, int> cntOptions;
        cntOptions.reserve(MaxSetCnt + 1);
		cntOptions[__unkwnTax] = 0;
		bool goDeeper(false);
		//create set of counts at current level
        for (int pos = 0; pos < MaxSetCnt; ++pos) {
			if (!DO[pos]) { continue; }
			const string& curTax = TO[pos]->get(DL);
			auto fnd = cntOptions.find(curTax);
			if (fnd == cntOptions.end()) {
				cntOptions[curTax] = 1;
			} else {
				fnd->second++;
			}
#ifdef LCAdebg
			cout << curTax << " ";
#endif // LCAdebg
		}
		//some static variables
		int unknwCnt(cntOptions[__unkwnTax]);
		double MaxPossCons(double (FullSetCnt - unknwCnt));
        const string* consens = &__unkwnTax;

		//check if all unkown
		if (unknwCnt >= MaxSetCnt) {
			break;
		}
		//eval if consensus count was present

      int bestCnt = -1;
		for (auto ic = cntOptions.begin(); ic != cntOptions.end(); ++ic) {
			if (ic->first == __unkwnTax) { continue; }
			if (ic->second > bestCnt || (ic->second == bestCnt && ic->first < *consens)) {
				bestCnt = ic->second;
				consens = &ic->first;
			}
		}

		//best count only needs to be classified tax and reaching the e.g. 90% of max possible count
		if (bestCnt >= (MaxPossCons * LCAfrac)) {
			FullSetCnt = bestCnt; // reset full set cnt to current consensus set size for deeper levels
			goDeeper = true;
			ret->set(DL, *consens);
#ifdef LCAdebg
			cout << " : " << *consens;
#endif
		}

#ifdef LCAdebg
		cout << endl;
#endif
		if (!goDeeper) { break; }

		if (FullSetCnt == (MaxSetCnt - unknwCnt) ) {
			continue;
		}

		//while loop again, to exclude non consens assignments
     for (int pos = 0; pos < MaxSetCnt; ++pos) {
			if (!DO[pos]) { continue; }
			if (TO[pos]->get(DL) == *consens) { continue; }
			DO[pos] = 0;
		}
	}
	return ret;
}

double filterBlastPrimary(vector<BlastRes>& BR, options* opt, double& bestID) {
	
	if (!opt->BLfilter) {
        bestID = BR.front().perID;
	 return BR.front().perID * (double)BR.front().alLen;
	}

	float minCov = opt->minCover;
	int maxL(0);
	
   for (auto it = BR.begin(); it != BR.end();it++) {
		if (it->perID > bestID && it->alLen >= maxL * 0.95) { bestID = it->perID; maxL = it->alLen; }
		else if (it->perID > bestID * 0.9 && it->alLen >= maxL * 1.2) { bestID = it->perID; maxL = it->alLen; }
	}
	
	//filter parameters
	double lengthToleranceF(0.85f);
	double tolerance(1.5);
	if (opt->reportBestHit) {
       tolerance = 0.f; lengthToleranceF= 0.95f;
	}
	else if (bestID >= 100) { tolerance = 0.05f; }
	else if (bestID >= 99.5) { tolerance = 0.15f; }
	else if (bestID >= 99) { tolerance = 0.25f; }
	else if (bestID >= 98) { tolerance = 0.75f; }
	else if (bestID >= 97) { tolerance = 1.0f; }
	
	//precalc maxL
   maxL = static_cast<int>(std::ceil(static_cast<double>(maxL) * lengthToleranceF));
	const int minAliLen = static_cast<int>(std::ceil(opt->minAliLen));
	if (maxL < minAliLen) {
		maxL = minAliLen;
	}

	size_t writePos = 0;
	for (size_t readPos = 0; readPos < BR.size(); ++readPos) {
		BlastRes& cur = BR[readPos];
		if ((cur.perID + tolerance) < bestID ||
			cur.alLen < maxL ||
			cur.Qcoverage < minCov) {
			continue;
		}
		if (writePos != readPos) {
			BR[writePos] = cur;
		}
		++writePos;
	}
	BR.resize(writePos);
	//cerr << "\n";

	//filter done
	return bestID*(double)maxL;
}


vector<TaxObj*> BlastToTax(const vector<BlastRes>& BR, RefTax* RT, options* opt, float& consPerID) {
	vector<TaxObj*> ret;
	ret.reserve(BR.size());
	vector<double> &thr = opt->idThr;


	int depth = RT->depth();
	bool anySpeciesCertain(false);
  for (auto it = BR.begin(); it != BR.end(); it++) {
		double curID(it->perID);
		consPerID += (float)curID;
        auto fnd = RT->find(it->Sbj);
		if (fnd != RT->end()) {
			TaxObj* F = new TaxObj(fnd->second);
			int maxD(0);
			//assign max depth based on % id to subject
			while (maxD < depth && thr[maxD] < curID ) { maxD++; }
			F->depth = maxD;
           ret.emplace_back(F);
			if (!F->speciesUncertain) { anySpeciesCertain = true; }
		}	else {
            cerr << "Could not find tax for Subject " << it->Sbj << endl;
			exit(74);
		}
	}
	consPerID /= BR.size();

	//remove uncertain species, in case of enough good hits
  if (anySpeciesCertain && ret.size() > 1) {
		for (auto* t : ret) {
            if (t->speciesUncertain) {
				t->makeSpeciesUnknown();
			}
		}
	}


	return ret;
}