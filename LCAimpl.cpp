#include "LCAimpl.h"

TaxObj* LCA(list<BlastRes*> BR, RefTax* RT, options* opt) {
	//hard coded -> replace later for dynamic allocation

	if (BR.size() == 0) { return NULL; } // empty, just return empty string
	
	
	//primary blast filter
	double bestPerID(0.f);
	double bestScore = filterBlastPrimary(BR,opt, bestPerID);
	if (BR.size() == 0) { return NULL; } // empty, just return empty string

	float avgPerID(0.f);
	list<TaxObj*> allTax = BlastToTax(BR, RT, opt, avgPerID);
	
	//routine that performs actual LCA matching etc
	//string ret((*BR.begin())->Query + "\t");
	
#ifdef LCAdebg
	cout << (*BR.begin())->Query << endl;
#endif
	bool singularHit(false);
	
	TaxObj* ret = LCAcore(allTax,singularHit,opt->LCAfract,opt->taxDepth) ;
	ret->setRepID(opt->reportID);
	ret->Subj = (*BR.begin())->Query;
	ret->perID = avgPerID;

	if (opt->hitRD ){
		ret->addHitDB((*BR.begin())->Sbj);
		ret->perID = (float) bestPerID;
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
	//for (auto&& child : BR) {
	for (list<BlastRes*>::iterator x = BR.begin(); x != BR.end();x++){
		delete *x;
	}
	BR.clear();
	for (list<TaxObj*>::iterator x = allTax.begin(); x != allTax.end(); x++) {
	//for (auto&& child : allTax) {
		delete *x;
	}
	allTax.clear();
	return ret;
}

TaxObj* LCAcore(list<TaxObj*> TO, bool &hitRd, double LCAfrac, int tdepth) {
	TaxObj* ret = new TaxObj(tdepth);
	
	int FullSetCnt((int)TO.size());
	int MaxSetCnt(FullSetCnt);
	if ( FullSetCnt == 1 ){//early return
		ret->copy_vals( *TO.begin());
		hitRd = true;
		return ret;
	}

	vector<bool> DO(FullSetCnt, true);

	for (int DL(0); DL < tdepth; DL++) {
		unordered_map<string, int> cntOptions;
		cntOptions[__unkwnTax] = 0;
		bool goDeeper(false);
		//create set of counts at current level
		int pos(-1);
		for (auto it = TO.begin(); it != TO.end(); it++) {
			pos++;  if (!DO[pos]) { continue; }
			string curTax = (*it)->get(DL);
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
		string consens (__unkwnTax);

		//check if all unkown
		if (unknwCnt >= MaxSetCnt) {
			break;
		}
		//eval if consensus count was present

		for (auto ic = cntOptions.begin(); ic != cntOptions.end(); ic++) {
			if (ic->first == __unkwnTax) { continue; }
			//best count only needs to be classified tax and reaching the e.g. 90% of max possible count
			if (ic->second >= (MaxPossCons* LCAfrac) ) {
				FullSetCnt = ic->second - unknwCnt; //reset full set cnt to the current best, to allow for lower levels to get to this
				goDeeper = true;
				ret->set(DL, ic->first);
				consens = ic->first;
#ifdef LCAdebg
				cout << " : " << consens ;
#endif
				break;
			}
		}

#ifdef LCAdebg
		cout << endl;
#endif
		if (!goDeeper) { break; }

		if (FullSetCnt == (MaxSetCnt - unknwCnt) ) {
			continue;
		}

		//while loop again, to exclude non consens assignments
		pos = -1;
		for (auto ix = TO.begin(); ix != TO.end(); ix++) {
			pos++;
			if ( (*ix)->get(DL) == consens) { continue; } //(*ix)->get(DL) == __unkwnTax ||
			DO[pos] = false;
		}
	}
	return ret;
}

double filterBlastPrimary(list<BlastRes*>& BR, options* opt, double& bestID) {
	
	if (!opt->BLfilter) {
		return (*BR.begin())->perID * (double)(*BR.begin())->alLen;
	}

	float minCov = opt->minCover;
	int maxL(0);
	
	for (auto it = BR.begin(); it != BR.end();it++) {
		if ((*it)->perID > bestID && (*it)->alLen >= maxL*0.95) { bestID = (*it)->perID; maxL = (*it)->alLen;}
		else if ((*it)->perID > bestID*0.9 && (*it)->alLen >= maxL*1.2) { bestID = (*it)->perID; maxL = (*it)->alLen; }
	}
	
	//filter parameters
	double lengthToleranceF(0.85f);
	double tolerance(1.5);
	if (opt->reportBestHit) {
		tolerance = 0.f; lengthToleranceF= 1.05f;
	}
	else if (bestID >= 100) { tolerance = 0.05f; }
	else if (bestID >= 99.5) { tolerance = 0.15f; }
	else if (bestID >= 99) { tolerance = 0.25f; }
	else if (bestID >= 98) { tolerance = 0.75f; }
	else if (bestID >= 97) { tolerance = 1.0f; }
	
	std::list<BlastRes*>::iterator i = BR.begin();

	//precalc maxL
	maxL *= lengthToleranceF;
	if (maxL < opt->minAliLen) {
		maxL = opt->minAliLen;
	}

	while (i != BR.end())	{

		/*
		if (((*i)->alLen * lengthToleranceF) < maxL) {
			int xxx = 1; cerr << "L";
		}
		if (((*i)->Qcoverage < minCov)) {
			int xxx = 1; cerr << "C";
		}
		if (((*i)->perID + tolerance) < bestID) {
			int xxx = 1; cerr << "I";
		}
		*/
		
		if ( ((*i)->perID + tolerance) < bestID ||
			((*i)->alLen ) < maxL || 
			(*i)->Qcoverage < minCov) {
			BR.erase(i++);
		} else {
			++i;
		}
	}
	//cerr << "\n";

	//filter done
	return bestID*(double)maxL;
}


list<TaxObj*> BlastToTax(list<BlastRes*>& BR, RefTax* RT, options* opt, float& consPerID) {
	list<TaxObj*> ret(0);
	vector<double> &thr = opt->idThr;


	int depth = RT->depth();
	bool anySpeciesCertain(false);
	for (auto it = BR.begin(); it != BR.end(); it++) {
		double curID((*it)->perID);
		consPerID += (float)curID;
		auto fnd = RT->find((*it)->Sbj);
		if (fnd != RT->end()) {
			TaxObj* F = new TaxObj(fnd->second);
			int maxD(0);
			//assign max depth based on % id to subject
			while (maxD < depth && thr[maxD] < curID ) { maxD++; }
			F->depth = maxD;
			ret.push_back(F);
			if (!F->speciesUncertain) { anySpeciesCertain = true; }
		}	else {
			cerr << "Could not find tax for Subject " << (*it)->Sbj << endl;
			exit(74);
		}
	}
	consPerID /= BR.size();

	//remove uncertain species, in case of enough good hits
	if (anySpeciesCertain && ret.size() > 1) {
		for (auto it = ret.begin(); it != ret.end(); it++) {
			(*it)->makeSpeciesUnknown();
		}
	}


	return ret;
}