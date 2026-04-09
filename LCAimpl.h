#pragma once
#include "RefTax.h"

//routine that calls various subfunctions like blast filter, LCA etc
TaxObj* LCA(vector<BlastRes>&, RefTax*, options*);

//filtering of the blast results
double filterBlastPrimary(vector<BlastRes>&, options* opt,double&);
//routine that performs actual LCA matching etc
TaxObj* LCAcore(const vector<TaxObj*>&, bool &hitRd, double LCAfrac=0.9f, int tdepth=__default_depth);


vector<TaxObj*> BlastToTax(const vector<BlastRes>& BR, RefTax* RT, options*, float& consPerID);
