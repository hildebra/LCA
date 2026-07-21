//libload.h
//contains constants and the std libs being loaded
#pragma once
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <iterator>
#include <map>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <cassert>
#include <unordered_map>
#include <unordered_set>



#if !defined(LCA_NO_GZIP) && !((defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__))
#define _gzipread
#endif

#ifdef _gzipread
#include "gzstream.h"
#endif

using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;


//some static vars..
enum { __default_depth = 7 };
inline const string __unkwnTax = "";
inline const string __unkwnTaxWR = "?";
inline const string __defaultTaxSep = "\t";


