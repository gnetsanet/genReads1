#ifndef MISCFUNC_H
#define MISCFUNC_H

#include "StringRef.h"

#include <vector>
#include <string>
#include <tr1/unordered_map>

using namespace std;

typedef tr1::unordered_map<string, unsigned long int> countMap;

vector<StringRef> split3( string const& str, char delimiter );

void defaultHumanChr( countMap &refLens );


#endif