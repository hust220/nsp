#ifndef JIAN_MLIB_H
#define JIAN_MLIB_H

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

namespace jian {

void tokenize(const string &, vector<string> &, const string & = " ");
string upper(string);
string lower(string);
int die(string str);
std::string env(std::string);

} /// namespace jian

#endif





