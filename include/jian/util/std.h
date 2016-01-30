#ifndef JIAN_STD_H
#define JIAN_STD_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cassert>
#include <regex>

// STL
#include <bitset>
#include <vector>
#include <array>
#include <queue>
#include <numeric>
#include <ctime>
#include <list>
#include <set>
#include <deque>
#include <unordered_map>
#include <map>
#include <memory>
#include <iterator>
#include <utility>
#include <algorithm>
#include <functional>
#include <type_traits>
using namespace std;

// boost
#include <boost/timer.hpp>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/preprocessor.hpp>

// constant
#define MASS_H 1.00794
#define MASS_C 12.0107
#define MASS_O 15.9994
#define MASS_N 14.0067
#define MASS_P 30.973762

namespace jian {
const static double PI = 3.1415927;
} // namespace jian

// Macro
#define CLASS_NAME(T) #T

#endif

