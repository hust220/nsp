#pragma once

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include "../utils/traits.hpp"

BEGIN_JN

namespace DihBoundImpl {

	using key_t = std::vector<int>;

	using val_t = double;

	struct hash { int operator ()(const key_t &v) const; };

	struct equal_to { bool operator ()(const key_t &vec1, const key_t &vec2) const; };

}

using DihBound = std::unordered_map<DihBoundImpl::key_t, DihBoundImpl::val_t, DihBoundImpl::hash, DihBoundImpl::equal_to>;

std::ostream &operator <<(std::ostream &out, const DihBound &dih_bound);

END_JN

