#include <algorithm>
#include <numeric>
#include "DihBound.hpp"
#include <functional>
#include "../pp.hpp"
#include "../utils/traits.hpp"

BEGIN_JN

namespace DihBoundImpl {

	int hash::operator ()(const key_t &v) const {
		std::hash<int> h;
		//return JN_FOLD(_1 ^ (h(_2) << 1), (int, 0), v);
		return std::accumulate(v.begin(), v.end(), 0, [&h](int n, auto && m) {return n ^ (h(m) << 1); });
	}

	bool equal_to::operator ()(const key_t &vec1, const key_t &vec2) const {
		if (vec1.size() != vec2.size()) {
			return false;
		}
		for (int i = 0; i < vec1.size(); i++) {
			if (vec1[i] != vec2[i]) {
				return false;
			}
		}
		return true;
	}

} // namespace DihBoundImpl

std::ostream &operator <<(std::ostream &out, const DihBound &dih_bound) {
	for (auto && i : dih_bound) {
		for (auto && n : i.first) {
			std::cout << ':' << i.second << std::endl;
		}
	}
	return out;
}

END_JN

