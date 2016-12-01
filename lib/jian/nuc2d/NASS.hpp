#pragma once

#include "../pp.hpp"
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <algorithm>

namespace jian {

class NASS {
public:
	using Key = char;
	using Keys = std::vector<Key>;
	using PairedKey = std::pair<Key, Key>;
	using PairedKeys = std::vector<PairedKey>;

    PairedKeys paired_keys;
    Keys unpaired_keys;
    Keys break_keys;

	static std::vector<int> get_bps(const str_t &ss);
    static const NASS &instance();
	static bool is_char_ss(char c);
	static bool check_ss(const str_t &ss, str_t &info_errors);
	static bool check_ss(const str_t &ss);
    static int len_ss(const str_t &ss);
    static str_t pure_ss(const str_t &ss);
    static str_t lower_ss(const str_t &ss, int n = 2);
    static str_t hinge_ss(const str_t &ss);
    static bool seq_match_ss(const str_t &seq, const str_t &ss);

private:
	NASS();
};

} // namespace jian

