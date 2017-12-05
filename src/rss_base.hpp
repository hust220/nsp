#pragma once

#include "string.hpp"
#include "pp.hpp"
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <algorithm>

BEGIN_JN

class NASS {
public:
	using Key = char;
	using Keys = std::vector<Key>;
	using PairedKey = std::pair<Key, Key>;
	using PairedKeys = std::vector<PairedKey>;

    PairedKeys paired_keys;
    Keys unpaired_keys;
    Keys break_keys;

	static std::vector<int> get_bps(const Str &ss);
    static const NASS &instance();
	static bool is_char_ss(char c);
	static bool check_ss(const Str &ss, Str &info_errors);
	static bool check_ss(const Str &ss);
    static int len_ss(const Str &ss);
    static Str pure_ss(const Str &ss);
    static Str lower_ss(const Str &ss, int n = 2);
    static Str hinge_ss(const Str &ss);
    static bool seq_match_ss(const Str &seq, const Str &ss);

private:
	NASS();
};

END_JN

