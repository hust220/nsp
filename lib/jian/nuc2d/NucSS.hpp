#pragma once

#include "../pp.hpp"
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <algorithm>

namespace jian {

class NucSS {
public:
    std::vector<std::pair<char, char>> paired_keys;
    std::vector<char> unpaired_keys;
    std::vector<char> break_keys;

    NucSS();
    static NucSS &instance();
	static bool is_char_ss(char c);
	static bool check_ss(const std::string &ss, std::string &info_errors);
	static bool check_ss(const std::string &ss);
    static int len_ss(const std::string &ss);
    static std::string pure_ss(const std::string &ss);
    static std::string lower_ss(const std::string &ss, int n = 2);
    static std::string hinge_ss(const std::string &ss);
    static bool seq_match_ss(const std::string &seq, const std::string &ss);
};

} // namespace jian

