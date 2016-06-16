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
    std::vector<std::pair<char, char>> paired_keys {
        {'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}, 
        {'A', 'a'}, {'B', 'b'}, {'C', 'c'}, {'D', 'd'}, 
        {'E', 'e'}, {'F', 'f'}, {'G', 'g'}, {'H', 'h'},
        {'I', 'i'}, {'J', 'j'}, {'K', 'k'}, {'L', 'l'},
        {'M', 'm'}, {'N', 'n'}, {'O', 'o'}
    };
    std::vector<char> unpaired_keys { '.', ':' };
    std::vector<char> break_keys { '&', '-' };

    static NucSS &instance() {
        static NucSS nuc_ss;
        return nuc_ss;
    }

    static std::map<char, char> get_map_keys() {
        std::map<char, char> map_keys;
        for (auto &&pair: instance().paired_keys) {
            map_keys[pair.first] = pair.second;
            map_keys[pair.second] = pair.first;
        }
        return map_keys;
    }

    static std::map<char, int> get_pos_keys() {
        std::map<char, int> pos_keys;
        int index = 1;
        for (auto &&pair: instance().paired_keys) {
            pos_keys[pair.first] = index;
            pos_keys[pair.second] = -index;
            index++;
        }
        return pos_keys;
    }

    static bool check_ss(const std::string &ss) {
        static std::map<char, char> map_keys = get_map_keys();
        static std::map<char, int> pos_keys = get_pos_keys();

        for (auto &&s: ss) {
            if (not std::count_if(instance().paired_keys.begin(), instance().paired_keys.end(), [&](const std::pair<char, char> &pair) {
                    return pair.first == s || pair.second == s;
                }) and not std::count(instance().unpaired_keys.begin(), instance().unpaired_keys.end(), s) 
                   and not std::count(instance().break_keys.begin(), instance().break_keys.end(), s)) {
                return false;
            }
        }
        std::map<char, int> map;
        for (auto &&s: ss) {
            if (pos_keys.count(s)) {
                if (pos_keys[s] > 0) {
                    if (map.count(s)) map[s]++; else map[s] = 1;
                } else {
                    char c = map_keys[s];
                    if (map.count(c)) {
                        map[c]--;    
                        if (map[c] < 0) return false;
                    } else {
                        return false;    
                    }
                }
            }
        }

        for (auto &&pair: map) if (pair.second != 0) return false;

        return true;
    }

    static int len_ss(const std::string &ss) {
        return std::count_if(ss.begin(), ss.end(), [&](const char &c){
            return std::count_if(instance().paired_keys.begin(), instance().paired_keys.end(), [&](const std::pair<char, char> &pair){
                return pair.first == c || pair.second == c;
            }) || std::count(instance().unpaired_keys.begin(), instance().unpaired_keys.end(), c);
        });
    }

    static std::string pure_ss(const std::string &ss) {
        std::string p_ss;
        std::copy_if(ss.begin(), ss.end(), std::back_inserter(p_ss), [](const char &c){return c != '&';});
        return p_ss;
    }
        
    static std::string lower_ss(const std::string &ss, int n = 2) {
        static std::map<char, int> pos_keys = get_pos_keys();
        std::string l_ss;
        std::transform(ss.begin(), ss.end(), std::back_inserter(l_ss), [&](const char &c){
            if (pos_keys.count(c) && (pos_keys[c] > n || pos_keys[c] < -n)) {
                return '.';    
            } else {
                return c;
            }
        });
        return l_ss;
    }

    static std::string hinge_ss(const std::string &ss) {
        std::string h_ss; EACH(i, ss, if (i == '(' || i == ')') h_ss += i); return h_ss;
    }

    static bool seq_match_ss(const std::string &seq, const std::string &ss) {
        int len_ss = 0; EACH(i, ss, if (i != '&') len_ss++);
        return seq.size() == len_ss;
    }

};

} // namespace jian

