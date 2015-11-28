#ifndef JIAN_NUC2D_UTIL
#define JIAN_NUC2D_UTIL

#include <util/std.h>

namespace jian {
namespace nuc2d {

static std::map<char, char> paired_keys = {
    {'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}, 
    {'\\', '/'}, {'A', 'a'}, {'B', 'b'}, {'C', 'c'},
    {'D', 'd'}, {'E', 'e'}, {'F', 'f'}, {'G', 'g'}
};

static std::vector<char> unpaired_keys = {
    '.', ':' 
};

static std::vector<char> break_keys = {
    '&', '-'
};

std::map<char, char> get_map_keys();
std::map<char, int> get_pos_keys();

template<typename SS> bool check_ss(SS &&ss) {
    static std::map<char, char> map_keys = get_map_keys();
    static std::map<char, int> pos_keys = get_pos_keys();

    for (auto &&s: ss) {
        if (!std::count_if(paired_keys.begin(), paired_keys.end(), [&](const std::pair<char, char> &pair) {
                return pair.first == s || pair.second == s;
            }) && !std::count(unpaired_keys.begin(), unpaired_keys.end(), s) 
               && !std::count(break_keys.begin(), break_keys.end(), s)) {
            return false;
        }
    }

    std::map<char, int> map;
    for (auto &&s: ss) {
        if (pos_keys.count(s)) {
            if (pos_keys[s] > 0) {
                if (map.count(s)) {
                    map[s]++;    
                } else {
                    map[s] = 1;    
                }
            } else {
                char c = map_keys[s];
                if (map.count(c)) {
                    map[c]--;    
                    if (map[c] < 0) {
                        return false;    
                    } else {
                        // pass
                    }
                } else {
                    return false;    
                }
            }
        } else {
            // pass
        }
    }

    for (auto &&pair: map) {
        if (pair.second != 0) {
            return false;
        } else {
            // pass
        }
    }

    return true;
}

template<typename SS> int len_ss(SS &&ss) {
    return std::count_if(ss.begin(), ss.end(), [&](const char &c){
        return std::count_if(paired_keys.begin(), paired_keys.end(), [&](const std::pair<char, char> &pair){
            return pair.first == c || pair.second == c;
        }) || std::count(unpaired_keys.begin(), unpaired_keys.end(), c);
    });
}
    
} // namespace nuc2d
} // namespace jian

#endif

