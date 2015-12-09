#ifndef JIAN_NUC2D_UTIL
#define JIAN_NUC2D_UTIL

#include <util/std.h>

namespace jian {
namespace nuc2d {

static std::vector<std::pair<char, char>> paired_keys = {
    {'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}, 
    {'A', 'a'}, {'B', 'b'}, {'C', 'c'}, {'D', 'd'}, 
    {'E', 'e'}, {'F', 'f'}, {'G', 'g'}, {'H', 'h'},
    {'I', 'i'}, {'J', 'j'}, {'K', 'k'}, {'L', 'l'},
    {'M', 'm'}, {'N', 'n'}, {'O', 'o'}
};

static std::vector<char> unpaired_keys = {
    '.', ':' 
};

static std::vector<char> break_keys = {
    '&', '-'
};

inline std::map<char, char> get_map_keys() {
    std::map<char, char> map_keys;
    for (auto &&pair: paired_keys) {
        map_keys[pair.first] = pair.second;
        map_keys[pair.second] = pair.first;
    }
    return map_keys;
}

inline std::map<char, int> get_pos_keys() {
    std::map<char, int> pos_keys;
    int index = 1;
    for (auto &&pair: paired_keys) {
        pos_keys[pair.first] = index;
        pos_keys[pair.second] = -index;
        index++;
    }
    return pos_keys;
}

template<typename T> bool check_ss(T &&ss) {
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

template<typename T> int len_ss(T &&ss) {
    return std::count_if(ss.begin(), ss.end(), [&](const char &c){
        return std::count_if(paired_keys.begin(), paired_keys.end(), [&](const std::pair<char, char> &pair){
            return pair.first == c || pair.second == c;
        }) || std::count(unpaired_keys.begin(), unpaired_keys.end(), c);
    });
}

template<typename T> std::string pure_ss(T &&ss) {
    std::string p_ss;
    std::copy_if(ss.begin(), ss.end(), std::back_inserter(p_ss), [](const char &c){return c != '&';});
    return p_ss;
}
    
template<typename T> std::string lower_ss(T &&ss, int n = 2) {
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

} // namespace nuc2d
} // namespace jian

#endif

