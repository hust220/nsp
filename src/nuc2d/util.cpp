#include "util.h"

namespace jian {
namespace nuc2d {

std::map<char, char> get_map_keys() {
    std::map<char, char> map_keys;
    for (auto &&pair: paired_keys) {
        map_keys[pair.first] = pair.second;
        map_keys[pair.second] = pair.first;
    }
    return map_keys;
}

std::map<char, int> get_pos_keys() {
    std::map<char, int> pos_keys;
    int index = 1;
    for (auto &&pair: paired_keys) {
        pos_keys[pair.first] = index;
        pos_keys[pair.second] = -index;
        index++;
    }
    return pos_keys;
}

} // namespace nuc2d
} // namespace jian


