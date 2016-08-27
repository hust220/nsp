#include <mutex>
#include <iomanip>
#include <set>
#include "Chain.hpp"
#include "molstream.hpp"
#include "../utils/Debug.hpp"
#include "../utils/file.hpp"

namespace jian {

static std::mutex mt;

int num_atoms(const Chain &chain) {
    int n = 0;
    for (auto && res : chain) {
        for (auto && atom : res) {
            n++;
        }
    }
    return n;
}

} // namespace jian

