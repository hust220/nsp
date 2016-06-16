#include "DihBound.hpp"
#include <functional>
#include "../pp.hpp"

namespace jian {
namespace DihBoundImpl {

unsigned int hash::operator ()(const key_t &v) const {
    std::hash<int> h; 
    return JN_FOLD(_1 ^ (h(_2) << 1), (unsigned int, 0), v); 
}

bool equal_to::operator ()(const key_t &vec1, const key_t &vec2) const{
    return vec1.size() == vec2.size() && ! JN_EXISTS(_1 != _2, vec1, vec2);
}

} // namespace DihBoundImpl

std::ostream &operator <<(std::ostream &out, const DihBound &dih_bound) {
    EACH(i, dih_bound, 
        EACH(n, i.first, 
            std::cout << n << ' ';
        ); 
        std::cout << ':' << i.second << std::endl;
    );
}

} // namespace jian

