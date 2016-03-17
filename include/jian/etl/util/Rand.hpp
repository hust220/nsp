#ifndef JIAN_UTIL_RAND
#define JIAN_UTIL_RAND

#include "../std.hpp"

namespace jian {

class Rand {
public:
    std::mt19937 _rand_engine{11};
    std::uniform_real_distribution<double> _unif_distr{0, 1};

    double rand() {
        return _unif_distr(_rand_engine);
    }

    template<typename T>
    void seed(T &&t) {
        _rand_engine.seed(t);
    }

};

} // namespace jian




#endif

