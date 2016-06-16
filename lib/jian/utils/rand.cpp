#include <random>

namespace jian {

static std::mt19937 _rand_engine{11};
static std::uniform_real_distribution<double> _unif_distr{0, 1};

double rand() {
    return _unif_distr(_rand_engine);
}

void seed(const double &t) {
    _rand_engine.seed(t);
}

} // namespace jian

