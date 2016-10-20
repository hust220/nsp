#pragma once

#include <random>
#include <vector>

namespace jian {
namespace dg {

class Move {
public:
    std::random_device rd;
    std::uniform_real_distribution<double> distr;

    int _index;
    std::vector<double> _shift;

    Move() : distr(0, 1), _shift(3, 0) {}

    auto rand_shift() -> decltype(_shift) {
        return decltype(_shift) {distr(rd)-0.5,distr(rd)-0.5,distr(rd)-0.5};
    }

    int pick(int n) {
        _index = int(distr(rd) * n);
        return _index;
    }

    template<typename Mat> void move(Mat &&mat) {
        _shift = rand_shift();
        for (int i = 0; i < 3; i++) mat(_index, i) += _shift[i];
    }

    template<typename Mat> void back(Mat &&mat) {
        for (int i = 0; i < 3; i++) mat(_index, i) -= _shift[i];
    }
};


} // namespace dg
} // namespace jian


