#pragma once

#include <array>

namespace jian {

template<int N>
class Mat {
public:
    using shape_t = std::array<size_t, N>;

    shape_t shape;
    shape_t size;
    float *p;

    Mat() {
        init();
    }

    Mat(const shape_t &s) {
        init();
        resize(s);
    }

    Mat(Mat &&m) {
        init();
        std::swap(p, m.p);
    }

    Mat &operator=(Mat &&m) {
        std::swap(p, m.p);
    }

    void init() {
        for (auto && i : shape) i = 0;
        p = NULL;
    }

    void resize(const shape_t &s) {
        if (p != NULL) clear();
        shape = s;
        size[N-1] = shape[N-1];
        for (size_t i = N-2; i >= 0; i--) {
            size[i] = shape[i] * size[i+1];
        }
        p = new float[size[0]];
    }

    void clear() {
        if (p != NULL) {
            delete [] p;
            p = NULL;
        }
    }

    static Mat Zero(const shape_t &s) {
        Mat m(s);
        for (size_t i = 0; i < size[0]; i++) (*p)[i] = 0;
        return m;
    }

    ~Mat() {
        clear();
    }

    float &operator()(int a, int b, int c, int d) {
        return p[a][b][c][d];
    }
};

} // namespace jian


