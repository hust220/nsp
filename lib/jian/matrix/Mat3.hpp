#pragma once

#include <array>

namespace jian {

class Mat3 {
public:
    std::array<int, 3> size;
    float ***p;

    Mat3();
    Mat3(int a, int b, int c);
    Mat3(Mat3 &&m);
    Mat3 &operator=(Mat3 &&m);
    void init();
    void resize(int a, int b, int c);

    template<typename F>
    void each(F && f) {
        for (int i = 0; i < size[0]; i++) {
            for (int j = 0; j < size[1]; j++) {
                for (int k = 0; k < size[2]; k++) {
                    f(p[i][j][k], i, j, k);
                }
            }
        }
    }

    void clear();
    static Mat3 Zero(int a, int b, int c);
    static Mat3 Ones(int a, int b, int c);
    static Mat3 Constant(int a, int b, int c, float v);
    ~Mat3();
    float &operator()(int a, int b, int c);
    const float &operator()(int a, int b, int c) const;
    void print() const;
};

} // namespace jian


