#pragma once

#include <array>
#include "../utils/traits.hpp"

BEGIN_JN

class Mat4 {
public:
    std::array<int, 4> size;
    float ****p;

    Mat4();
    Mat4(int a, int b, int c, int d);
    Mat4(Mat4 &&m);
    Mat4 &operator=(Mat4 &&m);
    void init();
    void resize(int a, int b, int c, int d);

    template<typename F>
    void each(F && f) {
        for (int i = 0; i < size[0]; i++) {
            for (int j = 0; j < size[1]; j++) {
                for (int k = 0; k < size[2]; k++) {
                    for (int l = 0; l < size[3]; l++) {
                        f(p[i][j][k][l], i, j, k, l);
                    }
                }
            }
        }
    }

    void clear();
    static Mat4 Zero(int a, int b, int c, int d);
    static Mat4 Ones(int a, int b, int c, int d);
    ~Mat4();
    float &operator()(int a, int b, int c, int d);
    const float &operator()(int a, int b, int c, int d) const;
    void print() const;
};

END_JN


