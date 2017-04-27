#pragma once

#include <array>
#include <iostream>
#include "../utils/traits.hpp"

BEGIN_JN

class Mat5 {
public:
    std::array<int, 5> size;
    float *****p;

    Mat5();
    Mat5(int a, int b, int c, int d, int e);
    Mat5(Mat5 &&m);
    Mat5 &operator=(Mat5 &&m);
    void init();
    void resize(int a, int b, int c, int d, int e);

    template<typename F>
    void each(F && f) {
        for (int i = 0; i < size[0]; i++) {
            for (int j = 0; j < size[1]; j++) {
                for (int k = 0; k < size[2]; k++) {
                    for (int l = 0; l < size[3]; l++) {
                        for (int m = 0; m < size[4]; m++) {
                            f(p[i][j][k][l][m], i, j, k, l, m);
                        }
                    }
                }
            }
        }
    }

    void clear();
    static Mat5 Zero(int a, int b, int c, int d, int e);
    static Mat5 Ones(int a, int b, int c, int d, int e);
    ~Mat5();
    float &operator()(int a, int b, int c, int d, int e);
    const float &operator()(int a, int b, int c, int d, int e) const;
    void print() const;
};

std::ostream &operator <<(std::ostream &out, const Mat5 &mat);

END_JN


