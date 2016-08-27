#pragma once

#include <array>

namespace jian {

class Mat4 {
public:
    std::array<int, 4> size;
    double ****p;

    Mat4() {
        init();
    }

    Mat4(int a, int b, int c, int d) {
        init();
        resize(a, b, c, d);
    }

    Mat4(Mat4 &&m) {
        init();
        std::swap(p, m.p);
    }

    Mat4 &operator=(Mat4 &&m) {
        std::swap(p, m.p);
    }

    void init() {
        size = {0, 0, 0, 0};
        p = NULL;
    }

    void resize(int a, int b, int c, int d) {
        int i, j, k;

        if (p != NULL) clear();
        size[0] = a;
        size[1] = b;
        size[2] = c;
        size[3] = d;
        p = new double***[a];
        for (i = 0; i < a; i++) {
            p[i] = new double**[b];
            for (j = 0; j < b; j++) {
                p[i][j] = new double*[c];
                for (k = 0; k < c; k++) {
                    p[i][j][k] = new double[d];
                }
            }
        }
    }

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

    void clear() {
        if (p != NULL) {
            int i, j, k;
            for (i = 0; i < size[0]; i++) {
                for (j = 0; j < size[1]; j++) {
                    for (k = 0; k < size[2]; k++) {
                        delete [] p[i][j][k];
                    }
                    delete [] p[i][j];
                }
                delete [] p[i];
            }
            delete [] p;
            p = NULL;
        }
    }

    static Mat4 Zero(int a, int b, int c, int d) {
        Mat4 m(a, b, c, d);
        m.each([](double &n, int a, int b, int c, int d){
            n = 0;
        });
        return m;
    }

    ~Mat4() {
        clear();
    }

    double &operator()(int a, int b, int c, int d) {
        return p[a][b][c][d];
    }
};

} // namespace jian


