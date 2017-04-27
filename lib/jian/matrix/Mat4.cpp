#include "Mat4.hpp"

#include <iostream>

BEGIN_JN

Mat4::Mat4() {
    init();
}

Mat4::Mat4(int a, int b, int c, int d) {
    init();
    resize(a, b, c, d);
}

Mat4::Mat4(Mat4 &&m) {
    init();
    std::swap(size, m.size);
    std::swap(p, m.p);
}

Mat4 &Mat4::operator=(Mat4 &&m) {
    std::swap(size, m.size);
    std::swap(p, m.p);
	return *this;
}

void Mat4::init() {
    size = {0, 0, 0, 0};
    p = NULL;
}

void Mat4::resize(int a, int b, int c, int d) {
    int i, j, k;

    if (p != NULL) clear();
    size[0] = a;
    size[1] = b;
    size[2] = c;
    size[3] = d;
    p = new float***[a];
    for (i = 0; i < a; i++) {
        p[i] = new float**[b];
        for (j = 0; j < b; j++) {
            p[i][j] = new float*[c];
            for (k = 0; k < c; k++) {
                p[i][j][k] = new float[d];
            }
        }
    }
}

void Mat4::clear() {
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

Mat4 Mat4::Zero(int a, int b, int c, int d) {
    Mat4 m(a, b, c, d);
    m.each([](float &n, int a, int b, int c, int d){
        n = 0;
    });
    return m;
}

Mat4 Mat4::Ones(int a, int b, int c, int d) {
    Mat4 m(a, b, c, d);
    m.each([](float &n, int a, int b, int c, int d){
        n = 1;
    });
    return m;
}

Mat4::~Mat4() {
    clear();
}

float &Mat4::operator()(int a, int b, int c, int d) {
    return p[a][b][c][d];
}

const float &Mat4::operator()(int a, int b, int c, int d) const {
    return p[a][b][c][d];
}

void Mat4::print() const {
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                for (int l = 0; l < size[3]; l++) {
                    std::cout << p[i][j][k][l] << ' ';
                }
                std::cout << std::endl;
            }
        }
    }
}

END_JN


