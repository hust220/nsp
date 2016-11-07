#include "Mat3.hpp"

#include <iostream>

namespace jian {

Mat3::Mat3() {
    init();
}

Mat3::Mat3(int a, int b, int c) {
    init();
    resize(a, b, c);
}

Mat3::Mat3(Mat3 &&m) {
    init();
    std::swap(size, m.size);
    std::swap(p, m.p);
}

Mat3 &Mat3::operator=(Mat3 &&m) {
    std::swap(size, m.size);
    std::swap(p, m.p);
	return *this;
}

void Mat3::init() {
    size = {0, 0, 0};
    p = NULL;
}

void Mat3::resize(int a, int b, int c) {
    int i, j;

    if (p != NULL) clear();
    size[0] = a;
    size[1] = b;
    size[2] = c;
    p = new float**[a];
    for (i = 0; i < a; i++) {
        p[i] = new float*[b];
        for (j = 0; j < b; j++) {
            p[i][j] = new float[c];
        }
    }
}

void Mat3::clear() {
    if (p != NULL) {
        int i, j;
        for (i = 0; i < size[0]; i++) {
            for (j = 0; j < size[1]; j++) {
                delete [] p[i][j];
            }
            delete [] p[i];
        }
        delete [] p;
        p = NULL;
    }
}

Mat3 Mat3::Zero(int a, int b, int c) {
    Mat3 m(a, b, c);
    m.each([](float &n, int a, int b, int c){
        n = 0;
    });
    return m;
}

Mat3 Mat3::Ones(int a, int b, int c) {
    Mat3 m(a, b, c);
    m.each([](float &n, int a, int b, int c){
        n = 1;
    });
    return m;
}

Mat3 Mat3::Constant(int a, int b, int c, float v) {
    Mat3 m(a, b, c);
    m.each([&v](float &n, int a, int b, int c){
        n = v;
    });
    return m;
}

Mat3::~Mat3() {
    clear();
}

float &Mat3::operator()(int a, int b, int c) {
    return p[a][b][c];
}

const float &Mat3::operator()(int a, int b, int c) const {
    return p[a][b][c];
}

void Mat3::print() const {
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                std::cout << p[i][j][k] << ' ';
            }
            std::cout << std::endl;
        }
    }
}

} // namespace jian


