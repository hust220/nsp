#include "Mat5.hpp"

BEGIN_JN

Mat5::Mat5() {
    init();
}

Mat5::Mat5(int a, int b, int c, int d, int e) {
    init();
    resize(a, b, c, d, e);
}

Mat5::Mat5(Mat5 &&m) {
    init();
    std::swap(size, m.size);
    std::swap(p, m.p);
}

Mat5 &Mat5::operator=(Mat5 &&m) {
    std::swap(size, m.size);
    std::swap(p, m.p);
	return *this;
}

void Mat5::init() {
    size = {0, 0, 0, 0, 0};
    p = NULL;
}

void Mat5::resize(int a, int b, int c, int d, int e) {
    int i, j, k, l;

    if (p != NULL) clear();
    size[0] = a;
    size[1] = b;
    size[2] = c;
    size[3] = d;
    size[4] = e;
    p = new float****[a];
    for (i = 0; i < a; i++) {
        p[i] = new float***[b];
        for (j = 0; j < b; j++) {
            p[i][j] = new float**[c];
            for (k = 0; k < c; k++) {
                p[i][j][k] = new float*[d];
                for (l = 0; l < d; l++) {
                    p[i][j][k][l] = new float[e];
                }
            }
        }
    }
}

void Mat5::clear() {
    if (p != NULL) {
        int i, j, k, l;
        for (i = 0; i < size[0]; i++) {
            for (j = 0; j < size[1]; j++) {
                for (k = 0; k < size[2]; k++) {
                    for (l = 0; l < size[3]; l++) {
                        delete [] p[i][j][k][l];
                    }
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

Mat5 Mat5::Zero(int a, int b, int c, int d, int e) {
    Mat5 m(a, b, c, d, e);
    m.each([](float &n, int a, int b, int c, int d, int e){
        n = 0;
    });
    return m;
}

Mat5 Mat5::Ones(int a, int b, int c, int d, int e) {
    Mat5 m(a, b, c, d, e);
    m.each([](float &n, int a, int b, int c, int d, int e){
        n = 1;
    });
    return m;
}

Mat5::~Mat5() {
    clear();
}

float &Mat5::operator()(int a, int b, int c, int d, int e) {
    return p[a][b][c][d][e];
}

const float &Mat5::operator()(int a, int b, int c, int d, int e) const {
    return p[a][b][c][d][e];
}

void Mat5::print() const {
    std::cout << "size: ";
    for (auto && i : size) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                for (int l = 0; l < size[3]; l++) {
                    for (int m = 0; m < size[4]; m++) {
                        std::cout << p[i][j][k][l][m] << ' ';
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
}

std::ostream &operator <<(std::ostream &out, const Mat5 &mat) {
	return out;
}

END_JN


