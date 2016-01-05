#ifndef JIAN_GEOM_ROTATE_H
#define JIAN_GEOM_ROTATE_H

#include "../fpl.h"
#include "../util/util.h"

namespace jian {
namespace geom {

template<typename T, typename Mat> 
void rotate(T &&t, Mat &&mat) {
    double x = t[0] * mat::ref(mat, 0, 0) + t[1] * mat::ref(mat, 1, 0) + t[2] * mat::ref(mat, 2, 0);
    double y = t[0] * mat::ref(mat, 0, 1) + t[1] * mat::ref(mat, 1, 1) + t[2] * mat::ref(mat, 2, 1);
    double z = t[0] * mat::ref(mat, 0, 2) + t[1] * mat::ref(mat, 1, 2) + t[2] * mat::ref(mat, 2, 2);
    t[0] = x; t[1] = y; t[2] = z;
}

template<class N, class L> 
N times(const MatrixXf &mat, const L &l) {
    N n = l;
    for (int i = 0; i < mat.rows(); i++) {
        n[i] = 0;
        for (int j = 0; j < mat.cols(); j++) {
            n[i] += mat(i, j);
        }
    }
    return n;
} /// times

template<class C1 = double, class C2 = double> 
MatrixXf x_rot_mat(C1 c, C2 s) {
    MatrixXf rot_mat(3, 3);
    rot_mat << 1, 0, 0,
               0, c, s,
               0,-s, c;
    return rot_mat;   
} /// x_rot_mat

template<class C1 = double, class C2 = double> 
MatrixXf y_rot_mat(C1 c, C2 s) {
    MatrixXf rot_mat(3, 3);
    rot_mat << c, 0,-s,
               0, 1, 0,
               s, 0, c;
    return rot_mat;   
} /// y_rot_mat

template<class C1 = double, class C2 = double> 
MatrixXf z_rot_mat(C1 c, C2 s) {
    MatrixXf rot_mat(3, 3);
    rot_mat << c, s, 0,
              -s, c, 0,
               0, 0, 1;
    return rot_mat;   
} /// z_rot_mat

template<typename L> 
void rotate(MatrixXf &src, const L &begin, const L &end, const double &angle) {
    L l = end - begin;
    for (int i = 0; i < src.rows(); i++) {
        for (int j = 0; j < src.cols(); j++) {
            src(i, j) -= begin[j];
        }
    }

    // rotate with z to put l on y-z plane
    double r, c, s, x_, y_, z_;
    r = sqrt(l[0] * l[0] + l[1] * l[1]);
    if (r != 0) {
        c = l[1] / r;
        s = l[0] / r;
        auto mat = z_rot_mat(c, s);
        src = src * mat;
    }

    // rotate with x to put l on z
    r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    if (r != 0) {
        c = l[2] / r;
        s = sqrt(l[0] * l[0] + l[1] * l[1]) / r;
        auto mat = x_rot_mat(c, s);
        src = src * mat;
    }

    // rotate with z
    c = cos(angle);
    s = sin(angle);
    auto mat = z_rot_mat(c, s);
    src = src * mat;

    // rotate with x
    r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    if (r != 0) {
        c = l[2] / r;
        s = -sqrt(l[0] * l[0] + l[1] * l[1]) / r;
        auto mat = x_rot_mat(c, s);
        src = src * mat;
    }

    // rotate with z
    r = sqrt(l[0] * l[0] + l[1] * l[1]);
    if (r != 0) {
        c = l[1] / r;
        s = -l[0] / r;
        auto mat = z_rot_mat(c, s);
        src = src * mat;
    }

    // Move src to the initial position
    for (int i = 0; i < src.rows(); i++) {
        for (int j = 0; j < src.cols(); j++) {
            src(i, j) += begin[j];
        }
    }
} // rotate

} // namespace geometry
} // namespace jian

#endif

