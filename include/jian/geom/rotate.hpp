#pragma once

#include "../etl.hpp"

namespace jian {
namespace geom {

using val_t = double;

template<typename T, typename U> 
void translate(T &t, const U &u) {
    for (int i = 0; i < 3; i++) t[i] += u[i];
}

template<typename T, typename Mat> 
void rotate(T &t, const Mat &mat) {
    val_t x = t[0] * ref(mat, 0, 0) + t[1] * ref(mat, 1, 0) + t[2] * ref(mat, 2, 0);
    val_t y = t[0] * ref(mat, 0, 1) + t[1] * ref(mat, 1, 1) + t[2] * ref(mat, 2, 1);
    val_t z = t[0] * ref(mat, 0, 2) + t[1] * ref(mat, 1, 2) + t[2] * ref(mat, 2, 2);
    t[0] = x; t[1] = y; t[2] = z;
}

template<typename T, typename U, typename Mat> 
void rotate(T &t, const U &origin, const Mat &mat) {
    for (int i = 0; i < 3; i++) t[i] -= origin[i];
    val_t x = t[0] * ref(mat, 0, 0) + t[1] * ref(mat, 1, 0) + t[2] * ref(mat, 2, 0);
    val_t y = t[0] * ref(mat, 0, 1) + t[1] * ref(mat, 1, 1) + t[2] * ref(mat, 2, 1);
    val_t z = t[0] * ref(mat, 0, 2) + t[1] * ref(mat, 1, 2) + t[2] * ref(mat, 2, 2);
    t[0] = x; t[1] = y; t[2] = z;
    for (int i = 0; i < 3; i++) t[i] += origin[i];
}

template<typename Mat = MatrixXd, typename T>
Mat rot_mat(int i, T &&v) {
    Mat mat(3, 3);
    double c = std::cos(v);
    double s = std::sin(v);
    if (i == 0) {
        mat << 1, 0, 0,
               0, c, s,
               0,-s, c;
    } else if (i == 1) {
        mat << c, 0,-s,
               0, 1, 0,
               s, 0, c;
    } else if (i == 2) {
        mat << c, s, 0,
              -s, c, 0,
               0, 0, 1;
    }
    return mat;
}

template<typename Mat = MatrixXd, class C1 = val_t, class C2 = val_t> 
Mat x_rot_mat(C1 c, C2 s) {
    Mat rot_mat(3, 3);
    rot_mat << 1, 0, 0,
               0, c, s,
               0,-s, c;
    return rot_mat;   
} // x_rot_mat

template<typename Mat = MatrixXd, class C1 = val_t, class C2 = val_t> 
Mat y_rot_mat(C1 c, C2 s) {
    Mat rot_mat(3, 3);
    rot_mat << c, 0,-s,
               0, 1, 0,
               s, 0, c;
    return rot_mat;   
} // y_rot_mat

template<typename Mat = MatrixXd, class C1 = val_t, class C2 = val_t> 
Mat z_rot_mat(C1 c, C2 s) {
    Mat rot_mat(3, 3);
    rot_mat << c, s, 0,
              -s, c, 0,
               0, 0, 1;
    return rot_mat;   
} // z_rot_mat

template<typename Mat, typename L> 
void rotate_angle_xyz(Mat &src, const L &begin, const L &end, const val_t &angle) {
    L l = end - begin;
    for (int i = 0; i < src.rows(); i++) for (int j = 0; j < src.cols(); j++) src(i, j) -= begin[j];

    // rotate with z to put l on y-z plane
    val_t r, c, s;
    r = sqrt(l[0] * l[0] + l[1] * l[1]);
    if (r != 0) src = src * z_rot_mat(l[1]/r, l[0]/r);

    // rotate with x to put l on z
    r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    if (r != 0) src = src * x_rot_mat(l[2]/r, std::sqrt(l[0]*l[0]+l[1]*l[1])/r);

    // rotate with z
    src = src * z_rot_mat(std::cos(angle), std::sin(angle));

    // rotate with x
    r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    if (r != 0) {
        c = l[2] / r;
        s = -sqrt(l[0] * l[0] + l[1] * l[1]) / r;
        src = src * x_rot_mat(c, s);
    }

    // rotate with z
    r = sqrt(l[0] * l[0] + l[1] * l[1]);
    if (r != 0) { c = l[1] / r; s = -l[0] / r; src = src * z_rot_mat(c, s); }

    // Move src to the initial position
    for (int i = 0; i < src.rows(); i++) for (int j = 0; j < src.cols(); j++) {
        src(i, j) += begin[j];
    }
} // rotate

} // namespace geometry
} // namespace jian

