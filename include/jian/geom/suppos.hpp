#pragma once

#include "../etl.hpp"
#include "rotate.hpp"

namespace jian {
namespace geom {

#define INIT_SUPPOS(sp) \
    auto &&sp##C1 = -(sp.c1); auto &&sp##ROT = sp.rot; auto &&sp##C2 = sp.c2;

#define APPLY_SUPPOS(m, sp) \
    geom::translate(m, sp##C1); geom::rotate(m, sp##ROT); geom::translate(m, sp##C2);

template<typename T, typename U>
void apply_suppos(T &&m, U &&sp) {
    auto &&c1 = -(sp.c1); auto &&rot = sp.rot; auto &&c2 = sp.c2;
    translate(m, c1); rotate(m, rot); translate(m, c2);
}

template<typename T, typename U>
auto suppos(const T &m, const U &n) {
    using val_t = double;
    using Vec = VectorXd;
    using Mat = MatrixXd;
    using Result = struct {Mat rot; Vec c1; Vec c2; val_t rmsd;};

    if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) {
        throw "jian::geom::suppos error! ("s + JN_STR(m.rows()) + ' ' + JN_STR(m.cols()) + ") (" + JN_STR(n.rows()) + ' ' + JN_STR(n.cols()) + ")\n";
    }
    int len = m.rows(); Mat x = m, y = n;
    Vec c1 = Vec::Zero(3), c2 = Vec::Zero(3);
    for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) { c1[j] += x(i, j); c2[j] += y(i, j); }
    for (int i = 0; i < 3; i++) { c1[i] = c1[i] / len; c2[i] = c2[i] / len; }
    for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) { x(i, j) -= c1[j]; y(i, j) -= c2[j]; }

    Mat g = x.transpose() * y;
    JacobiSVD<Mat> svd(g, ComputeFullU|ComputeFullV);
    Mat u = svd.matrixU(), v = svd.matrixV();

    Mat I = MatrixXd::Identity(3, 3);
    if (g.determinant() < 0) I(2, 2) = -1;
    Mat rot = u * I * v.transpose();
    
    Mat d = x * rot - y;
    val_t rmsd = 0; for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) rmsd += d(i, j) * d(i, j);
    rmsd = std::sqrt(rmsd / len);

    return Result{rot, c1, c2, rmsd};
}

template<typename T, typename U, typename F>
auto suppos(T &t, const U &m, const F &n) {
    auto sp = suppos(m, n);
    for (int i = 0; i < t.rows(); i++) for (int j = 0; j < 3; j++) t(i, j) -= sp.c1[j];
    t = t * sp.rot;
    for (int i = 0; i < t.rows(); i++) for (int j = 0; j < 3; j++) t(i, j) += sp.c2[j];
    return sp;
}

template<typename T, typename U, typename F, typename V>
auto suppos(T &src, const U &src_indices, const F &tgt, const V &tgt_indices) {
    MatrixXd m(src_indices.size(), 3), n(tgt_indices.size(), 3);
    for (int i = 0; i < src_indices.size(); i++) for (int j = 0; j < 3; j++) {
        m(i, j) = src(src_indices[i], j);
        n(i, j) = tgt(tgt_indices[i], j);
    }

    auto sp = suppos(m, n);

    for (int i = 0; i < src.rows(); i++) for (int j = 0; j < 3; j++) src(i, j) -= sp.c1[j];
    src = src * sp.rot;
    for (int i = 0; i < src.rows(); i++) for (int j = 0; j < 3; j++) src(i, j) += sp.c2[j];
    return sp;
}

template<typename T, typename U>
auto rmsd(const T &x, const U &y) {
    return suppos(x, y).rmsd;
}

template<typename Mat>
Mat suppos_axis_polar(val_t theta_o, val_t phi_o, val_t theta_n, val_t phi_n) {
    Mat m = Mat::Identity(3, 3);
    val_t ang=PI/2-phi_o;      if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
          ang=theta_o-theta_n; if (ang != 0) m *= x_rot_mat(std::cos(ang), std::sin(ang));
          ang=phi_n-PI/2;      if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
    return m;
}

template<typename Mat, typename O, typename N>
Mat suppos_axis_xyz(const O &o, const N &n) {
    Mat m = Mat::Identity(3, 3);
    val_t r, x, y, r1, x1, y1, r2, x2, y2, c, s, c1, c2, s1, s2;
    r=std::sqrt(o[0]*o[0]+o[1]*o[1]); x=o[0]; y=o[1];
    if (r != 0) {c=y/r; s=x/r; if (s != 0) m *= z_rot_mat(c, s);}
    r1=std::sqrt(o[0]*o[0]+o[1]*o[1]+o[2]*o[2]); x1=std::sqrt(o[0]*o[0]+o[1]*o[1]); y1=o[2];
    r2=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); x2=std::sqrt(n[0]*n[0]+n[1]*n[1]); y2=n[2];
    if(r1 != 0 && r2 != 0) {
        c1=y1/r1; s1=x1/r1; c2=y2/r2; s2=x2/r2; c=c1*c2+s1*s2; s=s1*c2-s2*c1;
        if (s != 0) m *= x_rot_mat(c, s);
    }
    r=std::sqrt(n[0]*n[0]+n[1]*n[1]); x=n[0]; y=n[1];
    if (r != 0) {c=y/r; s=-x/r; if (s != 0) m *= z_rot_mat(c, s);}
    return m;
}

} // namespace geom
} // namespace jian

