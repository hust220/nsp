#ifndef JIAN_GEOM_SUPPOS
#define JIAN_GEOM_SUPPOS

#include "../util/std.h"
#include "../util/mat.h"

namespace jian {
namespace geom {

template<typename T, typename U>
auto suppos(const T &m, const U &n) {
    using val_t = double;
    using Vec = VectorXd;
    using Mat = MatrixXd;
    using Result = struct {Mat rot; Vec c1; Vec c2; val_t rmsd;};

    if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) throw "jian::geom::suppos error!";
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

} // namespace geom
} // namespace jian

#endif

