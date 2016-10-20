#include "suppos.hpp"

namespace jian {
namespace geom {

suppos_t suppos(const Mat &m, const Mat &n) {
    if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) {
        throw std::string("jian::geom::suppos error! (") + JN_STR(m.rows()) + ' ' + JN_STR(m.cols()) + ") (" + JN_STR(n.rows()) + ' ' + JN_STR(n.cols()) + ")\n";
    }
    int len = m.rows(); Mat x = m, y = n;
    Vec c1 = Vec::Zero(3), c2 = Vec::Zero(3);
    for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) { c1[j] += x(i, j); c2[j] += y(i, j); }
    for (int i = 0; i < 3; i++) { c1[i] = c1[i] / len; c2[i] = c2[i] / len; }
    for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) { x(i, j) -= c1[j]; y(i, j) -= c2[j]; }

    Mat g = x.transpose() * y;
    Eigen::JacobiSVD<Mat> svd(g, Eigen::ComputeFullU|Eigen::ComputeFullV);
    Mat u = svd.matrixU(), v = svd.matrixV();

    Mat I = Eigen::MatrixXd::Identity(3, 3);
    if (g.determinant() < 0) I(2, 2) = -1;
    Mat rot = u * I * v.transpose();
    
    Mat d = x * rot - y;
    val_t rmsd = 0; for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) rmsd += d(i, j) * d(i, j);
    rmsd = std::sqrt(rmsd / len);

    return suppos_t{rot, c1, c2, rmsd};
}

val_t rmsd(const Mat &m, const Mat &n) {
    int len = m.rows(); 
    Mat x = m, y = n;
    Vec c1 = Vec::Zero(3), c2 = Vec::Zero(3);
    for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) { c1[j] += x(i, j); c2[j] += y(i, j); }
    for (int i = 0; i < 3; i++) { c1[i] = c1[i] / len; c2[i] = c2[i] / len; }
    for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) { x(i, j) -= c1[j]; y(i, j) -= c2[j]; }

    Mat g = x.transpose() * y;
    Eigen::JacobiSVD<Mat> svd(g, Eigen::ComputeFullU|Eigen::ComputeFullV);

    Mat I = Eigen::MatrixXd::Identity(3, 3);
    if (g.determinant() < 0) I(2, 2) = -1;

    Mat d = x * (svd.matrixU() * I * svd.matrixV().transpose()) - y;
    val_t rmsd = 0; 
    for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) rmsd += d(i, j) * d(i, j);
    rmsd = std::sqrt(rmsd / len);

    return rmsd;
}

Mat suppos_axis_polar(double theta_o, double phi_o, double theta_n, double phi_n) {
    Mat m = Mat::Identity(3, 3);
    double ang=PI/2-phi_o;      if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
          ang=theta_o-theta_n; if (ang != 0) m *= x_rot_mat(std::cos(ang), std::sin(ang));
          ang=phi_n-PI/2;      if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
    return m;
}

} // namespace geom
} // namespace jian

