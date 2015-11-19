#include "SupPos.h"

namespace jian {

tuple<Matrix3f, Point, Point> SupPos::operator ()(const MatrixXf &m, const MatrixXf &n) {
    if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) {
        cerr << "SupPos::SupPos error!" << endl;
        exit(1);
    }
    len = m.rows();
    x.resize(len, 3);
    y.resize(len, 3);
    for (int i = 0; i < 3; i++) {
        c1[i] = 0;
        c2[i] = 0;
    }
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < 3; j++) {
            c1[j] += m(i, j);
            c2[j] += n(i, j);
        }
    }
    for (int i = 0; i < 3; i++) {
        c1[i] = c1[i] / len;
        c2[i] = c2[i] / len;
    }
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < 3; j++) {
            x(i, j) = m(i, j) - c1[j];
            y(i, j) = n(i, j) - c2[j];
        }
    }
    Matrix3f g = x.transpose() * y;

    JacobiSVD<Matrix3f> svd(g, ComputeFullU|ComputeFullV);
    Matrix3f u = svd.matrixU();
    Matrix3f v = svd.matrixV();

    double det = g.determinant();
    Matrix3f I;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                I(i, j) = 1;
            } else {
                I(i, j) = 0;
            }
        }
    }
    if (det < 0) {
        I(2, 2) = -1;
    }
    rot = u * I * v.transpose();
    
    MatrixXf x_ = x * rot;
    x_ = x_ - y;
    rmsd = 0;
    for (int i = 0; i < len; i++) {
        rmsd += x_(i, 0) * x_(i, 0) + x_(i, 1) * x_(i, 1) + x_(i, 2) * x_(i, 2);
    }
    rmsd /= (double) len;
    rmsd = sqrt(rmsd);

    return make_tuple(rot, c1, c2);
}

tuple<Matrix3f, Point, Point> SupPos::operator ()(MatrixXf &target,
                                                  const MatrixXf &m, 
                                                  const MatrixXf &n) {
    (*this)(m, n);
    for (int i = 0; i < target.rows(); i++) {
        for (int j = 0; j < target.cols(); j++) {
            target(i, j) -= c1[j];
        }
    }
    target = target * rot;
    for (int i = 0; i < target.rows(); i++) {
        for (int j = 0; j < target.cols(); j++) {
            target(i, j) += c2[j];
        }
    }
    return make_tuple(rot, c1, c2);
}

Matrix3f SupPos::operator ()(Matr_ *m, Matr_ *n) {
    if (m->row != n->row || m->col != 3 || n->col != 3) {
        cerr << "SupPos::SupPos error!" << endl;
        exit(1);
    }
    len = m->row;
    MatrixXf x(len, 3);
    MatrixXf y(len, 3);
    for (int i = 0; i < 3; i++) {
        c1[i] = 0;
        c2[i] = 0;
    }
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < 3; j++) {
            c1[j] += m->data[i][j];
            c2[j] += n->data[i][j];
        }
    }
    for (int i = 0; i < 3; i++) {
        c1[i] = c1[i] / len;
        c2[i] = c2[i] / len;
    }
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < 3; j++) {
            x(i, j) = m->data[i][j] - c1[j];
            y(i, j) = n->data[i][j] - c2[j];
        }
    }
    Matrix3f g = x.transpose() * y;

    JacobiSVD<Matrix3f> svd(g, ComputeFullU|ComputeFullV);
    Matrix3f u = svd.matrixU();
    Matrix3f v = svd.matrixV();

    double det = g.determinant();
    Matrix3f I;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                I(i, j) = 1;
            } else {
                I(i, j) = 0;
            }
        }
    }
    if (det < 0) {
        I(2, 2) = -1;
    }
    rot = u * I * v.transpose();
    
    MatrixXf x_(len, 3);
    x_ = x * rot;
    x_ = x_ - y;
    rmsd = 0;
    for (int i = 0; i < len; i++) {
        rmsd += x_(i, 0) * x_(i, 0) + x_(i, 1) * x_(i, 1) + x_(i, 2) * x_(i, 2);
    }
    rmsd /= (double) len;
    rmsd = sqrt(rmsd);

    return rot;
}

} /// namespace jian
