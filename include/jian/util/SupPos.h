#ifndef JIAN_SUPPOS_H
#define JIAN_SUPPOS_H

#include "std.h"
#include "Point.h"
#include "mat.h"

namespace jian {

class SupPos {
public:
    MatrixXf x;
    MatrixXf y;
    int len = 0;
    Matrix3f rot;
    Point c1;
    Point c2;
    double rmsd = 0;

    std::tuple<Matrix3f, Point, Point> operator ()(const MatrixXf &m, const MatrixXf &n) {
        if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) throw "SupPos::SupPos error!";
        len = m.rows();
        x.resize(len, 3); y.resize(len, 3);
        for (int i = 0; i < 3; i++) c1[i] = c2[i] = 0;
        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
                c1[j] += m(i, j);
                c2[j] += n(i, j);
        }
        for (int i = 0; i < 3; i++) {
            c1[i] = c1[i] / len;
            c2[i] = c2[i] / len;
        }
        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
                x(i, j) = m(i, j) - c1[j];
                y(i, j) = n(i, j) - c2[j];
        }
        Matrix3f g = x.transpose() * y;

        JacobiSVD<Matrix3f> svd(g, ComputeFullU|ComputeFullV);
        Matrix3f u = svd.matrixU();
        Matrix3f v = svd.matrixV();

        double det = g.determinant();
        Matrix3f I;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
                if (i == j) I(i, j) = 1; else I(i, j) = 0;
        }
        if (det < 0) I(2, 2) = -1;
        rot = u * I * v.transpose();
        
        MatrixXf x_ = x * rot;
        x_ = x_ - y;
        rmsd = 0;
        for (int i = 0; i < len; i++) {
            rmsd += x_(i, 0) * x_(i, 0) + x_(i, 1) * x_(i, 1) + x_(i, 2) * x_(i, 2);
        }
        rmsd = sqrt(rmsd / len);

        return std::make_tuple(rot, c1, c2);
    }

    std::tuple<Matrix3f, Point, Point> operator ()(MatrixXf &target, const MatrixXf &m, const MatrixXf &n) {
        (*this)(m, n);
        for (int i = 0; i < target.rows(); i++) for (int j = 0; j < target.cols(); j++) target(i, j) -= c1[j];
        target = target * rot;
        for (int i = 0; i < target.rows(); i++) for (int j = 0; j < target.cols(); j++) target(i, j) += c2[j];
        return std::make_tuple(rot, c1, c2);
    }

//    Matrix3f operator ()(Matr_ *m, Matr_ *n) {
//        if (m->row != n->row || m->col != 3 || n->col != 3) throw "SupPos::SupPos error";
//        len = m->row;
//        MatrixXf x(len, 3);
//        MatrixXf y(len, 3);
//        for (int i = 0; i < 3; i++) c1[i] = c2[i] = 0;
//        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
//                c1[j] += m->data[i][j];
//                c2[j] += n->data[i][j];
//        }
//        for (int i = 0; i < 3; i++) {
//            c1[i] = c1[i] / len;
//            c2[i] = c2[i] / len;
//        }
//        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
//                x(i, j) = m->data[i][j] - c1[j];
//                y(i, j) = n->data[i][j] - c2[j];
//        }
//        Matrix3f g = x.transpose() * y;
//
//        JacobiSVD<Matrix3f> svd(g, ComputeFullU|ComputeFullV);
//        Matrix3f u = svd.matrixU();
//        Matrix3f v = svd.matrixV();
//
//        double det = g.determinant();
//        Matrix3f I;
//        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
//                if (i == j) I(i, j) = 1; else I(i, j) = 0;
//        }
//        if (det < 0) I(2, 2) = -1;
//        rot = u * I * v.transpose();
//        
//        MatrixXf x_(len, 3);
//        x_ = x * rot; x_ = x_ - y; rmsd = 0;
//        for (int i = 0; i < len; i++) {
//            rmsd += x_(i, 0) * x_(i, 0) + x_(i, 1) * x_(i, 1) + x_(i, 2) * x_(i, 2);
//        }
//        rmsd = sqrt(rmsd / len);
//
//        return rot;
//    }

};

} /// namespace jian

#endif

