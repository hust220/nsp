#pragma once

#include "../utils/Debug.hpp"
#include "Job.hpp"
#include "smooth.hpp"
#include "CG.hpp"
#include "MC.hpp"
#include "../geom.hpp"

namespace jian {

class DGImpl : public dg::MC {
public:
    DGImpl() = default;
    DGImpl(const DGImpl &dg) =  default;

    DGImpl(const DistBoundType &dist_bound) {
        bound = dist_bound;
        if (bound.rows() != bound.cols()) throw "DGImpl initialize error!";
        len = bound.rows();
        dg_smooth(bound, 3);
    }

    DGImpl(const DistBoundType &dist_bound, const DihBound &dih_bound) : DGImpl(dist_bound) {
        _dih_bound = dih_bound;
    }

    Mat operator ()() {
        Debug::print("### Start DG...\n");
        b2d(); metric(); d2c(); cg();
        if (_dist_en > 20 || _dih_en > 20) mc();
        Debug::print("DG Energy (", _dist_en, ' ', _dih_en, ").\n");
        return c;
    }

    Mat operator ()(const Mat &b) {
//        std::cout << bound.rows() << ' ' << bound.cols() << ' ' << b.rows() << ' ' << b.cols() << std::endl;
//        bound.resize(b.rows(), b.cols());
        bound = b;
        if (b.rows() != b.cols()) throw "operator()(const Mat &) error!";
        len = bound.rows();
        dg_smooth(bound, _min_dist);
        return (*this)();
    }

    Mat operator ()(const Mat &b, const DihBound &d) {
//        bound.resize(b.rows(), b.cols());
        bound = b;
        _dih_bound = d;
        if (b.rows() != b.cols()) throw "operator()(const Mat &) error!";
        len = bound.rows();
        dg_smooth(bound, _min_dist);
        return (*this)();
    }

    void b2d() {
        /* randomly assign values to the boundance matrix */
        d.resize(len, len);
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (j == i) d(i, j) = 0;
            else d(i, j) = d(j, i) = rand() * (bound(i, j) - bound(j, i)) + bound(j, i);
        }
    }

    void metric() {
        /* set metric matrix */
        double temp = 0;
        for (int j = 0; j < len; j++) for (int k = j + 1; k < len; k++) temp += d(j, k) * d(j, k);
        temp /= double(len * len);
        double a[len];
        for (int i = 0; i < len; i++) {
            a[i] = 0;
            for (int j = 0; j < len; j++) a[i] += d(i, j) * d(i, j);
            a[i] /= double(len);
            a[i] -= temp;
        }
        m.resize(len, len);
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            m(i, j) = m(j, i) = (a[i] + a[j] - d(i, j) * d(i, j)) / 2.;
        }
    }

    void d2c() {
        /* transform bounds matrix to coordinates matrix */
        Mat A(len, len);
        for (int i = 0; i < len; i++) for (int j = 0; j < len; j++) A(i, j) = m(i, j);
        Eigen::SelfAdjointEigenSolver<Mat> eigensolver(A);
        double max1 = 0, max2 = 0, max3 = 0;
        int m1 = 0, m2 = 0, m3 = 0;
        for (int i = 0; i < len; i++) {
            double temp = abs(eigensolver.eigenvalues()[i]);
            if (temp > max1) {
                max3 = max2; m3 = m2; max2 = max1; m2 = m1; max1 = temp; m1 = i;
            } else if (temp > max2) {
                max3 = max2; m3 = m2; max2 = temp; m2 = i;
            } else if (temp > max3) {
                max3 = temp; m3 = i;
            }
        }
        c.resize(len, 3);
        for (int i = 0; i < len; i++) {
            double temp = eigensolver.eigenvalues()[m1];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            c(i, 0) = temp * eigensolver.eigenvectors()(i, m1);
            temp = eigensolver.eigenvalues()[m2];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            c(i, 1) = temp * eigensolver.eigenvectors()(i, m2);
            temp = eigensolver.eigenvalues()[m3];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            c(i, 2) = temp * eigensolver.eigenvectors()(i, m3);
        }
    }

    Mat c2d(const Mat &coord) {
        Mat dist(coord.rows(), coord.cols());
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (i == j) dist(i, j) = 0;
            else dist(i, j) = dist(j, i) = geom::distance(coord.row(i), coord.row(j));
        }
        return dist;
    }

};

} // namespace jian

