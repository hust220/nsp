#pragma once

#include <Eigen/Dense>
#include "../utils/traits.hpp"

BEGIN_JN

namespace geom {

template<typename T>
inline auto square(const T &t) {
    return t * t;
}

template<class P1, class P2> 
inline auto distance(const P1 &p1, const P2 &p2) {
    return std::sqrt(square(p1[0]-p2[0])+square(p1[1]-p2[1])+square(p1[2]-p2[2]));
}

template<class P1, class P2> 
inline auto dist2(const P1 &p1, const P2 &p2) {
    return square(p1[0]-p2[0])+square(p1[1]-p2[1])+square(p1[2]-p2[2]);
}

template<typename P1, typename P2>
inline double norm(P1 &&p1, P2 &&p2, int n) {
    double sum = 0; for (int i = 0; i < n; i++) sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    return std::sqrt(sum);
}

template<class P1, class P2, class P3>
inline double angle(const P1 &p1, const P2 &p2, const P3 &p3) {
    double a1 = p1[0] - p2[0], a2 = p1[1] - p2[1], a3 = p1[2] - p2[2];
    double b1 = p3[0] - p2[0], b2 = p3[1] - p2[1], b3 = p3[2] - p2[2];
    double ra = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    double rb = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
    return acos((a1 * b1 + a2 * b2 + a3 * b3) / (ra * rb));
}

template<class T1, class T2, class T3, class T4>
inline double chirality(const T1 &p1, const T2 &p2, const T3 &p3, const T4 &p4) {
    double a[3] = {p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2]};
    double b[3] = {p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2]};
    double c[3] = {p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2]};
    double d[3] = {b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0]};
    return a[0]*d[0]+a[1]*d[1]+a[2]*d[2];
}

template<class P = Eigen::Vector3d, class P1, class P2, class P3>
P normal_vector(const P1 &p1, const P2 &p2, const P3 &p3) {
    double a1, a2, a3, b1, b2, b3, r;
    P p;
    a1 = p2[0] - p1[0]; a2 = p2[1] - p1[1]; a3 = p2[2] - p1[2];
    b1 = p3[0] - p2[0]; b2 = p3[1] - p2[1]; b3 = p3[2] - p2[2];
    p[0] = a2 * b3 - a3 * b2;
    p[1] = b1 * a3 - a1 * b3;
    p[2] = a1 * b2 - a2 * b1;
    r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    p[0] /= r; p[1] /= r; p[2] /= r;
    return p;
}

template<class P1, class P2, class P3, class P4>
inline double dihedral(const P1 &p1, const P2 &p2, const P3 &p3, const P4 &p4) {
    double a1, a2, a3, b1, b2, b3, c1, c2, c3;
    double a1_, a2_, a3_, b1_, b2_, b3_, c1_, c2_, c3_;
    double r, c, s;

    a1 = p1[0] - p2[0]; a2 = p1[1] - p2[1]; a3 = p1[2] - p2[2];
    b1 = p3[0] - p2[0]; b2 = p3[1] - p2[1]; b3 = p3[2] - p2[2];
    c1 = p4[0] - p2[0]; c2 = p4[1] - p2[1]; c3 = p4[2] - p2[2];
    if (b1 * b1 + b2 * b2 != 0) {
        r = sqrt(b1 * b1 + b2 * b2);
        c = b1 / r; s = -b2 / r;
        a1_ = c * a1 - s * a2; a2_ = s * a1 + c * a2; a3_ = a3;
        b1_ = c * b1 - s * b2; b2_ = s * b1 + c * b2; b3_ = b3;
        c1_ = c * c1 - s * c2; c2_ = s * c1 + c * c2; c3_ = c3;
    } else {
        a1_ = a1; a2_ = a2; a3_ = a3;
        b1_ = b1; b2_ = b2; b3_ = b3;
        c1_ = c1; c2_ = c2; c3_ = c3;
    }
    if (b1_ * b1_ + b3_ * b3_ != 0) {
        r = sqrt(b1_ * b1_ + b3_ * b3_);
        c = b1_ / r; s = b3_ / r;
        a1 = c * a1_ + s * a3_; a2 = a2_; a3 = -s * a1_ + c * a3_;
        b1 = c * b1_ + s * b3_; b2 = b2_; b3 = -s * b1_ + c * b3_;
        c1 = c * c1_ + s * c3_; c2 = c2_; c3 = -s * c1_ + c * c3_;
    } else {
        a1 = a1_; a2 = a2_; a3 = a3_;
        b1 = b1_; b2 = b2_; b3 = b3_;
        c1 = c1_; c2 = c2_; c3 = c3_;
    }
    if (a2 * a2 + a3 * a3 != 0) {
        r = sqrt(a2 * a2 + a3 * a3);
        c = a3 / r; s = a2 / r;
        a1_ = a1; a2_ = c * a2 - s * a3; a3_ = s * a2 + c * a3;
        b1_ = b1; b2_ = c * b2 - s * b3; b3_ = s * b2 + c * b3;
        c1_ = c1; c2_ = c * c2 - s * c3; c3_ = s * c2 + c * c3;
    }
//    if (c2_ * c2_ + c3_ * c3_ == 0) throw "jian::geom::dihedral error!";
    if (c2_ * c2_ + c3_ * c3_ == 0) return 0;
    double temp = acos(c3_ / sqrt(c2_ * c2_ + c3_ * c3_));
    if (c2_ > 0) {
        temp = -temp;
    }
    return temp;
} 

} // namespace geometry
END_JN

