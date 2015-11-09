#include "Point.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

namespace jian {

Point::Point(double xx, double yy, double zz, int i) {
    x = xx;
    y = yy;
    z = zz;
    type = i;
}

Point::Point(const Point &p) {
    x = p.x;
    y = p.y;
    z = p.z;
    type = p.type;
}

Point::Point(Point *p) {
    x = p->x;
    y = p->y;
    z = p->z;
    type = p->type;
}

double &Point::operator[](int n) {
    if (n == 0) {
        return x;
    } else if (n == 1) {
        return y;
    } else if (n == 2) {
        return z;
    } else {
        cerr << "Point::operator[] error!" << endl;
        exit(1);
    }
}

const double &Point::operator[](int n) const {
    if (n == 0) {
        return x;
    } else if (n == 1) {
        return y;
    } else if (n == 2) {
        return z;
    } else {
        cerr << "Point::operator[] error!" << endl;
        exit(1);
    }
}

Point *Point::normalVector(Point *p1, Point *p2, Point *p3) {
    double a1, a2, a3, b1, b2, b3, r;
    Point *p = new Point;
    a1 = p2->x - p1->x; a2 = p2->y - p1->y; a3 = p2->z - p1->z;
    b1 = p3->x - p2->x; b2 = p3->y - p2->y; b3 = p3->z - p2->z;
    p->x = a2 * b3 - a3 * b2;
    p->y = b1 * a3 - a1 * b3;
    p->z = a1 * b2 - a2 * b1;
    r = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);
    p->x /= r; p->y /= r; p->z /= r;
    return p;
}

Point *Point::normalVector(Point &p1, Point &p2, Point &p3) {
    double a1, a2, a3, b1, b2, b3, r;
    Point *p = new Point;
    a1 = p2.x - p1.x; a2 = p2.y - p1.y; a3 = p2.z - p1.z;
    b1 = p3.x - p2.x; b2 = p3.y - p2.y; b3 = p3.z - p2.z;
    p->x = a2 * b3 - a3 * b2;
    p->y = b1 * a3 - a1 * b3;
    p->z = a1 * b2 - a2 * b1;
    r = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);
    p->x /= r; p->y /= r; p->z /= r;
    return p;
}

double Point::dist(Obj<Point> p) {
    return dist(p.get());
}

double Point::dist(Point *p) {
    double temp;
    double i, j, k;

    temp = 0;
    i = x - p->x;
    j = y - p->y;
    k = z - p->z;
    temp = i * i + j * j + k * k;
    if (temp > 0) {
        temp = sqrt(temp);
    }
    return temp;
}

double Point::dist(const Point &p) {
    double temp;
    double i, j, k;

    temp = 0;
    i = x - p.x;
    j = y - p.y;
    k = z - p.z;
    temp = i * i + j * j + k * k;
    if (temp > 0) {
        temp = sqrt(temp);
    }
    return temp;
}

double Point::angle(Obj<Point> p1, Obj<Point> p2, Obj<Point> p3) {
    double a1, a2, a3, b1, b2, b3;
    double ra, rb;
    a1 = p1->x - p2->x; a2 = p1->y - p2->y; a3 = p1->z - p2->z;
    b1 = p3->x - p2->x; b2 = p3->y - p2->y; b3 = p3->z - p2->z;
    ra = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    rb = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
    return acos((a1 * b1 + a2 * b2 + a3 * b3) / (ra * rb)) * 180 / 3.1415926;
}

double Point::angle(Point *p1, Point *p2, Point *p3) {
    double a1, a2, a3, b1, b2, b3;
    double ra, rb;
    a1 = p1->x - p2->x; a2 = p1->y - p2->y; a3 = p1->z - p2->z;
    b1 = p3->x - p2->x; b2 = p3->y - p2->y; b3 = p3->z - p2->z;
    ra = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    rb = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
    return acos((a1 * b1 + a2 * b2 + a3 * b3) / (ra * rb)) * 180 / 3.1415926;
}

double Point::angle(Point &p1, Point &p2, Point &p3) {
    double a1, a2, a3, b1, b2, b3;
    double ra, rb;
    a1 = p1.x - p2.x; a2 = p1.y - p2.y; a3 = p1.z - p2.z;
    b1 = p3.x - p2.x; b2 = p3.y - p2.y; b3 = p3.z - p2.z;
    ra = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    rb = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
    return acos((a1 * b1 + a2 * b2 + a3 * b3) / (ra * rb)) * 180 / 3.1415926;
}

double Point::dihedral(Point &p1, Point &p2, Point &p3, Point &p4) {
    return dihedral(&p1, &p2, &p3, &p4);
}

double Point::dihedral(Obj<Point> p1, Obj<Point> p2, Obj<Point> p3, Obj<Point> p4) {
    return dihedral(p1.get(), p2.get(), p3.get(), p4.get());
}

double Point::dihedral(Point *p1, Point *p2, Point *p3, Point *p4) {
    double a1, a2, a3, b1, b2, b3, c1, c2, c3;
    double a1_, a2_, a3_, b1_, b2_, b3_, c1_, c2_, c3_;
    double r, c, s;

    if (p1 == NULL || p2 == NULL || p3 == NULL || p4 == NULL) {
        cerr << "Point::dihedral: error!" << endl;
        exit(1);
    }

    a1 = p1->x - p2->x; a2 = p1->y - p2->y; a3 = p1->z - p2->z;
    b1 = p3->x - p2->x; b2 = p3->y - p2->y; b3 = p3->z - p2->z;
    c1 = p4->x - p2->x; c2 = p4->y - p2->y; c3 = p4->z - p2->z;
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
    if (c2_ * c2_ + c3_ * c3_ == 0) {
        return 0;
//        cerr << "dihedral error" << endl;
//        exit(1);
    }
    double temp = acos(c3_ / sqrt(c2_ * c2_ + c3_ * c3_)) * 180 / 3.1415926;
    if (c2_ > 0) {
        temp = 360 - temp;
    }
    return temp;
}

Point *Point::grow(Point *pos1, Point *pos2, double angle, double dihedral, double bondLen) {
    double x = pos2->x - pos1->x;
    double y = pos2->y - pos1->y;
    double z = pos2->z - pos1->z;
    double l = sqrt(x * x + y * y + z * z);
    angle = (180 - angle) * 3.1415926 / 180;
    dihedral = dihedral * 3.1415926 / 180;
    double x0 = 0;
    double y0 = l + bondLen * cos(angle);
    double z0 = bondLen * sin(angle);

    double c = cos(dihedral);
    double s = sin(dihedral);
    double x1 = c * x0 + s * z0;
    double y1 = y0;
    double z1 = -s * x0 + c * z0;

    double r = sqrt(x * x + y * y + z * z);
    double rxy = sqrt(x * x + y * y);
    c = rxy / r;
    s = z / r;
    double x2 = x1;
    double y2 = c * y1 - s * z1;
    double z2 = s * y1 + c * z1;
    
    c = y / rxy;
    s = -x / rxy;
    double x3 = c * x2 - s * y2;
    double y3 = s * x2 + c * y2;
    double z3 = z2;
    
    x3 += pos1->x;
    y3 += pos1->y;
    z3 += pos1->z;
    Point *pos = new Point(x3, y3, z3);
    return pos;
}

Point *Point::rotate(Point *pos1, Point *pos2, Point *pos3, double dihedral, double bondLen1, double bondLen2) {
    double x2 = pos2->x - pos1->x;
    double y2 = pos2->y - pos1->y;
    double z2 = pos2->z - pos1->z;
    double x3 = pos3->x - pos1->x;
    double y3 = pos3->y - pos1->y;
    double z3 = pos3->z - pos1->z;
    
    double r3xy = sqrt(x3 * x3 + y3 * y3);
    double c3, s3;
    if (r3xy == 0) {
        c3 = 1;
        s3 = 0;
    } else {
        c3 = y3 / r3xy;
        s3 = x3 / r3xy;
    }
    double x2_1 = c3 * x2 - s3 * y2;
    double y2_1 = s3 * x2 + c3 * y2;
    double z2_1 = z2;
    double x3_1 = c3 * x3 - s3 * y3;
    double y3_1 = s3 * x3 + c3 * y3;
    double z3_1 = z3;

    double r3 = sqrt(x3 * x3 + y3 * y3 + z3 * z3);
    double c3_ = z3 / r3;
    double s3_ = r3xy / r3;
    double x2_2 = x2_1;
    double y2_2 = c3_ * y2_1 - s3_ * z2_1;
    double z2_2 = s3_ * y2_1 + c3_ * z2_1;
    double x3_2 = x3_1;
    double y3_2 = c3_ * y3_1 - s3_ * z3_1;
    double z3_2 = s3_ * y3_1 + c3_ * z3_1;

    dihedral = dihedral * 3.1415926 / 180;
    double c = cos(dihedral);
    double s = sin(dihedral);
    double x2_3 = c * x2_2 - s * y2_2;
    double y2_3 = s * x2_2 + c * y2_2;
    double z2_3 = z2_2;

    if (pos2->dist(pos1) > bondLen1 || pos2->dist(pos3) > bondLen2) {
        double r_ = sqrt(x2_3 * x2_3 + y2_3 * y2_3);
        double c_ = x2_3 / r_;
        double s_ = y2_3 / r_;
        x2_3 = z3_2 / 2 * cos(50 * 3.1415926 / 180) / sin(50 * 3.1415926 / 180) * c_;
        y2_3 = z3_2 / 2 * cos(50 * 3.1415926 / 180) / sin(50 * 3.1415926 / 180) * s_;
        z2_3 = z3_2 / 2;
    }

    double x2_4 = x2_3;
    double y2_4 = c3_ * y2_3 + s3_ * z2_3;
    double z2_4 = -s3_ * y2_3 + c3_ * z2_3;

    double x2_5 = c3 * x2_4 + s3 * y2_4;
    double y2_5 = -s3 * x2_4 + c3 * y2_4;
    double z2_5 = z2_4;

    Point *p = new Point(x2_5, y2_5, z2_5);
    return p;
}

void Point::minmax(double a, double b, double c, double d, double e, double &min, double &max) {
    double alpha = acos((b * b + d * d - a * a) / (2 * b * d));
    double beta = acos((b * b + c * c - e * e) / (2 * b * c));
    min = sqrt(c * c + d * d - cos(beta - alpha) * 2 * c * d);
    max = sqrt(c * c + d * d - cos(beta + alpha) * 2 * c * d);
}

double Point::gmMin(double a, double b, double c, double d, double e) {
    double alpha = acos((b * b + d * d - a * a) / (2 * b * d));
    double beta = acos((b * b + c * c - e * e) / (2 * b * c));
    return sqrt(c * c + d * d - cos(beta - alpha) * 2 * c * d);
}

double Point::gmMin(double a, double b, double c, double d1, double d2, double e1, double e2) {
    double temp1 = sqrt(b * b + c * c - 2 * b * c * (b * b + d1 * d1 - a * a) / (2 * b * d1));
    double temp2 = sqrt(b * b + a * a - 2 * b * a * (b * b + e1 * e1 - c * c) / (2 * b * e1));
    cout << temp1 << ' ' << temp2 << ' ' << e1 << ' ' << d1 << endl;
    if (e1 <= temp1 && d1 <= temp2) {
        return 0;
    } else if (e1 <= temp1 && d1 > temp2) {
        return fabs(d1 - c);
    } else if (e1 > temp1 && d1 <= temp2) {
        return fabs(e1 - a);
    } else {
        return Point::gmMin(a, b, c, d1, e1);
    }
    /*
    if (a + c <= b) {
        return Point::gmMin(a, b, c, d1, e1);
    } else {
        if (d1 > c && e1 > a) {
            return Point::gmMin(a, b, c, d1, e1);
        } else if (d1 > c && e1 <= a) {
            double temp = sqrt(b * b + c * c - 2 * b * c * (b * b + d1 * d1 - a * a) / (2 * b * d1));
            if (e1 <= temp) {
                return fabs(d1 - c);
            } else {
                return Point::gmMin(a, b, c, d1, e1);
            }
        } else if (d1 <= c && e1 > a) {
            double temp = sqrt(b * b + a * a - 2 * b * a * (b * b + e1 * e1 - c * c) / (2 * b * e1));
            if (d1 <= temp) {
                return fabs(e1 - a);
            } else {
                return Point::gmMin(a, b, c, d1, e1);
            }
        } else {
            return 0;
        }
    }
    */
}

double Point::gmMax(double a, double b, double c, double d, double e) {
    double alpha = acos((b * b + d * d - a * a) / (2 * b * d));
    double beta = acos((b * b + c * c - e * e) / (2 * b * c));
    return sqrt(c * c + d * d - cos(beta + alpha) * 2 * c * d);
}

Point *Point::rotate(Point *l1, Point *l2, double t) {
    Point *p = new Point;
    p->x = x - l1->x; p->y = y - l1->y; p->z = z - l1->z;
    Point *l = new Point(l2->x - l1->x, l2->y - l1->y, l2->z - l1->z);

    double r, c, s, x_, y_, z_;
    // rotate with z to put l on y-z plane
    r = sqrt(l->x * l->x + l->y * l->y);
    if (r != 0) {
        c = l->y / r;
        s = l->x / r;
        x_ = c * p->x - s * p->y;
        y_ = s * p->x + c * p->y;
        z_ = p->z;
        p->x = x_; p->y = y_; p->z = z_;
    }
    // rotate with x to put l on z
    r = sqrt(l->x * l->x + l->y * l->y + l->z * l->z);
    if (r != 0) {
        c = l->z / r;
        s = sqrt(l->x * l->x + l->y * l->y) / r;
        x_ = p->x;
        y_ = c * p->y - s * p->z;
        z_ = s * p->y + c * p->z;
        p->x = x_; p->y = y_; p->z = z_;
    }
    // rotate with z
    c = cos(t / 180.0 * 3.1415927);
    s = sin(t / 180.0 * 3.1415927);
    x_ = c * p->x - s * p->y;
    y_ = s * p->x + c * p->y;
    z_ = p->z;
    p->x = x_; p->y = y_; p->z = z_;
    // rotate with x
    r = sqrt(l->x * l->x + l->y * l->y + l->z * l->z);
    if (r != 0) {
        c = l->z / r;
        s = -sqrt(l->x * l->x + l->y * l->y) / r;
        x_ = p->x;
        y_ = c * p->y - s * p->z;
        z_ = s * p->y + c * p->z;
        p->x = x_; p->y = y_; p->z = z_;
    }
    // rotate with z
    r = sqrt(l->x * l->x + l->y * l->y);
    if (r != 0) {
        c = l->y / r;
        s = -l->x / r;
        x_ = c * p->x - s * p->y;
        y_ = s * p->x + c * p->y;
        z_ = p->z;
        p->x = x_; p->y = y_; p->z = z_;
    }
    p->x += l1->x; p->y += l1->y; p->z += l1->z;
    delete l;
    return p;
}

Point *Point::rotate(Point &l1, Point &l2, double t) {
    return rotate(&l1, &l2, t);
}

Point *Point::rotate(Point *l, double t) {
    Point *p = new Point;
    return rotate(p, l, t);
}

Point *Point::t2c(double *l, double *a, double *t, int n) {
    if (n < 4) {
        cerr << "Too few points!" << endl;
        exit(1);
    }
    Point *p = new Point[n];
    p[0].x = 0; p[0].y = 0; p[0].z = 0;
    p[1].x = 0; p[1].y = 0; p[1].z = l[0];
    p[2].x = l[1] * cos((a[0] - 90) / 180.0 * 3.1415927); p[2].y = 0; p[2].z = p[1].z + l[1] * sin((a[0] - 90) / 180.0 * 3.1415927);
    for (int i = 0; i < n - 3; i++) {
        Point *n = normalVector(p[i], p[i + 1], p[i + 2]);
        Point m;
        m.x = p[i + 2].x + n->x;
        m.y = p[i + 2].y + n->y;
        m.z = p[i + 2].z + n->z;
        Point *q = new Point;
        double a1 = p[i + 1].x - p[i + 2].x;
        double b1 = p[i + 1].y - p[i + 2].y;
        double c1 = p[i + 1].z - p[i + 2].z;
        double ratio = l[i + 2] / l[i + 1];
        q->x = ratio * a1 + p[i + 2].x;
        q->y = ratio * b1 + p[i + 2].y;
        q->z = ratio * c1 + p[i + 2].z;
        Point *r = q->rotate(p[i + 2], m, a[i + 1]);
        q = r->rotate(p[i + 1], p[i + 2], t[i]);
        p[i + 3].x = q->x;
        p[i + 3].y = q->y;
        p[i + 3].z = q->z;
        delete n;
        delete q;
        delete r;
    }
    return p;
}

Point *Point::t2c(double *l, double *a, double *t, Point *p_o, Point *p_n, int n) {
    if (n < 4) {
        cerr << "Too few points!" << endl;
        exit(1);
    }
    // p[2].x = l[1] * cos((a[0] - 90) / 180.0 * 3.1415927); p[2].y = 0; p[2].z = p[1].z + l[1] * sin((a[0] - 90) / 180.0 * 3.1415927);
    Point *p_temp = new Point[n + 3];
    p_temp[0].x = p_o[0].x; p_temp[0].y = p_o[0].y; p_temp[0].z = p_o[0].z;
    p_temp[1].x = p_o[1].x; p_temp[1].y = p_o[1].y; p_temp[1].z = p_o[1].z;
    p_temp[2].x = p_o[2].x; p_temp[2].y = p_o[2].y; p_temp[2].z = p_o[2].z;
    double *l_temp = new double[n + 1];
    for (int i = 0; i < n; i++) {
        l_temp[i] = l[i];
    }
    l_temp[n] = l[0];
    for (int i = 0; i < n; i++) {
        Point *n = normalVector(p_temp[i], p_temp[i + 1], p_temp[i + 2]);
        Point m;
        m.x = p_temp[i + 2].x + n->x;
        m.y = p_temp[i + 2].y + n->y;
        m.z = p_temp[i + 2].z + n->z;
        Point *q = new Point;
        double a1 = p_temp[i + 1].x - p_temp[i + 2].x;
        double b1 = p_temp[i + 1].y - p_temp[i + 2].y;
        double c1 = p_temp[i + 1].z - p_temp[i + 2].z;
        double ratio = l_temp[i + 1] / l_temp[i];
        q->x = ratio * a1 + p_temp[i + 2].x;
        q->y = ratio * b1 + p_temp[i + 2].y;
        q->z = ratio * c1 + p_temp[i + 2].z;
        Point *r = q->rotate(p_temp[i + 2], m, a[i]);
        q = r->rotate(p_temp[i + 1], p_temp[i + 2], t[i]);
        p_temp[i + 3].x = q->x;
        p_temp[i + 3].y = q->y;
        p_temp[i + 3].z = q->z;
        delete n;
        delete q;
        delete r;
    }
    Point *p = new Point[n];
    p[0].x = p_o[1].x; p[0].y = p_o[1].y; p[0].z = p_o[1].z;
    p[1].x = p_o[2].x; p[1].y = p_o[2].y; p[1].z = p_o[2].z;
    for (int i = 0; i < n - 2; i++) {
        p[i + 2].x = p_temp[i + 3].x;
        p[i + 2].y = p_temp[i + 3].y;
        p[i + 2].z = p_temp[i + 3].z;
    }
    p_n[0].x = p_temp[n].x; p_n[0].y = p_temp[n].y; p_n[0].z = p_temp[n].z;
    p_n[1].x = p_temp[n + 1].x; p_n[1].y = p_temp[n + 1].y; p_n[1].z = p_temp[n + 1].z;
    p_n[2].x = p_temp[n + 2].x; p_n[2].y = p_temp[n + 2].y; p_n[2].z = p_temp[n + 2].z;
    return p;
}

void Point::rotate(Point *p, int len, Point beg, Point end, double angle) {
    for (int i = 0; i < len; i++) {
        p[i].x -= beg.x;
        p[i].y -= beg.y;
        p[i].z -= beg.z;
    }
    Point temp(end.x - beg.x, end.y - beg.y, end.z - beg.z);

    double r, c, s, x_, y_, z_;
    r = sqrt(temp.x * temp.x + temp.y * temp.y);
    if (r != 0) {
        c = temp.y / r;
        s = temp.x / r;
        for (int i = 0; i < len; i++) {
            x_ = c * p[i].x - s * p[i].y;
            y_ = s * p[i].x + c * p[i].y;
            z_ = p[i].z;
            p[i].x = x_;
            p[i].y = y_;
            p[i].z = z_;
        }
    }
    r = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
    if (r != 0) {
        c = temp.z / r;
        s = sqrt(temp.x * temp.x + temp.y * temp.y) / r;
        for (int i = 0; i < len; i++) {
            x_ = p[i].x;
            y_ = c * p[i].y - s * p[i].z;
            z_ = s * p[i].y + c * p[i].z;
            p[i].x = x_;
            p[i].y = y_;
            p[i].z = z_;
        }
    }
    c = cos(angle / 180. * 3.1415927);
    s = sin(angle / 180. * 3.1415927);
    for (int i = 0; i < len; i++) {
        x_ = c * p[i].x - s * p[i].y;
        y_ = s * p[i].x + c * p[i].y;
        z_ = p[i].z;
        p[i].x = x_;
        p[i].y = y_;
        p[i].z = z_;
    }
    r = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
    if (r != 0) {
        c = temp.z / r;
        s = -sqrt(temp.x * temp.x + temp.y * temp.y) / r;
        for (int i = 0; i < len; i++) {
            x_ = p[i].x;
            y_ = c * p[i].y - s * p[i].z;
            z_ = s * p[i].y + c * p[i].z;
            p[i].x = x_;
            p[i].y = y_;
            p[i].z = z_;
        }
    }
    r = sqrt(temp.x * temp.x + temp.y * temp.y);
    if (r != 0) {
        c = temp.y / r;
        s = -temp.x / r;
        for (int i = 0; i < len; i++) {
            x_ = c * p[i].x - s * p[i].y;
            y_ = s * p[i].x + c * p[i].y;
            z_ = p[i].z;
            p[i].x = x_;
            p[i].y = y_;
            p[i].z = z_;
        }
    }

    for (int i = 0; i < len; i++) {
        p[i].x += beg.x;
        p[i].y += beg.y;
        p[i].z += beg.z;
    }
}

void Point::coincide(Point *p, int len, Point o, Point n) {
    double r, r2, c, s, x_, y_, z_;

    r2 = o.x * o.x + o.y * o.y;
    if (r != 0) {
        r = sqrt(r2);
        c = o.y / r;
        s = o.x / r;
        for (int i = 0; i < len; i++) {
            x_ = c * p[i].x - s * p[i].y;
            y_ = s * p[i].x + c * p[i].y;
            z_ = p[i].z;
            p[i].x = x_;
            p[i].y = y_;
            p[i].z = z_;
        }
    }

    r2 = o.x * o.x + o.y * o.y + o.z * o.z;
    if (r2 == 0) {
        cerr << "Point::coincide error1!" << endl;
        exit(1);
    }
    r = sqrt(r2);
    double a1 = acos(o.z / r);
    r2 = n.x * n.x + n.y * n.y + n.z * n.z;
    if (r2 == 0) {
        cerr << "Point::coincide error2!" << endl;
        exit(1);
    }
    r = sqrt(r2);
    double a2 = acos(n.z / r);
    double delta_a = a1 - a2;

    c = cos(delta_a);
    s = sin(delta_a);
    for (int i = 0; i < len; i++) {
        x_ = p[i].x;
        y_ = c * p[i].y - s * p[i].z;
        z_ = s * p[i].y + c * p[i].z;
        p[i].x = x_;
        p[i].y = y_;
        p[i].z = z_;
    }

    r2 = n.x * n.x + n.y * n.y;
    if (r2 != 0) {
        r = sqrt(r2);
        c = n.y / r;
        s = -n.x / r;
        for (int i = 0; i < len; i++) {
            x_ = c * p[i].x - s * p[i].y;
            y_ = s * p[i].x + c * p[i].y;
            z_ = p[i].z;
            p[i].x = x_;
            p[i].y = y_;
            p[i].z = z_;
        }
    }
}

double Point::chirality(Point *p1, Point *p2, Point *p3, Point *p4) {
    Point a(p1->x - p4->x, p1->y - p4->y, p1->z - p4->z);
    Point b(p2->x - p4->x, p2->y - p4->y, p2->z - p4->z);
    Point c(p3->x - p4->x, p3->y - p4->y, p3->z - p4->z);
    Point d(b.y * c.z - b.z * c.y, b.z * c.x - b.x * c.z, b.x * c.y - b.y * c.x);
    return a.x * d.x + a.y * d.y + a.z * d.z;
}

double Point::chirality(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {
    Point a(p1.x - p4.x, p1.y - p4.y, p1.z - p4.z);
    Point b(p2.x - p4.x, p2.y - p4.y, p2.z - p4.z);
    Point c(p3.x - p4.x, p3.y - p4.y, p3.z - p4.z);
    Point d(b.y * c.z - b.z * c.y, b.z * c.x - b.x * c.z, b.x * c.y - b.y * c.x);
    return a.x * d.x + a.y * d.y + a.z * d.z;
}

Point operator -(const Point &pt) {
    return Point(-pt.x, -pt.y, -pt.z);
}

} /// namespace jian


