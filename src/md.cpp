#include <jian/utils/log.hpp>
#include <jian/matrix.hpp>
#include <jian/utils/file.hpp>
#include <jian/utils/rand.hpp>
#include <boost/format.hpp>
#include <string>
#include <deque>
#include "nsp.hpp"

namespace jian {

class MD {
public:
    int _n;
    Mat _c, _v, _a;
    double _dt = 0.01;
    double _box_size = 15;
    double _rmin = 3;

    void read_coords(std::string f) {
        FILE_OPEN(ifile, f.c_str());
        ifile >> _n;
        _c.resize(_n, 3);
        for (int i = 0; i < _n; i++) {
            for (int j = 0; j < 3; j++) {
                ifile >> _c(i, j);
            }
        }
        ifile.close();
    }

    void init_coords() {
        _c.resize(_n, 3);
        for (int i = 0; i < _n; i++) {
            for (int j = 0; j < 3; j++) {
                _c(i, j) = 20*(jian::rand()-0.5);
            }
        }
    }

    void init_vels() {
        _v.resize(_n, 3);
        for (int i = 0; i < _n; i++) {
            for (int j = 0; j < 3; j++) {
                _v(i, j) = 10*(jian::rand()-0.5);
            }
        }
    }

    void dist(Vec &r, double &d, int m, int n) {
        r[0] = _c(n, 0) - _c(m, 0);
        r[1] = _c(n, 1) - _c(m, 1);
        r[2] = _c(n, 2) - _c(m, 2);
        d = r.norm();
        r /= d;
    }

    void set_a() {
        Vec r(3);
        double d, a, b;
        _a.setZero(_n, 3);
        for (int i = 0; i < _n; i++) {
            for (int j = i+1; j < _n; j++) {
                dist(r, d, i, j);
                if (d < 8) {
                    a = std::pow(_rmin,6);
                    b = std::pow(d,6);
                    a = 12*a*a/b/b/d-6*a/b/d;
                    for (int k = 0; k < 3; k++) {
                        b = a*r[k];
                        _a(i, k) += -b;
                        _a(j, k) += b;
                    }
                }
            }
        }
    }

    void cycle_c(double &c) {
        if (c>0) {
            c -= int((c+10)/20)*20;
        } else {
            c += int((-c+10)/20)*20;
        }
    }

    void move() {
        for (int i = 0; i < _n; i++) {
            for (int j = 0; j < 3; j++) {
                _c(i, j) += _v(i, j) * _dt + 0.5 * _a(i, j) * _dt * _dt;
                cycle_c(_c(i, j));
                _v(i, j) += _a(i, j) * _dt;
            }
        }
    }

    void run(int steps) {
        _a.resize(_n, 3);
        print_coords(0);
        for (int i = 0; i < steps; i++) {
            set_a();
//            LOGI << _a << std::endl;
            move();
            print_coords(i+1);
        }
    }

    void print_coords(int n) {
        LOGI << "MODEL " << n+1 << std::endl;
        for (int i = 0; i < _n; i++) {
            LOGI << (boost::format("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12s  \n") % 
                                    (i+1) % "X" % "X" % "X" % (i+1) % 
                                    _c(i,0) % _c(i,1) % _c(i,2) % 1.00 % 0.00 % "X");
        }
        LOGI << "ENDMDL" << std::endl;;
    }

};

REGISTER_NSP_COMPONENT(md) {
    MD md;
//    md.read_coords(par["c"][0]);
    md._n = std::stoi(par["n"][0]);;
    md.init_coords();
    md.init_vels();
    md.run(std::stoi(par["steps"][0]));
}

} // namespace jian

