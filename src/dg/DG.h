#ifndef DG_H
#define DG_H

#include <util/util.h>

namespace jian {

class DG {
public:
    /* random device */
    std::random_device rd;
    std::uniform_real_distribution<> unif_distr;

    int len = 0;

    MatrixXf bound;
    double _min_dist = 0;
    MatrixXf chir;
    MatrixXf d;
    MatrixXf m;
    MatrixXf c;
    MatrixXf g;
    double g2 = 0;
    double E = 0;
    double CH = 0;
    int view = 0;

    DG() : unif_distr(0, 1) {}

    DG(const MatrixXf &b, int min_dist = 0) : bound(b), _min_dist(min_dist), unif_distr(0, 1) {
        if (bound.rows() != bound.cols()) {
            std::cerr << "DG::DG error! The number of rows (" << bound.rows() << ") and that of columns (" << bound.cols() << ")in bound matrix must be equal. " << std::endl;
            exit(1);
        }
        len = bound.rows();
        smooth();
    }

    DG(const MatrixXf &b, const MatrixXf &c, int min_dist = 0) : 
      bound(b), chir(c) , _min_dist(min_dist), unif_distr(0, 1) {
        if (b.rows() != b.cols()) {
            std::cerr << "DG::DG error! The number of rows and that of columns in bound matrix must be equal. " << std::endl;
            exit(1);
        }
        len = bound.rows();
        smooth();
    }

    DG(const DG &dg): len(dg.len), bound(dg.bound), chir(dg.chir),
                      d(dg.d), m(dg.m), c(dg.c), g(dg.g), g2(dg.g2),
                      E(dg.E), CH(dg.CH), view(dg.view) {}

    DG &operator =(const DG &dg) {
        len = dg.len;
        bound = dg.bound;
        chir = dg.chir;
        d = dg.d;
        m = dg.m;
        c = dg.c;
        g = dg.g;
        g2 = dg.g2;
        E = dg.E;
        CH = dg.CH;
        view = dg.view;
        return (*this);
    }

    MatrixXf operator ()() {
        if (view) cerr << "bound: \n" << bound << endl;
        double E_min;
        MatrixXf c_min;
        for (int i = 0; i < 3; i++) {
            b2d();
            metric();
            d2c();
            if (view) cerr << c << endl;
            cg();
            if (view) cerr << c << endl;
            if (E > 20) {
                mc();
            }
            if (E <= 200) {
                return c;
            }
            if (i == 0) {
                E_min = E;
                c_min = c;
            } else {
                if (E_min > E) {
                    E_min = E;
                    c_min = c;
                }
            }
        }
        if (E != E_min) {
            E = E_min;
            c = c_min;
        }
        return c;
    }

    MatrixXf operator ()(const MatrixXf &b) {
        bound = b;
        if (b.rows() != b.cols()) throw "DG::operator()(const MatrixXf &) error!";
        len = bound.rows();
        smooth();
        return (*this)();
    }

    void smooth();
    void b2d();
    void metric();
    void d2c();
    void c2d();
    void gradient();
    double total_energy(const MatrixXf &coord);
    double atom_energy(const MatrixXf &coord, int i);
    void cg();
    void mc();
    double at(const MatrixXf &, int, int, int);
    void assign(MatrixXf &, int, int, double, int);

};

} /// namespace jian

#endif

