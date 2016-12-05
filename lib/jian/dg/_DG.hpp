#pragma once

#include "Job.hpp"
#include "smooth.hpp"
#include "CG.hpp"
#include "MC.hpp"
#include "../geom.hpp"

BEGIN_JN

class DG : public MC {
public:
    using Mat = MatrixXd;

    Mat _dg_b;
    Mat _dg_d;
    Mat _dg_m;
    Mat _dg_c;
    Mat _dg_g;
    DihBound _dg_dihs;
    double _dg_min_dist = 5;
    double _dg_g2 = 0;
    double _dg_en_dist = 0;
    double _dg_en_dih = 0;

    int _mc_index_atom;
    int _mc_index_coord;
    double _mc_shift;

    int dg_triangle_smoothing() {
        int f = 0;
        int len = _dg_b.rows(); 
        double d;
        FOR((i, len), FOR((j, i+1, len), FOR((k, j+1, len), 
           d = _dg_b(i, k) + _dg_b(j, k); if (_dg_b(i, j) > d) {_dg_b(i, j) = d; f++;} 
           if (_dg_b(j, i) > d) throw "Wrong bound matrix!";
           d = _dg_b(i, j) + _dg_b(j, k); if (_dg_b(i, k) > d) {_dg_b(i, k) = d; f++;} 
           if (_dg_b(k, i) > d) throw "Wrong bound matrix!";
           d = _dg_b(i, j) + _dg_b(i, k); if (_dg_b(j, k) > d) {_dg_b(j, k) = d; f++;} 
           if (_dg_b(k, j) > d) throw "Wrong bound matrix!";
           d = _dg_b(k, i) - _dg_b(j, k); if (_dg_b(j, i) < d) {_dg_b(j, i) = d; f++;}
           d = _dg_b(k, j) - _dg_b(i, k); if (_dg_b(j, i) < d) {_dg_b(j, i) = d; f++;}
           d = _dg_b(j, i) - _dg_b(j, k); if (_dg_b(k, i) < d) {_dg_b(k, i) = d; f++;}
           d = _dg_b(k, j) - _dg_b(i, j); if (_dg_b(k, i) < d) {_dg_b(k, i) = d; f++;}
           d = _dg_b(j, i) - _dg_b(i, k); if (_dg_b(k, j) < d) {_dg_b(k, j) = d; f++;}
           d = _dg_b(k, i) - _dg_b(i, j); if (_dg_b(k, j) < d) {_dg_b(k, j) = d; f++;}
        )));
        return f;
    }

    int dg_tetrangle_smoothing() {
        int f = 0;
        int len = _dg_b.rows(); 
        double d; 
    //    FOR((i, len), FOR((j, i+1, len), FOR((k, j+1, len), FOR((l, k+1, len),
    //        
    //    ))));
        return f;
    }

    void dg_smooth() {
        // check minimum distance
        int len = _dg_b.rows();
        FOR((i, _dg_b.rows()), FOR((j, _dg_b.cols()), 
            IF(i != j && _dg_b(i, j) < min_dist, _dg_b(i, j) = min_dist)));

        // Main Loop
        int max_step = 100;
        while (dg_triangle_smoothing() + dg_tetrangle_smoothing() != 0 && max_step-- > 0);
    }

    void dg_b2d() {
        _dg_d.resize(len, len);
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (j == i) _dg_d(i, j) = 0;
            else _dg_d(i, j) = _dg_d(j, i) = rand() * (_dg_b(i, j) - _dg_b(j, i)) + _dg_b(j, i);
        }
    }

    void dg_metric() {
        double temp = 0;
        for (int j = 0; j < len; j++) for (int k = j + 1; k < len; k++) temp += _dg_d(j, k) * _dg_d(j, k);
        temp /= double(len * len);
        double a[len];
        for (int i = 0; i < len; i++) {
            a[i] = 0;
            for (int j = 0; j < len; j++) a[i] += _dg_d(i, j) * _dg_d(i, j);
            a[i] /= double(len);
            a[i] -= temp;
        }
        _dg_m.resize(len, len);
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            _dg_m(i, j) = _dg_m(j, i) = (a[i] + a[j] - _dg_d(i, j) * _dg_d(i, j)) / 2.;
        }
    }

    void dg_d2c() {
        /* transform bounds matrix to coordinates matrix */
        Mat A(len, len);
        for (int i = 0; i < len; i++) for (int j = 0; j < len; j++) A(i, j) = _dg_m(i, j);
        SelfAdjointEigenSolver<Mat> eigensolver(A);
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
        _dg_c.resize(len, 3);
        for (int i = 0; i < len; i++) {
            double temp = eigensolver.eigenvalues()[m1];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            _dg_c(i, 0) = temp * eigensolver.eigenvectors()(i, m1);
            temp = eigensolver.eigenvalues()[m2];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            _dg_c(i, 1) = temp * eigensolver.eigenvectors()(i, m2);
            temp = eigensolver.eigenvalues()[m3];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            _dg_c(i, 2) = temp * eigensolver.eigenvectors()(i, m3);
        }
    }

    void dg_init() {
        Debug::println("### DG Initialize...");
        dg_smooth();
        Debug::println("### DG Initialization Done.");
    }

    void dg() {
        Debug::println("### DG Run...");
        dg_b2d(); 
        dg_metric(); 
        dg_d2c(); 
        mc();
        Debug::println("DG Done. Energy (", _dist_en, ' ', _dih_en, ").");
    }

    void mc_select() {
        _mc_index_atom = int(rand() * _dg_c.rows());
        _mc_index_coord = int(rand() * 3);
    }

    void mc_partial_energy() {}

    void mc_sample() {
        mc_backup();
        _dg_c(_mc_index_atom, _mc_index_coord) += (rand()-0.5) * 2;
    }

    void mc_backup() {
        _mc_old_coord = _dg_c(_mc_index_atom, _mc_index_coord);
    }

    void mc_back() {
        _dg_c(_mc_index_atom, _mc_index_coord) = _mc_old_coord;
    }

};

END_JN

