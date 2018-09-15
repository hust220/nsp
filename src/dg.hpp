#pragma once

#include "rand.hpp"
#include "matrix.hpp"
#include "dg_smooth.hpp"
#include "dg_move.hpp"
#include "geom.hpp"
#include "log.hpp"
#include "dg.hpp"
#include "jian.hpp"

namespace jian {

struct DgEnergy {
    double bound;
    double distance;
    double dihedral;

    double sum() const {
        return bound + distance + dihedral;
    }
};

struct DgConstraint {
    std::vector<int> key;
    std::vector<double> val;
};

class DgConstraints : public std::vector<DgConstraint> {
public:
    template<typename K1, typename K2>
    bool key_equal(const K1 &k1, const K2 &k2) {
        if (k1.size() != k2.size()) return false;
        return std::all_of(k1.begin(), k1.end(), [&k2](auto &&n1){
            return std::any_of(k2.begin(), k2.end(), [&n1](auto &&n2){
                return n1 == n2;
            });
        });
    }

    void append(const DgConstraint &constraint) {
        for (auto && c : *this) {
            if (key_equal(c.key, constraint.key)) {
                c.val = constraint.val;
                return;
            }
        }
        this->push_back(constraint);
    }
};

/**
 * Sampling by using Distance Geometry algorithm.
 *
 * Usage:
 *
 *   1) auto dg = std::make_unique<DG>(length, default_range, distance_constraints, dihedral_constraint);
 *      auto mat1 = dg->sample();
 *      auto mat2 = dg->sample();
 *
 *   2) auto dg = std::make_unique<DG>(distance_bound_matrix);
 *      auto mat1 = dg->sample();
 *      auto mat2 = dg->sample();
 *
 *   3) Set log file
 *      auto dg = std::make_unique<DG>(distance_bound_matrix);
 *      dg->set_log_file(std::string); // "std.out" or "std.err" are standard output and error pipe
 */
class Dg {
public:
    Log log;

    Mat bound; // distance bound matrix.
    double min_dist = 0;

    int len = 0; // must be initialized

    Mat d; // distance matrix
    Mat m; // metrix matrix
    Mat c; // coordinate matrix
    Mat g; // gradient matrix
    double g2 = 0;

    int num_mc_steps = 2000; // number mc steps
    double mc_tempr = 20; // mc temperature

    std::array<double, 2> default_range;

    std::vector<DgConstraint> distance_constraints;
    std::vector<DgConstraint> dihedral_constraints;
    std::vector<std::vector<int>> involved_distance_constraints;
    std::vector<std::vector<int>> involved_dihedral_constraints;

    Dg(const Mat &bound_, double min_dist_ = 0) : bound(bound_), min_dist(min_dist_) {
        // Set length
        len = bound.rows();

        // check minimum distance
        int len = bound.rows();
        for (int i = 0; i < bound.rows(); i++) for (int j = 0; j < bound.cols(); j++) {
            if (i != j && bound(i, j) < min_dist) {
                bound(i, j) = min_dist;
            }
        }
    }

    Dg(int len_,
       const std::array<double, 2> &default_range_,
       const std::vector<DgConstraint> &distance_constraints_,
       const std::vector<DgConstraint> &dihedral_constraints_ = {})
    {
        len = len_;
        default_range = default_range_;
        distance_constraints = distance_constraints_;
        dihedral_constraints = dihedral_constraints_;

        for (auto && constraint : distance_constraints) {
            if (constraint.val.size() == 1) constraint.val.push_back(constraint.val[0]);
        }

        sort_constraints(distance_constraints);
        sort_constraints(dihedral_constraints);

        set_involved_constraints(distance_constraints, involved_distance_constraints);
        set_involved_constraints(dihedral_constraints, involved_dihedral_constraints);

//        std::cout << bound << std::endl;

        log << "[*] DG: Set Bound." << std::endl;
        init_bound();

//        std::cout << bound << std::endl;

        // check minimum distance
        int len = bound.rows();
        for (int i = 0; i < bound.rows(); i++) for (int j = 0; j < bound.cols(); j++) {
            if (i != j && bound(i, j) < min_dist) {
                bound(i, j) = min_dist;
            }
        }
    }

    /**
     * Sample the distance once at one time.
     */
    void sample_distance(Mat &b, const std::array<int, 2> &pair) {
        int i = pair[0];
        int j = pair[1];
        double u = b(i, j);
        double l = b(j, i);
        b(i, j) = b(j, i) = rand() * (u - l) + l;
    }

    /**
     * Sample the matrix by a iterative method.
     */
    Mat sample_it() {
        // Make a copy of the bound
        Mat b = bound;

        log << "[*] Start DG Sampling by a iterative method..." << std::endl;

        std::array<int, 2> pair;

        while (!dg_smooth(b, pair)) {
//            std::cout << "Bound:" << std::endl;
//            std::cout << b << std::endl;
//
//            std::cout << "Maximum distance pair:" << std::endl;
//            std::cout << pair[0] << ' ' << pair[1] << std::endl;

            sample_distance(b, pair);

//            std::cout << "Sample distance:" << std::endl;
//            std::cout << b << std::endl;
        }

//        std::cout << "[*] DG: Randomly build a distance matrix." << std::endl;
        b2d(b);
//        std::cout << d << std::endl;

        log << "[*] DG: Construct metric matrix." << std::endl;
        metric();

//        std::cout << "[*] DG: Embed." << std::endl;
        embed();
//        std::cout << c << std::endl;

//        log << "[*] DG: Conjugate Gradient Optimization." << std::endl;
//        cg();
//
//        log << "[*] DG: Monte Carlo Optimization." << std::endl;
//        mc();

        return c;
    }

    Mat sample() {
        log << "[*] Start DG Sampling..." << std::endl;

        log << "[*] DG: Randomly build a distance matrix." << std::endl;
        b2d(bound);

        log << "[*] DG: Construct metric matrix." << std::endl;
        metric();

        log << "[*] DG: Embed." << std::endl;
        embed();

        log << "[*] DG: Conjugate Gradient Optimization." << std::endl;
        cg();

        log << "[*] DG: Monte Carlo Optimization." << std::endl;
        mc();

        return c;
    }

    void set_log_file(std::string filename)
    {
        log.file(filename);
    }

    /*
     * Transform coordiantes to distance matrix.
     */
    Mat c2d(const Mat &coord) {
        int len = coord.rows();

        Mat dist(len, len);
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (i == j) dist(i, j) = 0;
            else dist(i, j) = dist(j, i) = geom::distance(coord.row(i), coord.row(j));
        }
        return dist;
    }

private:
    void init_bound() {
        bound.resize(len, len);
        for (int i = 0; i < len; i++) {
            for (int j = i; j < len; j++) {
                if (i == j) {
                    bound(i, j) = 0;
                }
                else {
                    bound(i, j) = default_range[1];
                    bound(j, i) = default_range[0];
                }
            }
        }
        for (auto && c : distance_constraints) {
            bound(c.key[0], c.key[1]) = c.val[1];
            bound(c.key[1], c.key[0]) = c.val[0];
        }
    }

    template<typename T>
    void sort_constraints(T &&constraints) {
        for (auto && constraint : constraints) {
            std::sort(constraint.key.begin(), constraint.key.end(), [](int a, int b){
                return a < b;
            });
        }
        std::sort(constraints.begin(), constraints.end(), [](const DgConstraint &c1, const DgConstraint &c2) {
            return c1.key[0] < c2.key[0];
        });
    }

    template<typename T, typename U>
    void set_involved_constraints(T &&constraints, U &&involved_constraints) {
        involved_constraints.resize(len);
        for (int i = 0; i < constraints.size(); i++) {
            involved_constraints[constraints[i].key[0]].push_back(i);
            involved_constraints[constraints[i].key[1]].push_back(i);
        }
    }

    void cg() {
        auto en_t = total_energy(c);
        double oldE = en_t.sum();

        gradient();
        double oldG2 = g2;

        Mat d_o = Mat::Zero(len, 3);
        Mat d_n(len, 3);

        double a = 0.0001;
        double beta = 0;
        int upd = 0;
        int step = 0;

        while (upd < 500) {
            double tempG2 = g2, tempE = en_t.sum();
            step++;

            std::cout << format("Step: %05d, TotalEnergy: %8.3f, BoundEnergy: %8.3f, DistanceEnergy: %8.3f, DihedralEnergy: %8.3f",
                                step, en_t.sum(), en_t.bound, en_t.distance, en_t.dihedral) << std::endl;

            if (upd != 0) beta = g2 / oldG2;
            if (upd % 10 == 0) beta = 0;
            for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
                d_n(i, j) = -g(i, j) + beta * d_o(i, j);
                c(i, j) += a * d_n(i, j);
            }

            en_t = total_energy(c);
            gradient();

            if (en_t.sum() >= tempE || g2 >= tempG2) {
                for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) c(i, j) -= a * d_n(i, j);
                en_t = total_energy(c);
                gradient();
                a *= 0.5;
                if (a < 1.e-12) break;
            } else {
                upd++;
                a *= 2;
                for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) d_o(i, j) = d_n(i, j);
                oldG2 = tempG2;
                oldE = tempE;
                if (en_t.sum() < 1.e-12 || g2 < 1.e-12) break;
            }
        }
    }

    /**
     * Randomly assign values to the boundance matrix.
     */
    void b2d(Mat &b) {
        int len = b.rows();

        d.resize(len, len);
        for (int i = 0; i < len; i++) {
            for (int j = i; j < len; j++) {
                if (j == i) {
                    d(i, j) = 0;
                }
                else {
                    d(i, j) = d(j, i) = rand() * (b(i, j) - b(j, i)) + b(j, i);
                }
            }
        }
    }

    /*
     * Set the metric matrix.
     */
    void metric() {
        double temp = 0;
        for (int j = 0; j < len; j++) for (int k = j + 1; k < len; k++) temp += d(j, k) * d(j, k);
        temp /= double(len * len);
        std::vector<double> a(len);
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

    /*
     * transform metric matrix to coordinates matrix
     */
    void embed() {
        Mat A(len, len);
        for (int i = 0; i < len; i++) for (int j = 0; j < len; j++) A(i, j) = m(i, j);
//        Mat A = m;

        Eigen::SelfAdjointEigenSolver<Mat> eigensolver(A);

//        std::cout << "Eigen vectors:" << std::endl;
//        std::cout << eigensolver.eigenvectors() << std::endl;
//
//        std::cout << "Eigen values:" << std::endl;
//        std::cout << eigensolver.eigenvalues() << std::endl;

        auto & v = eigensolver.eigenvalues();

        std::vector<int> inds(len);
        std::iota(inds.begin(), inds.end(), 0);
//        std::sort(inds.begin(), inds.end(), [&eigensolver](int i, int j){ return std::fabs(eigensolver.eigenvalues[i]) > std::fabs(eigensolver.eigenvalues[j]); });
        std::sort(inds.begin(), inds.end(), [&v](int i, int j){ return std::fabs(v[i]) > std::fabs(v[j]); });

        std::cout << "EigenValues: ";
        for (auto && i : inds) {
            std::cout << std::fabs(v[i]) << ' ';
        }
        std::cout << std::endl;

//        double max1 = 0, max2 = 0, max3 = 0;
//        int m1 = 0, m2 = 0, m3 = 0;
//        for (int i = 0; i < len; i++) {
//            double temp = abs(eigensolver.eigenvalues()[i]);
//            if (temp > max1) {
//                max3 = max2; m3 = m2; max2 = max1; m2 = m1; max1 = temp; m1 = i;
//            } else if (temp > max2) {
//                max3 = max2; m3 = m2; max2 = temp; m2 = i;
//            } else if (temp > max3) {
//                max3 = temp; m3 = i;
//            }
//        }

//        std::cout << m1 << ' ' << m2 << ' ' << m3 << std::endl;

        c.resize(len, 3);
        for (int i = 0; i < len; i++) {
            double temp = eigensolver.eigenvalues()[inds[0]];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            c(i, 0) = temp * eigensolver.eigenvectors()(i, inds[0]);
            temp = eigensolver.eigenvalues()[inds[1]];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            c(i, 1) = temp * eigensolver.eigenvectors()(i, inds[1]);
            temp = eigensolver.eigenvalues()[inds[2]];
            if (temp > 0) temp = sqrt(temp); else if (temp < 0) temp = -sqrt(-temp);
            c(i, 2) = temp * eigensolver.eigenvectors()(i, inds[2]);
        }
    }

    /**
     * Calculate the gradient.
     */
    void gradient() {
        double err = 1.e-6;

        // Gradient
        g.resize(len, 3);
        g2 = 0;
        for (int i = 0; i < len; i++) {
            for (int t = 0; t < 3; t++) {
                c(i, t) -= err; double d1 = atom_energy(c, i); c(i, t) += err;
                c(i, t) += err; double d2 = atom_energy(c, i); c(i, t) -= err;
                g(i, t) = (d2 - d1) / (2 * err);
                g2 += square(g(i, t));
            }
        }
    }

    /**
     * Basic mc simulation.
     */
    template<typename Fn> 
    void base_mc(int steps, Fn &&ctrl_tempr) {
        DgMove mv;
        int cycle_steps = 10, len = c.rows(), local_succ_num = 0;
        double en = total_energy(c).sum();
        double min_en = en;
        auto min_coord = c;
        for (int i = 0; i < steps; i++) {
            auto index = mv.pick(len);
            auto en_atom_old = atom_energy(c, index);
            mv.move(c);
            auto en_atom_new = atom_energy(c, index);
            if (en_atom_new > en_atom_old && rand() > exp(-(en_atom_new - en_atom_old) / mc_tempr)) {
                mv.back(c);
            } else {
                en = en - en_atom_old + en_atom_new;
                local_succ_num++;
                if (en < min_en) (min_en = en, min_coord = c);
            }

            if (i % cycle_steps == cycle_steps - 1) {
                double local_succ_rate = double(local_succ_num) / cycle_steps;
                local_succ_num = 0;
                log << format("Step: %05d, temperature: %8.3f, success rate: %8.3f, energy: %8.3f", i+1, mc_tempr, local_succ_rate, en) << std::endl;
                if (! ctrl_tempr(local_succ_rate)) break;
            }
        }
        c = min_coord;
        gradient();
    }

    void heat() {
        log << "Step 1: heating..." << std::endl;
        base_mc(50, [&](double rate){if (rate < 0.8) {mc_tempr *= 2; return true;} else return false;});
    }

    void cool() {
        log << "Step 2: cooling..." << std::endl;
        base_mc(num_mc_steps, [&](double rate){
            if (rate > 0.5) {
                mc_tempr *= 0.9; 
            } else if (rate > 0.2) {
                mc_tempr *= 0.999;
            } else if (rate > 0.005) {
                mc_tempr *= 0.99999;
            }
            return true;
        });
    }

    void mc() {
        log << "Start MC.\n" << std::endl;

        heat();
        cool();

        log << "Finish MC.\n" << std::endl;
    }

    DgEnergy total_energy(const Mat &coord) {
        return {
            total_bound_energy(coord), total_distance_energy(coord), total_dihedral_energy(coord)
        };
    }

    double total_bound_energy(const Mat &coord) {
        double en = 0;
        for (int i = 0; i < coord.rows(); i++) {
            for (int j = i + 1; j < coord.rows(); j++) {
                double dist = geom::distance(coord.row(i), coord.row(j));
                en += relative_distance_energy(dist, default_range[0], default_range[1]);
            }
        }
        return en;
    }

    double total_distance_energy(const Mat &coord) {
        double en = 0;
        for (int i = 0; i < distance_constraints.size(); i++) {
            en += distance_constraint_energy(coord, i); 
        }
        return en;
    }

    double total_dihedral_energy(const Mat &coord) {
        double en = 0; 
        for (int i = 0; i < dihedral_constraints.size(); i++) {
            en += dihedral_constraint_energy(coord, i);
        }
        return en;
    }

    double atom_energy(const Mat &coord, int n) {
        return atom_bound_energy(coord, n) + atom_distance_energy(coord, n) + atom_dihedral_energy(coord, n);
    }

    double atom_bound_energy(const Mat &coord, int i) {
        double en = 0;
        for (int j = 0; j < coord.rows(); j++) {
            if (i != j) {
                double dist = geom::distance(coord.row(i), coord.row(j));
                en += relative_distance_energy(dist, default_range[0], default_range[1]);
            }
        }
        return en;
    }

    double atom_distance_energy(const Mat &coord, int n) {
        double en = 0;
        for (auto && i : involved_distance_constraints[n]) {
            en += distance_constraint_energy(coord, i);
        }
        return en;
    }

    double atom_dihedral_energy(const Mat &coord, int n) {
        double en = 0;
        for (auto && i : involved_dihedral_constraints[n]) {
            en += dihedral_constraint_energy(coord, i);
        }
        return en;
    }

    double distance_constraint_energy(const Mat &coord, int i) {
        auto &ct = distance_constraints[i];
        double dist = geom::distance(coord.row(ct.key[0]), coord.row(ct.key[1]));
        return relative_distance_energy(dist, ct.val[0], ct.val[1]);
    }

    double dihedral_constraint_energy(const Mat &coord, int i) {
        auto &ct = dihedral_constraints[i];
        double dih = geom::dihedral(coord.row(ct.key[0]), coord.row(ct.key[1]), coord.row(ct.key[2]), coord.row(ct.key[3]));
        return relative_dihedral_energy(dih, ct.val[0]);
    }

    void shrink_dihedral(double &n) {
        while (true) {
            if (n >= PI) n -= 2 * PI;
            else if (n < -PI) n += 2 * PI;
            else if (n < 0) n = -n;
            else break;
        }
    }

    double relative_distance_energy(double dist, double l, double u) {
        if (dist > u) return square(dist - u);
        else if (dist < l) return square(l - dist);
        else return 0;
    }

    double relative_dihedral_energy(double dih, double std) {
        double d = dih - std;
        shrink_dihedral(d);
        return 10 * d * d;
    }

};

}

