#include "rand.hpp"
#include "matrix.hpp"
#include "dg_smooth.hpp"
#include "dg_move.hpp"
#include "geom.hpp"
#include "log.hpp"
#include "dg.hpp"

namespace jian {

class DG {
public:
    int len = 0;

    Log log;

    Mat bound; // distance bound matrix.
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

    void read(int len_,
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

        std::cout << bound << std::endl;

        log << "[*] DG: Set Bound." << std::endl;
        init_bound();

        std::cout << bound << std::endl;

        log << "[*] DG: Smooth." << std::endl;
        dg_smooth(bound, 3);

        std::cout << bound << std::endl;
    }

    Mat sample() {
        log << "[*] Start DG Sampling..." << std::endl;

        log << "[*] DG: Randomly build a distance matrix." << std::endl;
        b2d();

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

private:
    void init_bound() {
        bound.resize(len, len);
        for (int i = 0; i < len; i++) {
            for (int j = i; j < len; j++) {
                if (i = j) {
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
        double en_t; // total energy
        en_t = total_energy(c);
        double oldE = en_t;

        gradient();
        double oldG2 = g2;

        Mat d_o = Mat::Zero(len, 3);
        Mat d_n(len, 3);

        double a = 0.0001;
        double beta = 0;
        int upd = 0;

        while (upd < 500) {
            double tempG2 = g2, tempE = en_t;

            if (upd != 0) beta = g2 / oldG2;
            if (upd % 10 == 0) beta = 0;
            for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) {
                d_n(i, j) = -g(i, j) + beta * d_o(i, j);
                c(i, j) += a * d_n(i, j);
            }

            en_t = total_energy(c);
            gradient();

            if (en_t >= tempE || g2 >= tempG2) {
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
                if (en_t < 1.e-12 || g2 < 1.e-12) break;
            }
        }
    }

    /**
     * Randomly assign values to the boundance matrix.
     */
    void b2d() {
        d.resize(len, len);
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (j == i) d(i, j) = 0;
            else d(i, j) = d(j, i) = rand() * (bound(i, j) - bound(j, i)) + bound(j, i);
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
     * transform bounds matrix to coordinates matrix
     */
    void embed() {
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

    /*
     * Transform coordiantes to distance matrix.
     */
    Mat c2d(const Mat &coord) {
        Mat dist(coord.rows(), coord.cols());
        for (int i = 0; i < len; i++) for (int j = i; j < len; j++) {
            if (i == j) dist(i, j) = 0;
            else dist(i, j) = dist(j, i) = geom::distance(coord.row(i), coord.row(j));
        }
        return dist;
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
        double en = total_energy(c);
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
                log << i + 1 << ": " << mc_tempr << "(temprature) " << local_succ_rate << "(success rate)" << std::endl;
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

    double total_energy(const Mat &coord) {
        return total_bound_energy(coord) + total_distance_energy(coord) + total_dihedral_energy(coord);
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

SP<DG> create_dg() {
    return std::make_shared<DG>();
}

}

