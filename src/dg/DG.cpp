#include "DG.h"

namespace jian {

MatrixXf DG::operator ()() {
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

void DG::smooth() {
    /// check minimum distance
    for (int i = 0; i < bound.rows(); i++) {
        for (int j = 0; j < bound.cols(); j++) {
            if (i == j) continue;
            if (bound(i, j) < _min_dist) {
                bound(i, j) = _min_dist;
            }
        }
    }

    /// smooth
    int flag = 0, temp = 999;
    int n = 0;
    while (temp && n < 5) {
        flag = 0;
        for (int i = 0; i < len - 2; i++) {
            for (int j = i + 1; j < len - 1; j++) {
                for (int k = j + 1; k < len; k++) {
                    /** max **/
                    double max;
                    /* i, j */
                    max = at(bound, j, k, len) + at(bound, i, k, len);
                    if (at(bound, i, j, len) > max && at(bound, j, i, len) < max) {
                        assign(bound, i, j, max, len);
                        flag++;
                    }
                    /* j, k */
                    max = at(bound, i, j, len) + at(bound, i, k, len);
                    if (at(bound, j, k, len) > max && at(bound, k, j, len) < max) {
                        assign(bound, j, k, max, len);
                        flag++;
                    }
                    /* i, k */
                    max = at(bound, i, j, len) + at(bound, j, k, len);
                    if (at(bound, i, k, len) > max && at(bound, k, i, len) < max) {
                        assign(bound, i, k, max, len);
                        flag++;
                    }

                    /** min **/
                    /* i, j */
                    double min = 0;
                    double min1 = at(bound, k, j, len) - at(bound, i, k, len);
                    double min2 = at(bound, k, i, len) - at(bound, j, k, len);
                    if (min1 > 0) {
                        min = min1;
                    } else if (min2 > 0) {
                        min = min2;
                    }
                    if (at(bound, j, i, len) < min && at(bound, i, j, len) > min) {
                        assign(bound, j, i, min, len);
                        flag++;
                    }
                    /* j, k */
                    min = 0;
                    min1 = at(bound, j, i, len) - at(bound, i, k, len);
                    min2 = at(bound, k, i, len) - at(bound, i, j, len);
                    if (min1 > 0) {
                        min = min1;
                    } else if (min2 > 0) {
                        min = min2;
                    }
                    if (at(bound, k, j, len) < min && at(bound, j, k, len) > min) {
                        assign(bound, k, j, min, len);
                        flag++;
                    }
                    /* i, k */
                    min = 0;
                    min1 = at(bound, j, i, len) - at(bound, j, k, len);
                    min2 = at(bound, k, j, len) - at(bound, i, j, len);
                    if (min1 > 0) {
                        min = min1;
                    } else if (min2 > 0) {
                        min = min2;
                    }
                    if (at(bound, k, i, len) < min && at(bound, i, k, len) > min) {
                        assign(bound, k, i, min, len);
                        flag++;
                    }
                }
            }
        }
        if (flag == temp) {
            n++;
        }
        temp = flag;
    }

}

void DG::b2d() {
    /* randomly assign values to the boundance matrix */
    d.resize(len, len);
    for (int i = 0; i < len; i++) {
        for (int j = i; j < len; j++) {
            if (j == i) {
                d(i, j) = 0;
            } else {
                d(i, j) = unif_distr(rd) * (bound(i, j) - bound(j, i)) + bound(j, i);
                d(j, i) = d(i, j);
            }
        }
    }
}

void DG::metric() {
    /* set metric matrix */
    double temp = 0;
    for (int j = 0; j < len; j++) {
        for (int k = j + 1; k < len; k++) {
            temp += d(j, k) * d(j, k);
        }
    }
    temp /= double(len * len);
    double a[len];
    for (int i = 0; i < len; i++) {
        a[i] = 0;
        for (int j = 0; j < len; j++) {
            a[i] += d(i, j) * d(i, j);
        }
        a[i] /= double(len);
        a[i] -= temp;
    }
    m.resize(len, len);
    for (int i = 0; i < len; i++) {
        for (int j = i; j < len; j++) {
            m(i, j) = (a[i] + a[j] - d(i, j) * d(i, j)) / 2.;
            m(j, i) = m(i, j);
        }
    }
}

void DG::d2c() {
    /* transform bounds matrix to coordinates matrix */
    MatrixXf A(len, len);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            A(i, j) = m(i, j);
        }
    }
    SelfAdjointEigenSolver<MatrixXf> eigensolver(A);
//    cout << eigensolver.eigenvalues() << endl;
//    cout << eigensolver.eigenvectors() << endl;
    double max1 = 0, max2 = 0, max3 = 0;
    int m1 = 0, m2 = 0, m3 = 0;
    for (int i = 0; i < len; i++) {
        double temp = abs(eigensolver.eigenvalues()[i]);
        if (temp > max1) {
            max3 = max2;
            m3 = m2;
            max2 = max1;
            m2 = m1;
            max1 = temp;
            m1 = i;
        } else if (temp > max2) {
            max3 = max2;
            m3 = m2;
            max2 = temp;
            m2 = i;
        } else if (temp > max3) {
            max3 = temp;
            m3 = i;
        }
    }
    c.resize(len, 3);
    //MatrixXf B(len, 3);
    for (int i = 0; i < len; i++) {
        double temp = eigensolver.eigenvalues()[m1];
        //double temp = eigensolver.eigenvalues()[m1] * eigensolver.eigenvectors()(i, m1);
        if (temp > 0) {
            temp = sqrt(temp);
        } else if (temp < 0) {
            temp = -sqrt(-temp);
        }
        c(i, 0) = temp * eigensolver.eigenvectors()(i, m1);
        //c(i, 0) = temp;
        //B(i, len - 1 - j) = temp;
        temp = eigensolver.eigenvalues()[m2];
        //temp = eigensolver.eigenvalues()[m2] * eigensolver.eigenvectors()(i, m2);
        if (temp > 0) {
            temp = sqrt(temp);
        } else if (temp < 0) {
            temp = -sqrt(-temp);
        }
        c(i, 1) = temp * eigensolver.eigenvectors()(i, m2);
        //c(i, 1) = temp;
        //B(i, len - 1 - j) = temp;
        temp = eigensolver.eigenvalues()[m3];
        //temp = eigensolver.eigenvalues()[m3] * eigensolver.eigenvectors()(i, m3);
        if (temp > 0) {
            temp = sqrt(temp);
        } else if (temp < 0) {
            temp = -sqrt(-temp);
        }
        c(i, 2) = temp * eigensolver.eigenvectors()(i, m3);
    }
}

void DG::c2d() {
    for (int i = 0; i < len; i++) {
        for (int j = i; j < len; j++) {
            if (i == j) {
                d(i, j) = 0;
            } else {
                d(i, j) = sqrt((c(i, 0) - c(j, 0)) * (c(i, 0) - c(j, 0)) + (c(i, 1) - c(j, 1)) * (c(i, 1) - c(j, 1)) + (c(i, 2) - c(j, 2)) * (c(i, 2) - c(j, 2)));
                d(j, i) = d(i, j);
            }
        }
    }
}

double DG::total_energy(const MatrixXf &coord) {
    double en = 0;
    double x, y, z, u2, l2, dist2;
    for (int i = 0; i < coord.rows(); i++) {
        for (int j = i + 1; j < coord.rows(); j++) {
            x = coord(i, 0) - coord(j, 0);
            y = coord(i, 1) - coord(j, 1);
            z = coord(i, 2) - coord(j, 2);
            dist2 = x * x + y * y + z * z;
            u2 = bound(i, j) * bound(i, j);
            l2 = bound(j, i) * bound(j, i);
            if (dist2 > u2) {
                en += (dist2 - u2) * (dist2 - u2) / u2 / u2;
            } else if (dist2 < l2) {
                en += (dist2 - l2) * (dist2 - l2) / dist2 / dist2;
            }
        }
    }
    return en;
}

double DG::atom_energy(const MatrixXf &coord, int n) {
    double en = 0;
    double x, y, z, u2, l2, dist2;
    for (int i = 0; i < coord.rows(); i++) {
        if (i != n) {
            x = coord(i, 0) - coord(n, 0);
            y = coord(i, 1) - coord(n, 1);
            z = coord(i, 2) - coord(n, 2);
            dist2 = x * x + y * y + z * z;
            if (i < n) {
                u2 = bound(i, n) * bound(i, n);
                l2 = bound(n, i) * bound(n, i);
            } else {
                u2 = bound(n, i) * bound(n, i);
                l2 = bound(i, n) * bound(i, n);
            }
            if (dist2 > u2) {
                en += (dist2 - u2) * (dist2 - u2) / u2 / u2;
            } else if (dist2 < l2) {
                en += (dist2 - l2) * (dist2 - l2) / dist2 / dist2;
            }
        }
    }
    return en;
}

void DG::gradient() {
    E = 0;
    MatrixXf C(3 * len, 3);
    double err = 1.e-6;
    for (int i = 0; i < len; i++) {
        for (int t = 0; t < 3; t++) {
            C(3 * i + t, 0) = 0;
            C(3 * i + t, 1) = 0;
            C(3 * i + t, 2) = 0;
        }
        for (int j = 0; j < len; j++) {
            if (i == j) continue;
            double d2_l, d2, d2_r, u2, l2;
            /* d2 */
            d2 = (c(i, 0) - c(j, 0)) * (c(i, 0) - c(j, 0)) + (c(i, 1) - c(j, 1)) * (c(i, 1) - c(j, 1)) + (c(i, 2) - c(j, 2)) * (c(i, 2) - c(j, 2));
            if (i < j) {
                u2 = bound(i, j) * bound(i, j);
                l2 = bound(j, i) * bound(j, i);
            } else{
                u2 = bound(j, i) * bound(j, i);
                l2 = bound(i, j) * bound(i, j);
            }
            for (int t = 0; t < 3; t++) {
                /* d2_l */
                if (t == 0) {
                    d2_l = (c(i, 0) - err - c(j, 0)) * (c(i, 0) - err - c(j, 0)) + (c(i, 1) - c(j, 1)) * (c(i, 1) - c(j, 1)) + (c(i, 2) - c(j, 2)) * (c(i, 2) - c(j, 2));
                } else if (t == 1) {
                    d2_l = (c(i, 0) - c(j, 0)) * (c(i, 0) - c(j, 0)) + (c(i, 1) - err - c(j, 1)) * (c(i, 1) - err - c(j, 1)) + (c(i, 2) - c(j, 2)) * (c(i, 2) - c(j, 2));
                } else {
                    d2_l = (c(i, 0) - c(j, 0)) * (c(i, 0) - c(j, 0)) + (c(i, 1) - c(j, 1)) * (c(i, 1) - c(j, 1)) + (c(i, 2) - err - c(j, 2)) * (c(i, 2) - err - c(j, 2));
                }
                /* d2_r */
                if (t == 0) {
                    d2_r = (c(i, 0) + err - c(j, 0)) * (c(i, 0) + err - c(j, 0)) + (c(i, 1) - c(j, 1)) * (c(i, 1) - c(j, 1)) + (c(i, 2) - c(j, 2)) * (c(i, 2) - c(j, 2));
                } else if (t == 1) {
                    d2_r = (c(i, 0) - c(j, 0)) * (c(i, 0) - c(j, 0)) + (c(i, 1) + err - c(j, 1)) * (c(i, 1) + err - c(j, 1)) + (c(i, 2) - c(j, 2)) * (c(i, 2) - c(j, 2));
                } else {
                    d2_r = (c(i, 0) - c(j, 0)) * (c(i, 0) - c(j, 0)) + (c(i, 1) - c(j, 1)) * (c(i, 1) - c(j, 1)) + (c(i, 2) + err - c(j, 2)) * (c(i, 2) + err - c(j, 2));
                }

                /* d2_l */
                if (d2_l > u2) {
                    C(3 * i + t, 0) += (d2_l - u2) * (d2_l - u2) / u2 / u2;
                } else if (d2_l < u2 && d2_l > l2) {
                    C(3 * i + t, 0) += 0;
                } else {
                    C(3 * i + t, 0) += (l2 - d2_l) * (l2 - d2_l) / d2_l / d2_l;
                }
                /* d2 */
                if (d2 > u2) {
                    C(3 * i + t, 1) += (d2 - u2) * (d2 - u2) / u2 / u2;
                } else if (d2 < u2 && d2 > l2) {
                    C(3 * i + t, 1) += 0;
                } else {
                    C(3 * i + t, 1) += (l2 - d2) * (l2 - d2) / d2 / d2;
                }
                /* d2_r */
                if (d2_r > u2) {
                    C(3 * i + t, 2) += (d2_r - u2) * (d2_r - u2) / u2 / u2;
                } else if (d2_r < u2 && d2_r > l2) {
                    C(3 * i + t, 2) += 0;
                } else {
                    C(3 * i + t, 2) += (l2 - d2_r) * (l2 - d2_r) / d2_r / d2_r;
                }
            }
        }
        E += 0.5 * C(3 * i + 0, 1);
    }
    
    MatrixXf ch(len * 3, 3);
    for (int i = 0; i < len * 3; i++) {
        for (int j = 0; j < 3; j++) {
            ch(i, j) = 0;
        }
    }
    CH = 0;
    if (chir.rows() != 0) {
        for (int i = 0; i < chir.rows(); i++) {
            Point p[4];
            int k[4];
            for (int j = 0; j < 4; j++) {
                k[j] = int(chir(i, j));
                p[j].x = c(k[j], 0);
                p[j].y = c(k[j], 1);
                p[j].z = c(k[j], 2);
            }
            double center = chir(i, 4);
            for (int j = 0; j < 4; j++) {
                for (int t = 0; t < 3; t++) {
                    double temp;
                    if (t == 0) {
                        p[j].x -= err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 0) += (temp - center) * (temp - center);
                        p[j].x += err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 1) += (temp - center) * (temp - center);
                        CH += 0.25 * ch(k[j] * 3 + t, 1);
                        p[j].x += err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 2) += (temp - center) * (temp - center);
                    } else if (t == 1) {
                        p[j].y -= err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 0) += (temp - center) * (temp - center);
                        p[j].y += err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 1) += (temp - center) * (temp - center);
                        CH += 0.25 * ch(k[j] * 3 + t, 1);
                        p[j].y += err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 2) += (temp - center) * (temp - center);
                    } else {
                        p[j].z -= err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 0) += (temp - center) * (temp - center);
                        p[j].z += err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 1) += (temp - center) * (temp - center);
                        CH += 0.25 * ch(k[j] * 3 + t, 1);
                        p[j].z += err;
                        temp = Point::chirality(p[0], p[1], p[2], p[3]);
                        ch(k[j] * 3 + t, 2) += (temp - center) * (temp - center);
                    }
                }
            }
        }
    }

    g.resize(len, 3);
    g2 = 0;
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < 3; j++) {
            g(i, j) = (C(3 * i + j, 2) - C(3 * i + j, 0)) / (2 * err) + (ch(3 * i + j, 2) - ch(3 * i + j, 0)) / (2 * err);
            g2 += g(i, j) * g(i, j);
        }
    }
}

void DG::cg() {
    MatrixXf d_o(len, 3);
    MatrixXf d_n(len, 3);

    if (view) std::cerr << "=============== Conjugate Gradient ===============" << std::endl;
    if (view) std::cerr << setw(8) << "step" << setw(15) << "factor" << setw(15) << "energy" << setw(15) << "chirality" << std::endl;

    double a = 0.001;
    double beta = 0;
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < 3; j++) {
            d_o(i, j) = 0;
        }
    }
    gradient();
    double oldE = E + CH;
    double oldG2 = g2;
    int upd = 0;
    while (upd < 500) {
        double tempG2 = g2;
        double tempE = E + CH;

        if (upd != 0) beta = g2 / oldG2;
        if (upd % 10 == 0) beta = 0;
        for (int i = 0; i < len; i++) {
            for (int j = 0; j < 3; j++) {
                d_n(i, j) = -g(i, j) + beta * d_o(i, j);
                c(i, j) += a * d_n(i, j);
            }
        }

        gradient();

        if (E + CH >= tempE || g2 >= tempG2) {
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < 3; j++) {
                    c(i, j) -= a * d_n(i, j);
                }
            }
            gradient();
            a *= 0.5;
            if (a < 1.e-12) break;
        } else {
            upd++;
            a *= 2;
            for (int i = 0; i < len; i++) {
                for (int j = 0; j < 3; j++) {
                    d_o(i, j) = d_n(i, j);
                }
            }
            oldG2 = tempG2;
            oldE = tempE;
            if (view) std::cerr << setw(8) << upd << setw(15) << a << setw(15) << E << setw(15) << CH << std::endl;
            if (E + CH < 1.e-12 && g2 < 1.e-12) break;
        }
    }
    if (view) std::cerr << setw(8) << upd << setw(15) << a << setw(15) << E << setw(15) << CH << std::endl;
}

/** simulated annealing monte carlo method */
void DG::mc() {
    E = total_energy(c);
    double min_en = E;
    MatrixXf min_coord = c;
    int succ_num = 0;
    double succ_rate = 1;
    int total_steps = 10000;
    int echo_steps = 1000;
    int cycle_steps = 50;
    double tempr = 20;
    int local_succ_num = 0;
    double en_atom_old = 0;
    double en_atom_new = 0;

    /* estimating an appropriate maximum temprature */
    /*
    std::cout << "Step 1: maximum temprature estimation. " << std::endl;
    for (int i = 0; i < 50; i++) {
        int m = int(unif_distr(rd) * len);
        int n = int(unif_distr(rd) * 3);
        double x = (unif_distr(rd) - 0.5) * 2;
        c(m, n) += x;
        double en = total_energy(c);
        orig_tempr += fabs(en - E);
        c(m, n) -= x;
    }
    orig_tempr /= 50;
    std::cout << "maximum temprature: " << orig_tempr << std::endl;
    */

    /* print initial condition */
    if (view) {
        std::cerr << "=============== Simulated Annealing =================" << std::endl;
        std::cerr << "Initial condition." << std::endl;
        std::cerr << "Energy: " << total_energy(c) << std::endl;
        std::cerr << "Temprature: " << tempr << std::endl;
    }

    /* warming */
    if (view) {
        std::cerr << "Step 1: warming... " << std::endl;
        std::cerr << setw(8) << "step" << setw(15) << "energy" << setw(15) << "temprature" << setw(15) << "success rate" << std::endl;
    }
    int warming_steps = 1000;
    for (int i = 0;; i++) {
        int m = int(unif_distr(rd) * len);
        int n = int(unif_distr(rd) * 3);
        double x = (unif_distr(rd) - 0.5) * 2;
        en_atom_old = atom_energy(c, m);
        c(m, n) += x;
        en_atom_new = atom_energy(c, m);
        if (en_atom_new > en_atom_old) { // randomly reject or accept
            if (unif_distr(rd) > exp(-(en_atom_new - en_atom_old) / tempr)) { // reject
                c(m, n) -= x;
            } else { // accept
                E = E - en_atom_old + en_atom_new;
                local_succ_num++;
            }
        } else { // accept
            E = E - en_atom_old + en_atom_new;
            local_succ_num++;
        }

        if (i % cycle_steps == cycle_steps - 1) {
            double local_succ_rate = double(local_succ_num) / cycle_steps;
            if (local_succ_rate < 0.8) {
                tempr *= 2;
            } else {
                if (view) std::cerr << setw(8) << i + 1 << setw(15) << E << setw(15) << tempr << setw(15) << local_succ_rate << std::endl;
                break;
            }
            local_succ_num = 0;
            //succ_rate = succ_num / (i + 1.0);
            if (view) std::cerr << setw(8) << i + 1 << setw(15) << E << setw(15) << tempr << setw(15) << local_succ_rate << std::endl;
        }
    }

    /* cooling */
    if (view) std::cerr << "Step 2: cooling... " << std::endl;
    if (view) std::cerr << setw(8) << "step" << setw(15) << "energy" << setw(15) << "temprature" << setw(15) << "success rate" << std::endl;
    for (int i = 0; i < total_steps; i++) {
        int m = int(unif_distr(rd) * len);
        int n = int(unif_distr(rd) * 3);
        double x = (unif_distr(rd) - 0.5) * 2;
        en_atom_old = atom_energy(c, m);
        c(m, n) += x;
        en_atom_new = atom_energy(c, m);
        if (en_atom_new > en_atom_old) { // randomly reject or accept
            if (unif_distr(rd) > exp(-(en_atom_new - en_atom_old) / tempr)) { // reject
                c(m, n) -= x;
            } else { // accept
                E = E - en_atom_old + en_atom_new;
                local_succ_num++;
            }
        } else { // accept
            E = E - en_atom_old + en_atom_new;
            local_succ_num++;
            
            // store minimal energy and correspondingly coordinates
            if (E < min_en) {
                min_en = E;
                min_coord = c;
            }
        }

        if (i % cycle_steps == cycle_steps - 1) {
            double local_succ_rate = double(local_succ_num) / cycle_steps;
            //if (local_succ_rate > 0.5) {
            //    tempr *= 0.5;
            //} else {
                tempr *= 0.9;
            //}
            local_succ_num = 0;
            //succ_rate = succ_num / (i + 1.0);
            if (view) std::cerr << setw(8) << i + 1 << setw(15) << E << setw(15) << tempr << setw(15) << local_succ_rate << std::endl;
        }
    }

    // pick up the coordinate corresponding to the minimal energy
    c = min_coord;
    gradient();

    if (view) std::cerr << "minimal E: " << total_energy(min_coord) << ' ' << E << std::endl;
}

double DG::at(const MatrixXf &bound, int i, int j, int len) {
    if (j >= len) {
        j -= len;
        return at(bound, j, i, len);
    }
    if (i >= len) {
        i -= len;
        return at(bound, j, i, len);
    }
    return bound(i, j);
}

void DG::assign(MatrixXf &bound, int i, int j, double d, int len) {
    if (j >= len) {
        j -= len;
        assign(bound, j, i, d, len);
        return;
    }
    if (i >= len) {
        i -= len;
        assign(bound, j, i, d, len);
        return;
    }
    bound(i, j) = d;
}

}

