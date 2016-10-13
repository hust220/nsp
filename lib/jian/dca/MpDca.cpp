#include <algorithm>
#include "MpDca.hpp"

namespace jian {
namespace dca {

REG_DCA_FAC("mp", MpDca);

MpDca::MpDca() : Dca() {
    m_step_size = 0.01;
}

void MpDca::set_step(float n) {
    m_step_size = n;
}

void MpDca::init_val() {
    int i, ai;
    float v_sum;
    Vecf v(q);

    std::cout << "init: " << std::endl;
    std::cout << M << ' ' << q << std::endl;
    hi = Matf::Ones(M, q);
    for (i = 0; i < M; i++) {
        v_sum = 0;
        for (ai = 0; ai < q; ai++) {
            v[ai] = std::log(Pi(i, ai));
            v_sum += v[ai];
        }
        for (ai = 0; ai < q; ai++) {
            hi(i, ai) = v[ai] - v_sum / q;
        }
    }
    eij = Mat4::Zero(M, M, q, q);
    pi = Matf::Ones(M, q);
    pij = Mat4::Ones(M, M, q, q);
    mij = Mat3::Constant(M, M, q, 1.0/q);
    mijk = Mat5::Zero(M, M, M, q, q);
}

void MpDca::cal_pi() {
    int i, k, ai, ak;
    float v_sum, sum, prod;
    Vecf v(q);

    for (i = 0; i < M; i++) {
        v_sum = 0;
        for (ai = 0; ai < q; ai++) {
            prod = std::exp(hi(i, ai));
            for (k = 0; k < M; k++) {
                if (k != i) {
                    sum = 0;
                    for (ak = 0; ak < q; ak++) {
                        sum += std::exp(-eij(k, i, ak, ai)) * mij(k, i, ak);
                    }
                    prod *= sum;
                }
            }
            v[ai] = prod;
            v_sum += v[ai];
        }
        for (ai = 0; ai < q; ai++) {
            pi(i, ai) = v[ai] / v_sum;
        }
    }
}

void MpDca::cal_pij() {
    int i, j, l, ai, aj, al;
    float v_sum, pi_sum, sum, sum1, sum2;
    Vecf v(q);

    for (i = 0; i < M; i++) for (j = 0; j < M; j++) for (aj = 0; aj < q; aj++) {
        v_sum = 0;
        pi_sum = 0;
        for (ai = 0; ai < q; ai++) {
            sum = 0;
            for (l = 0; l < M; l++) {
                if (l != i) {
                    sum1 = 0;
                    sum2 = 0;
                    for (al = 0; al < q; al++) {
                        sum1 += std::exp(-eij(l, i, al, ai)) * mijk(l, i, j, al, aj);
                        sum2 += std::exp(-eij(l, i, al, ai)) * mij(l, i, al);
                    }
                    sum += sum1 / sum2;
                }
            }
            if (i == j && ai == aj) sum += 1;
            v[ai] = sum * Pi(i, ai);
            v_sum += v[ai];
            pi_sum += Pi(i, ai);
        }
        for (ai = 0; ai < q; ai++) {
            v[ai] -= Pi(i, ai) / pi_sum * v_sum;
            pij(i, j, ai, aj) = v[ai] + Pi(i, ai) * Pi(j, aj);
        }
    }
}

void MpDca::solve_pi() {
    int i, j, k, ai, aj, ak;
    float v_sum, sum, prod, diff, d;
    Vecf v(q);
    std::vector<int> ls1(M), ls2(M);

    std::iota(ls1.begin(), ls1.end(), 0);
    std::iota(ls2.begin(), ls2.end(), 0);

    do {
        diff = 0;
        std::random_shuffle(ls1.begin(), ls1.end());
        std::random_shuffle(ls2.begin(), ls2.end());
        for (auto && i : ls1) for (auto && j : ls2) {
            v_sum = 0;
            for (ai = 0; ai < q; ai++) {
                sum = 0;
                for (aj = 0; aj < q; aj++) {
                    sum += std::exp(-eij(i, j, ai, aj)) * mij(j, i, aj);
                }
                v[ai] = Pi(i, ai) / sum;
                v_sum += v[ai];
            }
            for (ai = 0; ai < q; ai++) {
                v[ai] /= v_sum;
                d = std::fabs(mij(i, j, ai) - v[ai]);
                if (d > diff) diff = d;
                mij(i, j, ai) = v[ai];
            }
        }
        /*
        for (auto && i : ls1) for (auto && j : ls2) {
            v_sum = 0;
            for (ai = 0; ai < q; ai++) {
                prod = std::exp(hi(i, ai));
                for (k = 0; k < M; k++) {
                    if (k != i && k != j) {
                        sum = 0;
                        for (ak = 0; ak < q; ak++) {
                            sum += std::exp(-eij(k, i, ak, ai)) * mij(k, i, ak);
                        }
                        prod *= sum;
                    }
                }
                v[ai] = prod;
                v_sum += prod;
            }
            for (ai = 0; ai < q; ai++) {
                v[ai] /= v_sum;
                d = std::fabs(mij(i, j, ai) - v[ai]);
                if (d > diff) diff = d;
                mij(i, j, ai) = v[ai];
            }
        }
        */
        std::cout << "solve pi: " << diff << std::endl;
    } while (diff > 0.001);
    cal_pi();
}

void MpDca::solve_pij() {
    int i, j, k, l, ai, ak, al;
    float v_sum, mij_sum, sum, sum1, sum2, diff, d;
    Vecf v(q);
    std::vector<int> ls1(M), ls2(M), ls3(M);

    std::iota(ls1.begin(), ls1.end(), 0);
    std::iota(ls2.begin(), ls2.end(), 0);
    std::iota(ls3.begin(), ls3.end(), 0);

    do {
        diff = 0;
        std::random_shuffle(ls1.begin(), ls1.end());
        std::random_shuffle(ls2.begin(), ls2.end());
        std::random_shuffle(ls3.begin(), ls3.end());
        for (auto && i : ls1) for (auto && j : ls2) for (auto && k : ls3) for (ak = 0; ak < q; ak++) {
        //for (i = 0; i < M; i++) for (j = 0; j < M; j++) for (k = 0; k < M; k++) for (ak = 0; ak < q; ak++) {
            v_sum = 0;
            mij_sum = 0;
            for (ai = 0; ai < q; ai++) {
                sum = 0;
                for (l = 0; l < M; l++) {
                    if (l != i && l != j) {
                        sum1 = 0;
                        sum2 = 0;
                        for (al = 0; al < q; al++) {
                            sum1 += std::exp(-eij(l, i, al, ai)) * mijk(l, i, k, al, ak);
                            sum2 += std::exp(-eij(l, i, al, ai)) * mij(l, i, al);
                        }
                        sum += sum1 / sum2;
                    }
                }
                if (i == k && ai == ak) sum += 1;
                v[ai] = mij(i, j, ai) * sum;
                v_sum += v[ai];
                mij_sum += mij(i, j, ai);
            }
            for (ai = 0; ai < q; ai++) {
                v[ai] -= v_sum * mij(i, j, ai) / mij_sum;
                d = std::fabs(v[ai] - mijk(i, j, k, ai, ak));
                if (d > diff) diff = d;
                mijk(i, j, k, ai, ak) = v[ai];
            }
        }
        std::cout << "solve pij: " << diff << std::endl;
    } while (diff > 0.001);
    cal_pij();
}

void MpDca::update_eij(float &diff) {
    int i, j, ai, aj;
    float d;

    for (i = 0; i < M; i++) for (j = 0; j < M; j++) for (ai = 0; ai < q; ai++) for (aj = 0; aj < q; aj++) {
        //d = m_step_size * (Pij(i, j, ai, aj) - pij(i, j, ai, aj) - (Pi(i, ai) + Pi(j, aj) - pi(i, ai) - pi(j, aj)) / q);
        d = m_step_size * (Pij(i, j, ai, aj) - pij(i, j, ai, aj));
        if (std::fabs(d) > diff) diff = std::fabs(d);
        eij(i, j, ai, aj) += d;
    }
}

void MpDca::update_hi(float &diff) {
    int i, j, ai, aj;
    float prod, sum, v_sum, d;
    Vecf v(q);

    for (i = 0; i < M; i++) for (ai = 0; ai < q; ai++) {
        d = m_step_size * (Pi(i, ai) - pi(i, ai));
        if (std::fabs(d) > diff) diff = std::fabs(d);
        hi(i, ai) += d;
    }

    /*
    for (i = 0; i < M; i++) {
        v_sum = 0;
        for (ai = 0; ai < q; ai++) {
            prod = 1;
            for (j = 0; j < M; j++) {
                if (j != i) {
                    sum = 0;
                    for (aj = 0; aj < q; aj++) {
                        sum += eij(i, j, ai, aj) * mij(j, i, aj);
                    }
                    prod *= sum;
                }
            }
            v[ai] = prod;
            v_sum += prod;
        }
        for (ai = 0; ai < q; ai++) {
            hi(i, ai) = v[ai] - v_sum / q;
        }
    }
    */
}

void MpDca::calculate_eij() {
    int i, j, a, b;
    float diff;

    init_val();
    i = 0;
    std::cout << "step: " << m_step_size << std::endl;
    std::cout << "Pi: " << std::endl;
    std::cout << Pi << std::endl;
    std::cout << "Pij: " << std::endl;
    Pij.print();
    std::cout << "pi: " << std::endl;
    std::cout << pi << std::endl;
    std::cout << "pij: " << std::endl;
    pij.print();
    do {
        diff = 0;
        solve_pi();
        solve_pij();
        std::cout << "mij: " << std::endl;
        mij.print();
        std::cout << "pi: " << std::endl;
        std::cout << pi << std::endl;
        std::cout << "mijk: " << std::endl;
        mijk.print();
        std::cout << "pij: " << std::endl;
        pij.print();
        update_eij(diff);
        //update_hi(diff);
        std::cout << "i: " << i << ", diff: " << diff << std::endl;
        std::cout << "hi:" << std::endl;
        std::cout << hi << std::endl;
        //eij.print();
        i++;
    } while (diff > 0.001);
}

float MpDca::cal_di(int i, int j) {
    int ai, aj, a1, a2;
    float DI, sum;
    Matf Pdir(q, q);

    for (ai = 0; ai < q; ai++) for (aj = 0; aj < q; aj++) {
        sum = 0;
        for (a1 = 0; a1 < q; a1++) for (a2 = 0; a2 < q; a2++) {
            sum += mij(i, j, a1) * std::exp(-eij(i, j, a1, a2)) * mij(j, i, a2);
        }
        Pdir(ai, aj) = mij(i, j, ai) * std::exp(-eij(i, j, ai, aj)) * mij(j, i, aj) / sum;
    }
    DI = 0;
    for (ai = 0; ai < q; ai++) for (aj = 0; aj < q; aj++) {
        DI += Pdir(ai, aj) * std::log(Pdir(ai, aj) / Pi(i, ai) / Pi(j, aj));
    }
    return DI;
}

} // namespace dca
} // namespace jian


