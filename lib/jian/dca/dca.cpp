#include <utility>
#include <algorithm>
#include <list>
#include <tuple>
#include <map>
#include <deque>
#include <string>
#include <jian/utils/Par.hpp>
#include <jian/utils/file.hpp>
#include <jian/utils/log.hpp>
#include <jian/matrix.hpp>
#include "Mat4.hpp"
#include "dca.hpp"

namespace jian {
namespace dca {

class DCA {
public:
    using seqs_t = std::deque<std::string>;

    int N, M, q;

    Mat align, seqids, fi, Pi, C, eij, DI;
    Vec ma;
    Mat4 fij, Pij;
    double theta = 0.8;
    double Meff;
    double pseudocount_weight = 0.5;
    std::string _seq;
    std::vector<int> _indices;

    //std::map<char, int> rnaNum = {{'-', 0}, {'A', 1}, {'U', 2}, {'G', 3}, {'C', 4}};
    std::vector<char> _symbols {'A', 'U', 'G', 'C'};
    std::map<char, int> _map_symbols;

    ~DCA() {
    }

    void fastaread(std::string fastafile, seqs_t &seqs) {
        std::string tmp;
        EACH_LINE(fastafile.c_str(),
            if (L[0] == '>') {
                seqs.push_back(tmp);
                tmp = "";
            } else {
                tmp += jian::trim_copy(L);
            }
        );
        seqs.push_back(tmp);
        seqs.pop_front();
    }

    void init(std::string Rfamfile, int n) {
        seqs_t seqs;
        char c;

        fastaread(Rfamfile, seqs);
        M = seqs[0].size();
        N = seqs.size();
        align = Mat::Zero(N, M);
        q = _symbols.size() + 1;
        for (int i = 0; i < q-1; i++) {
            _map_symbols[_symbols[i]] = i;
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                c = seqs[i][j];
                if (_map_symbols.count(c)) {
                    align(i, j) = _map_symbols[c];
                } else {
                    align(i, j) = q-1;
                }
            }
        }

        _seq = seqs[n];
        _indices.resize(M);
        int tmp = 0;
        for (int i = 0; i < M; i++) {
            char c = _seq[i];
            if (_map_symbols.count(c) && _map_symbols[c] < q-1) {
                _indices[i] = tmp;
                tmp++;
            } else {
                _indices[i] = -1;
            }
        }
    }

    double seqid(int a, int b) {
        int sum = 0;
        for (int i = 0; i < M; i++) {
            if (align(a, i) == align(b, i)) {
                sum++;
            }
        }
        return sum / double(M);
    }

    void cal_seqids() {
        int i, j;
        seqids = Mat::Zero(N, N);
        for (i = 0; i < N; i++) {
            seqids(i, i) = 1;
            for (j = i + 1; j < N; j++) {
                seqids(i, j) = seqids(j, i) = seqid(i, j);
            }
        }
    }

    void cal_ma() {
        ma = Vec::Zero(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (seqids(i, j) > theta) {
                    ma[i]++;
                }
            }
        }
    }

    void cal_meff() {
        Meff = 0;
        for (int i = 0; i < N; i++) {
            Meff += 1.0 / ma[i];
        }
    }

    void cal_fi() {
        fi = Mat::Zero(M, q);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                int a = align(i, j);
                fi(j, a) += 1.0 / ma[i];
            }
        }
    }

    void cal_fij() {
        fij = Mat4::Zero(M, M, q, q);
        int a, b;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                for (int k = 0; k < N; k++) {
                    a = align(k, i);
                    b = align(k, j);
                    fij(i, j, a, b) += 1.0 / ma[k];
                }
            }
        }
    }

    void calculate_f() {
        cal_seqids();
        cal_ma();
        cal_meff();
        for (int i = 0; i < N; i++) {
            ma[i] *= Meff;
        }
        cal_fi();
        cal_fij();
    }

    void calculate_P() {
        int i, j, k, l;

        Pi = Mat::Zero(M, q);
        for (i = 0; i < M; i++) {
            for (j = 0; j < q; j++) {
                Pi(i, j) = (1.0 - pseudocount_weight) * fi(i, j) + pseudocount_weight / q;
            }
        }
        Pij = Mat4::Zero(M, M, q, q);
        for (i = 0; i < M; i++) for (j = 0; j < M; j++) for (k = 0; k < q; k++) for (l = 0; l < q; l++) {
            Pij(i, j, k, l) = (1.0 - pseudocount_weight) * fij(i, j, k, l);
            if (i != j) {
                Pij(i, j, k, l) += pseudocount_weight / q / q;
            } else if (k == l) {
                Pij(i, j, k, l) += pseudocount_weight / q;
            }
        }
    }

    void calculate_C() {
        int i, j, k, l;
        C = Mat::Zero(M * (q - 1), M * (q-1));
        for (i = 0; i < M; i++) for (j = 0; j < M; j++) for (k = 0; k < q-1; k++) for (l = 0; l < q-1; l++) {
            //C(i, j, k-1, l-1) = Pij(i, j, k, l) + Pi(i, k) * Pi(j, l);
            C(i*(q-1)+k, j*(q-1)+l) = Pij(i, j, k, l) - Pi(i, k) * Pi(j, l);
        }
    }

    void calculate_eij() {
        eij = C.inverse();
        for (int i = 0; i < C.rows(); i++) {
            for (int j = 0; j < C.cols(); j++) {
                eij(i, j) = std::exp(-(eij(i, j)));
            }
        }
    }

    void set_mu(const Mat &m, const Vec &pi, const Vec &pj, Vec &mu1, Vec &mu2) {
        double epsilon = 1e-4;
        double diff = 1.0;
        Vec v1, v2;
        double sum1, sum2;
        int i;
        double d;

        while (diff > epsilon) {
            diff = 0;
            v1 = m * mu2;
            v2 = m.transpose() * mu1;
            sum1 = 0;
            sum2 = 0;
            for (i = 0; i < q-1; i++) {
                v1[i] = pi[i] / v1[i];
                v2[i] = pj[i] / v2[i];
                sum1 += v1[i];
                sum2 += v2[i];
            }
            for (i = 0; i < q-1; i++) {
                v1[i] /= sum1;
                d = std::fabs(mu1[i] - v1[i]);
                if (d > diff) diff = d;
                v2[i] /= sum2;
                d = std::fabs(mu2[i] - v2[i]);
                if (d > diff) diff = d;
            }
            mu1 = v1;
            mu2 = v2;
        }
    }

    double cal_di(int i, int j) {
        // set mu
        Vec mu1 = Vec::Constant(q-1, 1.0/(q-1));
        Vec mu2 = Vec::Constant(q-1, 1.0/(q-1));
        Mat Pdir = eij.block(i*(q-1), j*(q-1), q-1, q-1);
        Vec pi = Pi.row(i);
        Vec pj = Pi.row(j);
        set_mu(Pdir, pi, pj, mu1, mu2);

        // calculate Pdir
        int a, b;
        double sum = 0;
        for (a = 0; a < q-1; a++) for (b = 0; b < q-1; b++) {
            Pdir(a, b) *= mu1[a] * mu2[b];
            sum += Pdir(a, b);
        }
        for (a = 0; a < q-1; a++) for (b = 0; b < q-1; b++) {
            Pdir(a, b) /= sum;
        }

        // calulate DI
        double DI = 0;
        for (a = 0; a < q-1; a++) for (b = 0; b < q-1; b++) {
            DI += Pdir(a, b) * std::log(Pdir(a, b) / pi[a] / pj[b]);
        }
        return DI;
    }

    void calculate_MI(mis_t &mis) {
    }

    void calculate_DI(dis_t &dis) {
        int i, j;
        DI = Mat::Zero(M, M);
        dis.clear();
        double d;
        for (i = 0; i < M; i++) {
            for (j = i + 1; j < M; j++) {
                d = cal_di(i, j);
                dis.push_back({i, j, d});
            }
        }
        dis.sort([](auto && p1, auto && p2){
            return p1.di > p2.di;
        });
    }

    void del_gaps(dis_t &dis) {
        dis_t ls;
        for (auto && p : dis) {
            int a = _indices[p.n1];
            int b = _indices[p.n2];
            if (a != -1 && b != -1) {
                ls.push_back({a, b, p.di});
            }
        }
        dis = ls;
    }
    
    void analyze(std::string input, int n, result_t &rt) {
        LOGI << "Load FASTA alignment data ..." << input << std::endl;
        init(input, n);
        LOGI << "M = " << M << ", N = " << N << ", q = " << q << std::endl;

        LOGI << "Calculate fi, fij ... " << std::endl;
        calculate_f();
        LOGI << "Meff: " << Meff << std::endl;
        LOGI << "fi: " << std::endl;
        LOGI << fi << std::endl;

        LOGI << "Calculate Pi, Pij ..." << std::endl;
        calculate_P();
        LOGI << "Pi: " << std::endl;
        LOGI << Pi << std::endl;

        LOGI << "Calculate C ..." << std::endl;
        calculate_C();
        LOGI << "C: " << std::endl;
        LOGI << C << std::endl;

        LOGI << "Calculate eij ..." << std::endl;
        calculate_eij();
//        LOGI << "eij: " << std::endl;
//        LOGI << eij << std::endl;
    
        LOGI << "Calculate DI ..." << std::endl;
        LOGI << "DI: " << std::endl;
        calculate_MI(rt.mis);
        calculate_DI(rt.dis);
        del_gaps(rt.dis);
    
    }

    static DCA &instance() {
        static DCA dca;
        return dca;
    }

};

void analyze(std::string file, int n, result_t &rt) {
    DCA::instance().analyze(file, n, rt);
}

void print_mis(const mis_t &mis) {
    for (auto && p : mis) {
        LOGI << p.n1 << ' ' << p.n2 << ' ' << p.mi << std::endl;
    }
}

void print_dis(const dis_t &dis) {
    for (auto && p : dis) {
        LOGI << p.n1 << ' ' << p.n2 << ' ' << p.di << std::endl;
    }
}

} // namespace analyze
} // namespace jian


