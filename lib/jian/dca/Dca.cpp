#include "Dca.hpp"

namespace jian {
namespace dca {

Dca::~Dca() {}

void Dca::fastaread(std::string fastafile, Dca::seqs_t &seqs) {
    std::string tmp;
	BEGIN_READ_FILE(fastafile, " ") {
		if (L[0] == '>') {
			seqs.push_back(tmp);
			tmp = "";
		}
		else {
			tmp += jian::trim_copy(L);
		}
	} END_READ_FILE;
    seqs.push_back(tmp);
    seqs.pop_front();
}

void Dca::trim_seqs(Dca::seqs_t &seqs, int n) {
    int i = 0;
    int a = -1, b = -1;
    char c_old = '-';
    int l = seqs[n].size();
    for (auto && c : seqs[n]) {
        if (a == -1 && c != '-') {
            a = i;
        }
        if (c_old != '-' && c == '-') {
            b = i-1;
        } else if (c != '-') {
            b = -1;
        }
        c_old = c;
        i++;
    }
    if (b == -1) b = l - 1;
    for (auto && seq : seqs) {
        seq = seq.substr(a, b - a + 1);
    }
}

void Dca::init(std::string Rfamfile, int n) {
    seqs_t seqs;
    char c;
	unsigned i, j;

    fastaread(Rfamfile, seqs);
    trim_seqs(seqs, n);
    M = seqs[0].size();
    N = seqs.size();
    align = Mati::Zero(N, M);
    q = m_symbols.size() + 1;
    for (i = 0; i < q-1; i++) {
        m_map_symbols[m_symbols[i]] = i;
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            c = seqs[i][j];
            if (m_map_symbols.count(c)) {
                align(i, j) = m_map_symbols[c];
            } else {
                align(i, j) = q-1;
            }
        }
    }

    m_seq = seqs[n];
    m_indices.resize(M);
    int tmp = 0;
    for (i = 0; i < M; i++) {
        char c = m_seq[i];
        if (m_map_symbols.count(c) && m_map_symbols[c] < q-1) {
            m_indices[i] = tmp;
            tmp++;
        } else {
            m_indices[i] = -1;
        }
    }
}

float Dca::seqid(int a, int b) {
    int sum = 0;
    for (int i = 0; i < M; i++) {
        if (align(a, i) == align(b, i)) {
            sum++;
        }
    }
    return sum / float(M);
}

void Dca::cal_seqids() {
    int i, j;
    seqids = Matf::Zero(N, N);
    for (i = 0; i < N; i++) {
        seqids(i, i) = 1;
        for (j = i + 1; j < N; j++) {
            seqids(i, j) = seqids(j, i) = seqid(i, j);
        }
    }
}

void Dca::cal_ma() {
	unsigned i, j;
    ma = Vecf::Zero(N);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (seqids(i, j) > theta) {
                ma[i]++;
            }
        }
    }
}

void Dca::cal_meff() {
    Meff = 0;
    for (int i = 0; i < N; i++) {
        Meff += 1.0f / ma[i];
    }
}

void Dca::cal_fi() {
    fi = Matf::Zero(M, q);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            int a = align(i, j);
            fi(j, a) += 1.0f / ma[i];
        }
    }
}

void Dca::cal_fij() {
    fij = Mat4::Zero(M, M, q, q);
    int a, b;
	int i, j, k;
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            for (k = 0; k < N; k++) {
                a = align(k, i);
                b = align(k, j);
                fij(i, j, a, b) += 1.0f / ma[k];
            }
        }
    }
}

void Dca::calculate_f() {
    cal_seqids();
    cal_ma();
    cal_meff();
    for (int i = 0; i < N; i++) {
        ma[i] *= Meff;
    }
    cal_fi();
    cal_fij();
}

void Dca::calculate_P() {
	int i, j, k, l;

    Pi = Matf::Zero(M, q);
    for (i = 0; i < M; i++) {
        for (j = 0; j < q; j++) {
            Pi(i, j) = (1.0f - pseudocount_weight) * fi(i, j) + pseudocount_weight / q;
        }
    }
    Pij = Mat4::Zero(M, M, q, q);
    for (i = 0; i < M; i++) for (j = 0; j < M; j++) for (k = 0; k < q; k++) for (l = 0; l < q; l++) {
        Pij(i, j, k, l) = (1.0f - pseudocount_weight) * fij(i, j, k, l);
        if (i != j) {
            Pij(i, j, k, l) += pseudocount_weight / q / q;
        } else if (k == l) {
            Pij(i, j, k, l) += pseudocount_weight / q;
        }
    }
}

void Dca::calculate_DI(std::string out_file) {
    LOGI << "Calculate eij ..." << std::endl;
    calculate_eij();

    std::ofstream ofile(out_file.c_str());
	int i, j;
	int a, b;
    DI = Matf::Zero(M, M);
    for (i = 0; i < M; i++) {
        a = m_indices[i];
        if (a != -1) {
            for (j = i + 1; j < M; j++) {
                b = m_indices[j];
                if (b != -1) {
                    ofile << a+1 << ' ' << b+1 << ' ' << cal_di(i, j) << std::endl;
                }
            }
        }
    }
    ofile.close();
}

Dca &Dca::run(std::string input, std::string out_file, int n) {
    LOGI << "Load FASTA alignment file: " << input << std::endl;
    init(input, n);
    LOGI << "M = " << M << ", N = " << N << ", q = " << q << std::endl;

    LOGI << "Calculate fi, fij ... " << std::endl;
    calculate_f();
    LOGI << "Meff: " << Meff << std::endl;

    LOGI << "Calculate Pi, Pij ..." << std::endl;
    calculate_P();

    LOGI << "Calculate DI ..." << std::endl;
    calculate_DI(out_file);

    return *this;
}

} // namespace dca
} // namespace jian


