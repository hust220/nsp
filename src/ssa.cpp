#include <jian/utils/log.hpp>
#include <jian/matrix.hpp>
#include <string>
#include <deque>
#include "nsp.hpp"

namespace jian {

using align_t = std::deque<std::string>;

class SSA {
public:
    int M, N;
    Mat m;
    std::string _seq1, _seq2;
    align_t _align;

    static SSA &instance() {
        static SSA ssa;
        return ssa;
    }

    void scoring() {
        int i, j;
        for (i = 1; i < N+1; i++) {
            for (j = 1; j < M+1; j++) {
                if (_seq2[i-1] == _seq1[j-1]) {
                    m(i, j) = m(i-1, j-1) + 1;
                } else {
                    m(i, j) = std::max(std::max(m(i-1,j), m(i,j-1)),m(i-1,j-1));
                }
            }
        }
    }

    void backtrack(align_t &align) {
        align.resize(2);
//        LOGI << m << std::endl;
        int i = N, j = M;
        char a, b;
        while (i != 0 || j != 0) {
            a = _seq2[i-1];
            b = _seq1[j-1];
            if (a == b) {
                align[0].insert(0,1,b);
                align[1].insert(0,1,a);
                i--;
                j--;
            } else {
                if (m(i-1,j) > m(i,j-1)) {
                    align[0].insert(0,1,'-');
                    align[1].insert(0,1,a);
                    i--;
                } else {
                    align[0].insert(0,1,b);
                    align[1].insert(0,1,'-');
                    j--;
                }
            }
        }
        for (auto && s : align) {
            LOGI << s << std::endl;
        }
    }

    void ssa(std::string seq1, std::string seq2, align_t &aln) {
        _seq1 = seq1;
        _seq2 = seq2;
        M = seq1.size();
        N = seq2.size();
        m = Mat::Zero(N+1, M+1);
        scoring();
        backtrack(aln);
    }
};

REGISTER_NSP_COMPONENT(ssa) {
    std::string seq1 = par["seq"][0];
    std::string seq2 = par["seq"][1];
    align_t aln;
    SSA::instance().ssa(seq1, seq2, aln);
}

} // namespace jian

