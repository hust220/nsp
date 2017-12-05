#include "rtsp_en_bp.hpp"
#include "rtsp_bires.hpp"
#include "rtsp_bires_models.hpp"
#include "cg.hpp"
#include "geom.hpp"

BEGIN_JN

static Mat mat_bp(const Residue &res1, const Residue &res2) {
    auto cg = CG::fac_t::make("6p");
    Residue r1 = cg->to_cg(res1), r2 = cg->to_cg(res2);

    Mat m(12, 3);
    for (Int i = 0; i < 6; i++) for (Int j = 0; j < 3; j++) m(i, j) = r1[i][j];
    for (Int i = 0; i < 6; i++) for (Int j = 0; j < 3; j++) m(i+6, j) = r2[i][j];
    return m;
}

struct MatsBp {
    Array<Mat, 16> mats;

    MatsBp() {
        Int i = 0;
        for (auto && a : {'A', 'U', 'G', 'C'}) {
            for (auto && b : {'A', 'U', 'G', 'C'}) {
                const Chain & c = bires_chain(a, b, '(', ')');
                mats[i] = mat_bp(c[0], c[1]);
                i++;
            }
        }
    }
};

static Array<Mat, 16> & mats_bp() {
    static MatsBp mats;
    return mats.mats;
}

static Int seq_code(Char c) {
    static Map<Char, Int> m{{'A', 0}, {'U', 1}, {'T', 1}, {'G', 2}, {'C', 3}};
    return m.at(c);
}

Num en_bp(const Residue &r1, const Residue &r2) {
    Int code = seq_code(r1.name.back())*4+seq_code(r2.name.back());
    Num d = geom::rmsd(mat_bp(r1, r2), mats_bp()[code]);
    return d*d;
}

END_JN

