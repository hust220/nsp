#pragma once

#include "pdb.hpp"
#include "rtsp_bires.hpp"
#include "rtsp_bps.hpp"
#include "score_par_bp.hpp"
#include "score.hpp"
#include "mc.hpp"
#include "geom.hpp"
#include "rand.hpp"

BEGIN_JN

Chain build_loop(Str seq, Str ss);

Chain init_loop(Str seq, Str ss);

inline auto sp_bp(const Residue &res1, const Residue &res2) {
    Mat m1, m2;
    auto cg = CG::fac_t::make("6p");
    Residue r1 = cg->to_cg(res1), r2 = cg->to_cg(res2);
    Int l1 = size(r1), l2 = size(r2);
    m1.resize(l1, 3);
    m2.resize(l2, 3);
    for (Int i = 0; i < l1; i++) {
        for (Int j = 0; j < 3; j++) {
            m1(i, j) = r1[i][j];
            m2(i, j) = r2[i][j];
        }
    }
    return geom::suppos(m1, m2);
}

inline void sample_loop(Chain &chain, Str ss, Int n) {
    Str seq = JN_ seq(chain);
    Int l = size(seq);
    Chain c = bires_chain(seq[n], seq[n+1], ss[n], ss[n+1]);
    if (n < l/2) {
        auto sp = sp_bp(c[1], chain[n+1]);
        for (Int i = 0; i < 2; i++) for (auto && atom : c[i]) sp.apply(atom);
        sp = sp_bp(chain[n], c[0]);
        for (Int i = 0; i < n; i++) for (auto && atom : chain[i]) sp.apply(atom);
        chain[n] = c[0];
        chain[n+1] = c[1];
    }
    else {
        auto sp = sp_bp(c[0], chain[n]);
        for (Int i = 0; i < 2; i++) for (auto && atom : c[i]) sp.apply(atom);
        sp = sp_bp(chain[n+1], c[1]);
        for (Int i = n+2; i < l; i++) for (auto && atom : chain[i]) sp.apply(atom);
        chain[n] = c[0];
        chain[n+1] = c[1];
    }
}

inline void sample_loop_local(Chain &chain, Str ss, Int n) {
    Str seq = JN_ seq(chain);
    Int l = size(seq);
    Chain c = bires_chain(seq[n], seq[n+1], ss[n], ss[n+1]);

    auto sp = sp_bp(c[0], chain[n]);
    for (Int i = 0; i < 2; i++) for (auto && atom : c[i]) sp.apply(atom);
    chain[n] = c[0];
    chain[n+1] = c[1];
}

inline void sample_loop(Chain &chain, Str ss) {
    Int n = Int(JN_ rand() * (size(ss)-1));
    sample_loop(chain, ss, n);
}

inline void sample_loop_local(Chain &chain, Str ss) {
    Int n = Int(JN_ rand() * (size(ss)-1));
    sample_loop_local(chain, ss, n);
}

/**
 * Build and sample loop by fragment assembly
 */
class BuildLoopFA : public MC {
public:
    struct En {
        Num crash = 0, pairing = 0, len = 0;
    };

    Chain pred;
    Str ss;
    Int ind;
    Chain restored;
    SP<Score> scorer;
    Bps bps;
    En en;
    Map<Str, Array<Num, 6>> bp_dists;
    Str traj_name;
    Bool show_log;

    BuildLoopFA() {
        scorer = Score::fac_t::make("6p");
        scorer->init();
        load_bp_dists();
        show_log = false;
    }

    void load_bp_dists() {
        for (auto && it : FileLines(to_str(Env::lib(), "/RNA/pars/nuc3d/bp_distances.txt"))) {
            if (size(it.arr) == 7) {
                for (int i = 0; i < 6; i++) {
                    bp_dists[it.arr[0]][i] = JN_ lexical_cast<Num>(it.arr[i+1]);
                }
            }
        }
    }

    void run() {
        bps = ss_to_bps(ss);
        write_traj();
//        _mc_queue = "samc:50000:20-0:0.999";
        _mc_queue = "samc:50000";
        _mc_highest_tempr = ss.size() * 3;
        _mc_lowest_tempr = 0;
        _mc_dec_rate = 0.995;
        _mc_write_steps = 100;
        mc_run();
    }

    Num en_crash(const Residue &r1, const Residue &r2) {
        return 100*scorer->en_crash(r1, r2);
    }

//    Num en_bp(const Residue &res1, const Residue &res2) {
//        auto cg = CG::fac_t::make("6p");
//        Residue r1 = cg->to_cg(res1), r2 = cg->to_cg(res2);
//
//        Str seq = to_str(res1.name.back(), res2.name.back());
//        const auto & dists = bp_dists[seq];
//
//        Num d = geom::distance(r1[0], r2[0]);
//        Num e1 = JN_ square(d - dists[0]);
//
//        d = geom::distance(r1[2], r2[2]);
//        e1 += JN_ square(d - dists[2]);
//
//        scorer->en_bp(r1, r2);
//        Num e2 = scorer->m_en_pairing;
//
//        return e1 + 20 * e2;
//
//    }

    Num en_bp(const Residue &res1, const Residue &res2) {
        auto cg = CG::fac_t::make("6p");
        Residue r1 = cg->to_cg(res1), r2 = cg->to_cg(res2);
        Str seq = to_str(res1.name.back(), res2.name.back());
        const auto & dists = bp_dists[seq];
        Num e = 0;
        for (Int i = 0; i < 6; i++) {
            Num d = geom::distance(r1[i], r2[i]);
            e += JN_ square(d - dists[i]);
        }
        return 0.1*e;

    }

    virtual double mc_total_energy() {
        En e;
        Int l = size(ss);
        for (Int i = 0; i < l; i++) {
            for (Int j = i + 1; j < l; j++) {
                e.crash += en_crash(pred[i], pred[j]);
            }
        }
        for (auto && bp : bps) {
            e.pairing += en_bp(pred[bp[0]], pred[bp[1]]);
        }
        for (Int i = 0; i < l-1; i++) {
            e.len += scorer->en_len(pred, i);
        }
        en = e;
        return e.crash + e.pairing + e.len;
    }

    virtual double mc_partial_energy() {
        return mc_total_energy();

//        Int l = size(ss);
//        En e;
//        bool left = ind<l/2;
//        bool right = !left;
//        for (Int i = 0; i < l; i++) {
//            for (Int j = i + 1; j < l; j++) {
//                if ((right && i < ind && j < ind) || (left && i > ind+1 && j > ind+1)) continue;
//                e.crash += en_crash(pred[i], pred[j]);
//            }
//        }
//        for (auto && bp : bps) {
//            int i = bp[0];
//            int j = bp[1];
//            if ((right && i < ind && j < ind) || (left && i > ind+1 && j > ind+1)) continue;
//            e.pairing += en_bp(pred[i], pred[j]);
//        }
//        for (Int i = 0; i < l-1; i++) {
//            if (left && i <= ind+1 || (right && i >= ind-1)) {
//                e.len += scorer->en_len(pred, i);
//            }
//        }
//        return e.crash + e.pairing + e.len;
    }

    virtual void mc_write() {
        mc_total_energy();
        write_traj();
        if (show_log) {
            std::cout
                << "Step: " << _mc_step + 1
                << " total:" << en.crash+en.pairing+en.len
                << " crash:" << en.crash
                << " pairing:" << en.pairing
                << " len:" << en.len
                << " tempr:" << _mc_tempr
                << std::endl;
        }
    }

    void write_traj() {
        if (traj_name.empty()) return;

        static MolWriter writer;
        static Int n = 0;

        n++;
        if (n == 1) file::clean(traj_name);
        std::ofstream output(traj_name, std::ios::app);
        writer.bind_stream(output);
        writer.write_model([this](){writer.write(pred);});
        output.close();
    }

    virtual void mc_select() {
        ind = Int(JN_ rand() * (size(ss)-1));
    }

    virtual void mc_sample() {
        sample_loop(pred, ss, ind);
    }

    virtual void mc_backup() {
        Int l = size(ss);
        restored.clear();
        if (ind < l/2) for (Int i = 0; i <= ind+1; i++) restored.push_back(pred[i]);
        else for (Int i = ind; i < l; i++) restored.push_back(pred[i]);
    }

    virtual void mc_rollback() {
        Int l = size(ss);
        if (ind < l/2) for (Int i = 0; i <= ind+1; i++) pred[i] = restored[i];
        else for (Int i = ind; i < l; i++) pred[i] = restored[i-ind];
    }
};

END_JN

