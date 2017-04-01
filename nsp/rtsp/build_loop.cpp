#include "build_loop.hpp"
#include "bires.hpp"
#include "bps.hpp"
//#include "en_bp.hpp"
#include "../scoring/ParBp.hpp"
#include "../scoring/Score.hpp"
#include "../mc.hpp"
#include <jian/geom.hpp>
#include <jian/utils/rand.hpp>

BEGIN_JN

static auto sp_bp(const Residue &res1, const Residue &res2) {
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

Chain init_loop(Str seq, Str ss) {
    Int n = size(seq);
    Chain chain = bires_chain(seq[0], seq[1], ss[0], ss[1]);
    //mol_write(chain, "aa-1.pdb");
    for (Int i = 1; i < n-1; i++) {
        Chain c = bires_chain(seq[i], seq[i+1], ss[i], ss[i+1]);
        //mol_write(c, to_str("aa-", i+1, ".pdb"));
        auto sp = sp_bp(c[0], chain.back());
        for (auto && atom : c[1]) sp.apply(atom);
        chain.push_back(c[1]);
    }
    return chain;
}

/*
static Bool satisfy(const Chain &chain, Str ss) {
    auto bps = ss_to_bps(ss);
    for (auto && bp : bps) std::cout << bp[0] << '-' << bp[1] << ' '; std::cout << std::endl;
    for (auto && bp : bps) {
        if (!is_bp(chain[bp[0]], chain[bp[1]])) {
            return false;
        }
    }
    return true;
}
*/

static void sample_loop(Chain &chain, Str ss, Int n) {
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

void sample_loop(Chain &chain, Str ss) {
    Int n = Int(JN_ rand() * (size(ss)-1));
    sample_loop(chain, ss, n);
}

class SampleLoop : public MC {
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

        SampleLoop() {
            scorer = Score::fac_t::make("6p");
            scorer->init();
            load_bp_dists();
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
            //write_traj();
            _mc_queue = "samc:20000:10-0";
            mc_run();
        }

        Num en_crash(const Residue &r1, const Residue &r2) {
            return 1000*scorer->en_crash(r1, r2);
        }

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
            return 10*e;

        }

        /*
        Num en_bp(const Residue &r1, const Residue &r2) {
            Num d1 = geom::distance(r1["P"], r2["P"]);
            Num d2 = geom::distance(r1["C1*"], r2["C1*"]);
            if (d1 < 17 || d1 > 20) {
                return JN_ square(d1-18);
            }
            else if (d2 > 12) {
                return JN_ square(d2-12);
            }
            else {
                scorer->en_bp(r1, r2);
                return 50*scorer->m_en_pairing;
            }
        }
        */

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
            Int l = size(ss);
            En e;
            for (Int i = 0; i <= ind; i++) {
                for (Int j = ind+1; j < l; j++) {
                    e.crash += en_crash(pred[i], pred[j]);
                }
            }
            for (auto && bp : bps) {
                //std::cout << bp[0] << ' ' << bp[1] << std::endl;
                e.pairing += en_bp(pred[bp[0]], pred[bp[1]]);
            }
            e.len += scorer->en_len(pred, ind);
            return e.crash + e.pairing + e.len;
        }

        virtual void mc_write() {
            mc_total_energy();
            //write_traj();
            /*
            std::cout
                << "total:" << en.crash+en.pairing+en.len
                << " crash:" << en.crash
                << " pairing:" << en.pairing
                << " len:" << en.len
                << " tempr:" << _mc_tempr
                << std::endl;
            */
        }

        void write_traj() {
            static MolWriter writer;
            static Int n = 0;
            n++;
            if (n == 1) file::clean("aa.traj.pdb");
            std::ofstream output("aa.traj.pdb", std::ios::app);
            writer.bind_stream(output);
            writer.write_model([this](){writer.write(pred);});
            output.close();
        }

        virtual void mc_select() {
            ind = Int(JN_ rand() * (size(ss)-1));
        }

        virtual void mc_sample() {
            mc_backup();
            sample_loop(pred, ss, ind);
        }

        void mc_backup() {
            Int l = size(ss);
            restored.clear();
            if (ind < l/2) for (Int i = 0; i <= ind; i++) restored.push_back(pred[i]);
            else for (Int i = ind+1; i < l; i++) restored.push_back(pred[i]);
        }

        virtual void mc_back() {
            Int l = size(ss);
            if (ind < l/2) for (Int i = 0; i <= ind; i++) pred[i] = restored[i];
            else for (Int i = ind+1; i < l; i++) pred[i] = restored[i-ind-1];
        }
};

Chain build_loop(Str seq, Str ss) {
    SampleLoop sampler;
    sampler.pred = init_loop(seq, ss);
    sampler.ss = ss;
    sampler._mc_init_tempr = 1;
    sampler.run();
    return sampler.pred;
}

END_JN
