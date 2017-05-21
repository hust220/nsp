#include "nsp.hpp"
#include "geom.hpp"
#include "pdb.hpp"
#include "score_aa.hpp"
#include "score.hpp"
#include "file.hpp"
#include "rtsp_format.hpp"
#include "log.hpp"

BEGIN_JN

namespace {

    Format formater;

    void read_chain(Chain &chain, S filename) {
        chain_read_model(chain, filename);
        chain = formater(chain);
    }

    void sum_counts(S filename, int rows, int cols) {
        Veci v;
        int i, j, d;
        std::ifstream ifile;

        v = Veci::Zero(cols);
        FOPEN(ifile, filename);
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                ifile >> d;
                v[j] += d;
            }
        }
        FCLOSE(ifile);

        for (i = 0; i < cols; i++) {
            JN_OUT << v[i] << ' ';
        }
        JN_OUT << std::endl;
    }

    void score_res(Score * scoring, S filename, S score_type = "pairing") {
        Chain chain;
        int i, j, l;

        chain_read_model(chain, filename);
        l = chain.size();
        //JN_JN_OUT << l << std::endl;
        for (i = 0; i < l; i++) {
            for (j = 0; j < l; j++) {
                if (i == j) JN_OUT << 0 << "\t";
                else {
                    scoring->en_bp(chain[i], chain[j]);
                    if (score_type == "pairing") JN_OUT << scoring->m_en_pairing << "\t";
                    else if (score_type == "stacking") JN_OUT << scoring->m_en_stacking << "\t";
                    else if (score_type == "wc") JN_OUT << scoring->m_en_wc << "\t";
                    else if (score_type == "nwc") JN_OUT << scoring->m_en_nwc << "\t";
                    else {
                        LOG << "Illegal score type: " << score_type << std::endl;
                        throw "error!";
                    }
                }
            }
            JN_OUT << std::endl;
        }
    }

    void score_s(Score * scoring, S filename) {
        Chain chain;
        chain_read_model(chain, filename);
        scoring->run(chain);
        JN_OUT <<
            "Score of " << filename << ": " <<
            //scoring->m_score_dih << "(dih) " <<
            //scoring->m_score_dist << "(dist) " <<
            scoring->m_score << "(total)" <<
            std::endl;
    }

    void score_l(Score * scoring, S filename) {
        for (auto &&it : FileLines(filename)) {
            score_s(scoring, it.arr[0]);
        }
    }

    void train_s(Score * scoring, S filename) {
        Chain chain;

        LOG << "Train " << filename << " ..." << std::endl;
        chain_read_model(chain, filename);
        scoring->train(chain);
    }

    void train_l(Score * scoring, S filename) {
        for (auto &&it : FileLines(filename)) {
            train_s(scoring, it.arr[0]);
        }
    }

    double en_crash(const Residue &r1, const Residue &r2) {
        int i, j;
        double d, e;

        e = 0;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                d = geom::distance(r1[i], r2[j]);
                if (i == 0 || j == 0) {
                    if (d < 5) {
                        e += square(d - 4);
                    }
                }
                else if ((i == 1 || j == 1) && d < 5) {
                    e += square(d - 5);
                }
                else if (d < 3.5) {
                    e += square(d - 3.5);
                }
            }
        }
        return e;
    }

    REGISTER_NSP_COMPONENT(score) {
        //std::ofstream stream;
        S method;
        S score_type = "pairing";

        par.set(score_type, "score_type");

        CG *m_cg;

        method = "aa";
        par.set(method, "cg");

        m_cg = CG::fac_t::create(method);

        if (par.has("crash")) {
            double e, d;
            Chain chain;
            int i, j, l;

            chain_read_model(chain, par.get("s"));
            chain = m_cg->to_cg(chain);
            l = chain.size();
            e = 0;
            for (i = 0; i < l; i++) {
                for (j = 0; j < l; j++) {
                    d = (i == j ? 0 : en_crash(chain[i], chain[j]));
                    e += d;
                    JN_OUT << d << "\t";
                }
                JN_OUT << std::endl;
            }
            JN_OUT << e << std::endl;

        }
        else if (par.has("sum_counts")) {
            Par::pars_t & pars = par["sum_counts"];
            S filename = pars[0];
            int rows = std::stoi(pars[1]);
            int cols = std::stoi(pars[2]);
            sum_counts(filename, rows, cols);
        }
        else {
            Score *scoring = Score::fac_t::create(method);
            scoring->init();

            if (par.has("print_freqs")) {
                scoring->print_freqs(JN_OUT);
            }
            else if (par.has("print_counts")) {
                scoring->print_counts(JN_OUT);
            }
            else if (par.has("train")) {
                if (par.has("s")) {
                    train_s(scoring, par.get("s"));
                }
                else if (par.has("l")) {
                    train_l(scoring, par.get("l", "list"));
                }
                scoring->print_counts(JN_OUT);
            }
            else {
                if (par.has("s")) {
                    if (par.has("res")) {
                        score_res(scoring, par.get("s"), score_type);
                    }
                    else {
                        score_s(scoring, par.get("s"));
                    }
                }
                else if (par.has("l")) {
                    score_l(scoring, par.get("l", "list"));
                }
            }
            delete scoring;
        }
        delete m_cg;
    }
}

END_JN
















