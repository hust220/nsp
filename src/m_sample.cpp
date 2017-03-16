#include "nsp.hpp"
#include <nsp/nuc3d/Assemble.hpp>
#include <nsp/pdb/utils/cluster_chains.hpp>
#include <nsp/scoring/Score.hpp>

BEGIN_JN

namespace {

    void write_pred(const Chain &chain, Str dir, Str name, int n) {
        mol_write(chain, to_str(dir, '/', name, ".pred", n, ".pdb"));
    }

    void write_sample(const nuc3d::Assemble &ass, int n) {
        mol_write(ass._pred_chain, to_str(ass.m_out_dir, '/', ass._name, ".sample", n, ".pdb"));
    }

    REGISTER_NSP_COMPONENT(sample) {
        nuc3d::Assemble ass(par);

        std::ostringstream stream;
        int n;
        bool write_samplings = par.has("write_samplings");
        //bool write_samplings = true;

        int num = 1;
        par.set(num, "n", "num");

        int num_samplings = 100 * num;
        par.set(num_samplings, "num_samplings");

        ass.select_templates();
        ass.assemble();

        n = 1;
        std::deque<Chain> chains;
        chains.push_back(std::move(ass._pred_chain));

        if (write_samplings) write_sample(ass, n);
        //if (write_samplings) mol_write(ass._pred_chain, to_str(ass._name, ".", n, ".sample.pdb"));
        for (n = 2; n <= num_samplings; n++) {
            ass.sample_one_template();
            ass.assemble();
            chains.push_back(std::move(ass._pred_chain));
            ass.log << "# Writing sampling structure " << n << std::endl;
            if (write_samplings) write_sample(ass, n);
            //if (write_samplings) mol_write(ass._pred_chain, to_str(ass._name, ".", n, ".sample.pdb"));
        }

        ass.log << "# Clustering..." << std::endl;
        auto clusters = pdb::cluster_chains(chains, num);

        ass.log << "# Scoring..." << std::endl;
        auto scorer = Score::fac_t::make("aa");
        scorer->init();
        Double min_score = 9999;
        Int min_ind = -1;
        Int n_cluster = 0;
        for (auto && cluster : clusters) {
            Deque<Double> scores;

            ass.log << "# Scores of cluster " << n_cluster + 1 << ":";
            for (auto && i : cluster) {
                scorer->run(chains[i]);
                Double d = scorer->m_score;
                ass.log << " " << i << "(" << d << ")";
                scores.push_back(d);
            }
            ass.log << std::endl;

            auto it = std::min_element(scores.begin(), scores.end());
            Int ind = std::distance(scores.begin(), it);
            ass.log << "Minimum score: " << cluster[ind] << "(" << scores[ind] << ")" << std::endl;
            if (scores[ind] < min_score) {
                min_score = scores[ind];
                min_ind = cluster[ind];
            }

            ass.log << "# Writing prediction " << n_cluster + 1 << std::endl;
            write_pred(chains[cluster[ind]], ass.m_out_dir, ass._name, n_cluster+1);
            n_cluster++;
        }
        mol_write(chains[min_ind], to_str(ass.m_out_dir, '/', ass._name, ".pred.pdb"));

    }

}

END_JN

