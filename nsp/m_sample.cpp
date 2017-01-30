#include "nsp.hpp"
#include <jnbio/nuc3d/Assemble.hpp>
#include <jnbio/pdb/utils/cluster_chains.hpp>
#include <jnbio/scoring/Score.hpp>

BEGIN_JN

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

	if (write_samplings) mol_write(ass._pred_chain, to_str(ass._name, ".", n, ".sample.pdb"));
	for (n = 2; n <= num_samplings; n++) {
		ass.sample_one_template();
		ass.assemble();
		chains.push_back(std::move(ass._pred_chain));
		ass.log << "# Writing sampling structure " << n << std::endl;
		if (write_samplings) mol_write(ass._pred_chain, to_str(ass._name, ".", n, ".sample.pdb"));
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
		mol_write(chains[cluster[ind]], to_str(ass._name, ".", n_cluster + 1, ".pred.pdb"));
		n_cluster++;
	}
    mol_write(chains[min_ind], to_str(ass._name, ".pred.pdb"));

}

END_JN

