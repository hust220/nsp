#include "nsp.hpp"
#include <jian/nuc3d/Assemble.hpp>
#include <jian/pdb/utils/cluster_chains.hpp>

BEGIN_JN

REGISTER_NSP_COMPONENT(assemble) {
	nuc3d::Assemble ass(par);

	std::ostringstream stream;
	int n;

	int num = 1;
	par.set(num, "n", "num");

	ass.select_templates();

	n = 1;
	mol_write(ass._pred_chain, to_str(ass._name, ".pred.", n, ".pdb"));
	for (n = 2; n <= num; n++) {
		ass.sample_all_templates();
		ass.assemble();
		ass.log << "# Writing sampling structure " << n << std::endl;
		mol_write(ass._pred_chain, to_str(ass._name, ".pred.", n, ".pdb"));
	}

}

REGISTER_NSP_COMPONENT(sample) {
	nuc3d::Assemble ass(par);

	std::ostringstream stream;
	int n;
	bool write_samplings = par.has("write_samplings");

	int num = 1;
	par.set(num, "n", "num");

	int num_samplings = 100 * num;
	par.set(num_samplings, "num_samplings");

	ass.select_templates();

	n = 1;
	std::deque<Chain> chains;
	chains.push_back(std::move(ass._pred_chain));

	if (write_samplings) mol_write(ass._pred_chain, to_str(ass._name, ".sample.", n, ".pdb"));
	for (n = 2; n <= num_samplings; n++) {
		ass.sample_one_template();
		ass.assemble();
		chains.push_back(std::move(ass._pred_chain));
		ass.log << "# Writing sampling structure " << n << std::endl;
		if (write_samplings) mol_write(ass._pred_chain, to_str(ass._name, ".sample.", n, ".pdb"));
	}

	ass.log << "# Clustering..." << std::endl;
	auto result = pdb::cluster_chains(chains, n);

	for (int i = 0; i < n; i++) {
		mol_write(chains[result[i][0]], to_str(ass._name, ".pred.", i + 1, ".pdb"));
	}

}

END_JN

