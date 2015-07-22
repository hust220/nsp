#include "Model.h"

namespace jian {

Model::Model(std::string pdbfile) {
    /* set name */
    assert(pdbfile.substr(pdbfile.size() - 4, 4) == ".pdb");
    name = pdbfile.substr(0, pdbfile.size() - 4);

    /* set chains */
    ifstream ifile(pdbfile.c_str());
    assert(ifile.is_open());

    string line;
    vector<string> lines;
    int n = 0;
    while (getline(ifile, line, '\n')) {
        if (!line.compare(0, 4, "ATOM")) {
            n++;
            if (lines.size() != 0 && line[21] != lines.back()[21]) {
                chains.push_back(Chain(lines, name));
                lines.clear();
            }
            lines.push_back(line);
        }
    }

    n != 0 || die("The file '" + pdbfile + "' has nothing!");

    chains.push_back(Chain(lines, name));
    ifile.close();
}

Model::Model(vector<string> lines) {
	name = "none";
	vector<string> chain_lines;
	char chain_name = ' ';
	int atom_num = 0;

	for (auto &line: lines) {
		atom_num++;
		if (line[21] != chain_name && atom_num != 1) {
			push(Chain(chain_lines));
			chain_lines.clear();
		}
		chain_name = line[21];
		chain_lines.push_back(line);
	}
	push(Chain(chain_lines));
	chain_lines.clear();
}

Chain &Model::operator [](int n) {
	return chains[n];
}

const Chain &Model::operator [](int n) const {
	return chains[n];
}

void Model::push(const Chain &chain) {
	chains.push_back(chain);
}

ostream &operator <<(ostream &output, const Model &model) {
	int atom_num = 1;
	int residue_num = 1;
	output << fixed << setprecision(3);
	for (auto &chain: model.chains) {
		for (auto &residue: chain.residues) {
			for (auto &atom: residue.atoms) {
				output << "ATOM" 
				       << setw(7)  << atom_num << "  "
					   << left << setw(4)  << atom.name
					   << right << setw(3) << residue.name
					   << setw(2)  << chain.name 
					   << setw(4)  << residue_num 
					   << setw(12) << atom.x 
					   << setw(8)  << atom.y 
					   << setw(8)  << atom.z 
					   << "\n";
				atom_num++;
			}
			residue_num++;
		}
		output << "TER" << endl;
	}
	return output;
}

void Model::print() {
	cout << *this;
}

void Model::write(string pdbname) {
	ofstream ofile(pdbname.c_str());
	ofile << *this;
	ofile.close();
}

vector<Chain>::iterator Model::begin() {
    return chains.begin();
}
vector<Chain>::iterator Model::end() {
    return chains.end();
}

} /// namespace jian



