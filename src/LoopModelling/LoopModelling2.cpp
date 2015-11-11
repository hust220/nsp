#include "LoopModelling2.h"

RNA LoopModelling2::operator ()(string seq, string ss, int num) {
	_seq = seq;
	_ss = ss;
	assert(seq.size() == count_if(ss.begin(), ss.end(), [](char c) {
		return set<char>{'.', '(', ')', '[', ']'}.count(c);
	}));

	init();
	dg = DG(_bound);
	dg.view = _view;
	return to_all_atom(dg());
}

RNA LoopModelling2::to_all_atom(const MatrixXf &scaffold) {
	int res_nums = _seq.size();
	Chain chain;

	for (int i = 0; i < res_nums; i++) {
		chain.residues.push_back(
			build_nuc(
				string() + _seq[i],
				scaffold.block(i * _atom_nums_per_nuc, 0, _atom_nums_per_nuc, 3)
			)
		);
	}

	RNA rna;
	rna.chains.push_back(chain);
	return rna;
}

void LoopModelling2::init() {
	int res_nums = _seq.size();
	int atom_nums = res_nums * _atom_nums_per_nuc;

	/// begin initializing bound matrix
	_bound.resize(atom_nums, atom_nums);
	for (int i = 0; i < atom_nums; i++) {
		for (int j = i; j < atom_nums; j++) {
			if (i == j) {
				_bound(i, j) = 0;
			} else {
				_bound(j, i) = 3;
				_bound(i, j) = 999;
			}
		}
	}
	/// end initializing bound matrix

	/// begin mononucleotide
	for (int i = 0; i < res_nums; i++) {
		for (int j = 0; j < _atom_nums_per_nuc; j++) {
			for (int k = j + 1; k < _atom_nums_per_nuc; k++) {
				_bound(i * _atom_nums_per_nuc + j, i * _atom_nums_per_nuc + k)
					= _mono_nuc_pars[string() + _seq[i]](j, k);
				_bound(i * _atom_nums_per_nuc + k, i * _atom_nums_per_nuc + j)
					= _mono_nuc_pars[string() + _seq[i]](j, k);
			}
		}
	}
	/// end mononucleotide

	/// begin adjacent nucleotides
	for (int i = 0; i + 1 < res_nums; i++) {
		int j = i + 1;
		for (int m = 0; m < _atom_nums_per_nuc; m++) {
			for (int n = 0; n < _atom_nums_per_nuc; n++) {
				string str = string() + _seq[i] + _seq[j];
				_bound(i * _atom_nums_per_nuc + m, j * _atom_nums_per_nuc + n)
					= _adj_nuc_pars[str](m, n);
				_bound(j * _atom_nums_per_nuc + n, i * _atom_nums_per_nuc + m)
					= _adj_nuc_pars[str](m, n);
			}
		}
	}
	/// end adjacent nucleotides

	/// begin base pairs
	Mol2D mol2d(_ss, _seq);
	set_base_pairs(mol2d.head);

	string str = _ss;
	transform(str.begin(), str.end(), str.begin(), [](const char &c) {
		map<char, char> temp_map{{'(', '['}, {')', ']'}, {'[', '('}, 
		                         {']', ')'}, {'.', '.'}, {'&', '&'}};
		return temp_map[c];
	});
	Mol2D mol2d_2(str, _seq);
	set_base_pairs(mol2d_2.head);
	/// end base pairs

	/// chirality
//	int flag = 0;
//	for (int i = 3; i < len;) {
//		if (ss[i / 6] == '(' && i % 6 == 5 && i != 5) {
//			i += 4;
//		} else {
//			i++;
//		}
//		flag++;
//	}
//	delete chir;
//	chir = new Matr_(flag, 5);
//	for (int i = 3, flag = 0; i < len;) {
//		chir->data[flag][0] = i - 3;
//		chir->data[flag][1] = i - 2;
//		chir->data[flag][2] = i - 1;
//		chir->data[flag][3] = i;
//		chir->data[flag][4] = chirality[(i - 2) % 6];
//		if (ss[i / 6] == '(' && i % 6 == 5 && i != 5) {
//			i += 4;
//		} else {
//			i++;
//		}
//		flag++;
//	}

}

void LoopModelling2::set_base_pairs(loop *src) {
	if (src == NULL) {
		return;
	} else {
		if (src->s.head != NULL) {
			MatrixXf helix_par = get_helix_par(get_helix(src->s));
			int helix_len = src->s.getLen();
			int orig_pos_chain1 = (src->s.head->res1.num - 1) * _atom_nums_per_nuc;
			int orig_pos_chain2 = (src->s.head->res2.num - helix_len) * _atom_nums_per_nuc;
			int atom_nums_per_chain = _atom_nums_per_nuc * helix_len;
			for (int i: {0, 1}) {
				for (int j: {0, 1}) {
					for (int ii = 0; ii < atom_nums_per_chain; ii++) {
						for (int jj = 0; jj < atom_nums_per_chain; jj++) {
							if (i == j && ii == jj) {
								continue;
							}
							int a = orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * i + ii;
							int b = orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * j + jj;
							int c = i * (helix_par.rows() - atom_nums_per_chain) + ii;
							int d = j * (helix_par.rows() - atom_nums_per_chain) + jj;
							double dist = helix_par(c, d);
							if (a < b) {
								_bound(a, b) = dist + _err_radius;
							} else {
								_bound(a, b) = dist - _err_radius;
							}
						}
					}
				}
			}
		}
		set_base_pairs(src->son);
		set_base_pairs(src->brother);
	}
}

R5P LoopModelling2::get_helix(const helix &h) {
	string file_name = _lib + "/info/helix";
	ifstream ifile(file_name.c_str());
	assert(ifile);
	string helix_name, helix_seq, helix_ss, helix_src;
	int helix_len;
	string target_seq = h.seq();
	string pdb_name;
	while (ifile) {
		ifile >> helix_name >> helix_len >> helix_seq >> helix_ss >> helix_src;
		if (target_seq == helix_seq) {
			pdb_name = _lib + "/helix/" + helix_name + ".pdb";
			break;
		}
	}
	ifile.close();
	if (pdb_name == "") {
		return create_helix(target_seq);
	} else {
		return R5P(pdb_name);
	}
}

R5P LoopModelling2::create_helix(const string &seq) {
	if (seq.size() < 4) {
		std::cerr << "Assemble::createHelix error! The length of the helix"
		          << " to be create should not be less than 4!" << std::endl;
		exit(1);
	} else if (seq.size() % 2 == 1) {
		std::cerr << "Assemble::createHelix error! The length of the helix"
		          << " to be create should be even!" << std::endl;
		exit(1);
	} else if (seq.size() == 4 || seq.size() == 6) {
		string file_name = _lib + "/basepair/" + seq + ".pdb";
		ifstream ifile(file_name.c_str());
		if (!ifile) {
			file_name = _lib + "/basepair/XXXXXX.pdb";
		}
		ifile.close();
		return RNA(file_name);
	} else {
		string file_name = _lib + "/basepair/" + seq.substr(0, 3) + 
		                   seq.substr(seq.size() - 3, 3) + ".pdb";
		ifstream ifile(file_name.c_str());
		if (!ifile) {
			file_name = _lib + "/basepair/XXXXXX.pdb";
		}
		ifile.close();
		return Connect()(RNA(file_name), create_helix(seq.substr(1, seq.size() - 2)), 2, 3);
	}
}

MatrixXf LoopModelling2::get_helix_par(const R5P &r5p) {
	int atom_nums = r5p.atom_nums();
	MatrixXf helix_par(atom_nums, atom_nums);
	int num_1 = 0;
	for (auto &chain: r5p.chains) {
		for (auto &residue: chain.residues) {
			for (auto &atom: residue.atoms) {
				int num_2 = 0;
				for (auto &chain2: r5p.chains) {
					for (auto &residue2: chain2.residues) {
						for (auto &atom2: residue2.atoms) {
							helix_par(num_1, num_2) = Point(atom.x, atom.y, atom.z).dist(
									Point(atom2.x, atom2.y, atom2.z));
							num_2++;
						}
					}
				}
				num_1++;
			}
		}
	}
	return helix_par;
}

Residue LoopModelling2::build_nuc(const string &name, const MatrixXf &scaffold) {
	auto par = _aa_pars[name];
	auto vec = map<string, vector<int>>{{"A", {4, 8, 11, 17, 19}},
	                                    {"U", {4, 8, 11, 14, 17}},
										{"G", {4, 8, 11, 17, 20}},
										{"C", {4, 8, 11, 14, 17}}}[name];
	for (int i = 0; i < scaffold.rows(); i++) {
		for (int j = i + 1; j < scaffold.rows(); j++) {
			double x = scaffold(i, 0) - scaffold(j, 0);
			double y = scaffold(i, 1) - scaffold(j, 1);
			double z = scaffold(i, 2) - scaffold(j, 2);
			double dist = sqrt(x * x + y * y + z * z);
			par(vec[i], vec[j]) = dist + _err_radius;
			par(vec[j], vec[i]) = dist - _err_radius;
		}
	}

	DG dg(par);
	dg.view = _view;
	auto coords = dg();

	MatrixXf temp_coords(scaffold.rows(), 3);
	for (int i = 0; i < scaffold.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			temp_coords(i, j) = coords(vec[i], j);
		}
	}
	
	SupPos()(coords, temp_coords, scaffold);

	return make_residue(name, coords);
}

Residue LoopModelling2::make_residue(const string &name, const MatrixXf &coords) {
	Residue residue;
	residue.name = name;
	auto vec = map<string, vector<string>>{
		{"A", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"}},
		{"U", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"}},
		{"G", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"}},
		{"C", {"P", "O1P", "O2P", "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"}}
	}[name];
	for (int i = 0; i < coords.rows(); i++) {
		residue.atoms.push_back(Atom(vec[i], coords(i, 0), coords(i, 1), coords(i, 2)));
	}
	return residue;
}




