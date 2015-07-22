#include "../Score"
#include "../Assemble"

int main(int argc, char **argv) {
	if (!strcmp(argv[1], "-score:list")) {
		string str;
		Obj<Score> score;
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(argv[3]);
		}

		ifstream ifile(argv[2]);
		while (ifile >> str) {
			Obj<RNA> rna = new RNA(str);
			cout << str << ' ' << score->run(rna) << endl;
		}
		ifile.close();
	} else if (!strcmp(argv[1], "-score")) {
		Obj<Score> score;
		Obj<RNA> rna = new RNA(argv[2]);
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(argv[3]);
		}
		cout << score->run(rna) << endl;
		double *scores = score->distAnal->score;
		cerr << rna->name << ' ' << scores[0] << ' ' << scores[1] << ' ' << scores[2] << ' ' << scores[3] << ' ' << scores[4] << ' ' << score->dihAnal->getScore() << endl;
	} else if (!strcmp(argv[1], "-score:more")) {
		Obj<Score> score;
		Obj<RNA> rna = new RNA(argv[2]);
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(argv[3]);
		}
		score->run(rna);
		int len = rna->getLen();
		cout << "Nucleotide internal energy: " << endl;
		for (int i = 0; i < len; i++) {
			cout << score->distAnal->nuc_score_[i] << ' ';
		}
		cout << endl;
		cout << "Base-stacking energy: " << endl;
		for (int i = 0; i < len - 1; i++) {
			cout << score->distAnal->stacking_score_[i] << ' ';
		}
		cout << endl;
		cout << "Base-pairing energy: " << endl;
		score->distAnal->pairing_score_->print();
	} else if (!strcmp(argv[1], "-assemble")) {
		Assemble ass(argv[2]);
		ass.run();
		RNAs *rnas = ass.getResults();
		for (int i = 0; i < rnas->getLen(); i++) {
			stringstream ss;
			string name;
			ss << ass.output_path();
			ss << ass.job_name();
			ss << '-';
			ss << i + 1;
			ss << ".pdb";
			ss >> name;
			(*rnas)[i]->write(name);
		}
		delete rnas;
	} else if (!strcmp(argv[1], "-assemble:test")) {
		Assemble ass(argv[2]);
		ass.is_test = 1;
		ass.run();
		RNAs *rnas = ass.getResults();
		for (int i = 0; i < rnas->getLen(); i++) {
			stringstream ss;
			string name;
			ss << ass.output_path();
			ss << ass.job_name();
			ss << '-';
			ss << i + 1;
			ss << ".pdb";
			ss >> name;
			(*rnas)[i]->write(name);
		}
		delete rnas;
	} else if (!strcmp(argv[1], "-mc")) {
	} else if (!strcmp(argv[1], "-connect")) {
		string str1(argv[2]);
		string str2(argv[3]);
		RNA *rna1 = new RNA(str1);
		RNA *rna2 = new RNA(str2);
		Connect connect;
		RNA *rna = new RNA(connect(*rna1, *rna2, atoi(argv[4]), atoi(argv[5])));
		rna->print();
		delete rna1;
		delete rna2;
		delete rna;
	} else if (!strcmp(argv[1], "-updateChains")) {
		RNA *rna = new RNA(argv[2]);
		ifstream ifile(argv[3]);
		string ss;
		ifile >> ss;
		rna->updateChains(ss);
		rna->print();
		ifile.close();
		delete rna;
	} else if (!strcmp(argv[1], "-mutate")) {
		RNA *rna = new RNA(argv[2]);
		rna->mutate(argv[3]);
		rna->print();
		delete rna;
	} else if (!strcmp(argv[1], "-seq")) {
		Pdb pdb(argv[2]);
		std::string fill = (argc == 4 ? argv[3] : "");
		for (auto &chain: pdb[0].chains) {
			for (auto &residue: chain) {
				std::cout << residue.name << fill;
			}
		}
		std::cout << std::endl;
	} else if (!strcmp(argv[1], "-chain")) {
		Obj<RNA> rna = new RNA(argv[2]);
		cout << rna->getChain() << endl;
	} else if (!strcmp(argv[1], "-len")) {
		RNA *rna = new RNA(argv[2]);
		cout << rna->getLen() << endl;
		delete rna;
	} else if (!strcmp(argv[1], "-len:atoms")) {
		RNA *rna = new RNA(argv[2]);
		cout << rna->totalAtoms() << endl;
		delete rna;
	} else if (!strcmp(argv[1], "-split")) {
		Split split(argv[2]);
		split.run();
	} else if (!strcmp(argv[1], "-rmsd:p")) {
		RMSD rmsd(argv[2]);
		cout << rmsd.run() << endl;
	} else if (!strcmp(argv[1], "-rmsd")) {
		Obj<RNA> rna1(new RNA(argv[2]));
		Obj<RNA> rna2(new RNA(argv[3]));
		RMSD rmsd(rna1, rna2);
		cout << rmsd.run() << endl;
	} else if (!strcmp(argv[1], "-rmsd:xy")) {
		RMSD rmsd(argv[2], argv[3]);
		cout << rmsd.run() << endl;
//	} else if (!strcmp(argv[1], "-lm5")) {
//		string str1(argv[2]);
//		string str2(argv[3]);
//		LoopModelling5 lm5(str1, str2);
//		lm5.run();
	} else if (!strcmp(argv[1], "-loop")) {
		string str1(argv[2]);
		string str2(argv[3]);
		LoopModelling loopModelling(str1, str2);
		RNA *rna = loopModelling.run();
		rna->print();
//	} else if (!strcmp(argv[1], "-lm")) {
//		string str1(argv[2]);
//		string str2(argv[3]);
//		LM lm(str1, str2);
//		lm.run();
	} else if (!strcmp(argv[1], "-mol2d")) {
		Mol2D mol2d(argv[2], 1);
		mol2d.print();
//	} else if (!strcmp(argv[1], "-lm2")) {
//		string str1(argv[2]);
//		string str2(argv[3]);
//		LM2 lm(str1, str2);
//		lm.run();
//	} else if (!strcmp(argv[1], "-loop2")) {
//		string str1(argv[2]);
//		string str2(argv[3]);
//		LoopModelling2 loopModelling(str1, str2);
//		RNA *rna = loopModelling.run();
//		rna->print();
//	} else if (!strcmp(argv[1], "-loop3")) {
//		string str1(argv[2]);
//		string str2(argv[3]);
//		LoopModelling3 loopModelling3(str1, str2);
//		RNA *rna = loopModelling3.run();
//		rna->print();
//	} else if (!strcmp(argv[1], "-loop4")) {
//		string str1(argv[2]);
//		string str2(argv[3]);
//		LoopModelling4 loopModelling4(str1, str2);
//		RNA *rna = loopModelling4.run();
//		rna->print();
	} else if (!strcmp(argv[1], "-train:distance")) {
		Obj<Train> train;
		if (argc == 3) {
			train = new Train;
		} else if (argc == 4) {
			train = new Train(argv[3]);
		}
		train->dist(argv[2]);
	} else if (!strcmp(argv[1], "-train:dihedral")) {
		Train train;
		train.dih(argv[2]);
	} else if (!strcmp(argv[1], "-test")) {
		char *p = NULL;
		cout << (p ? p : "") << endl;
	} else if (!strcmp(argv[1], "-r6p")) {
		R6P r6p(argv[2]);
		cout << r6p.res_nums() << ' ' << 6 << endl;
		for (auto &chain: r6p.chains) {
			for (auto &residue: chain.residues) {
				for (auto &atom: residue.atoms) {
					for (auto &chain2: r6p.chains) {
						for (auto &residue2: chain2.residues) {
							for (auto &atom2: residue2.atoms) {
								cout << Point(atom.x, atom.y, atom.z).dist(
								        Point(atom2.x, atom2.y, atom2.z)) << ' ';
							}
						}
					}
					cout << endl;
				}
			}
		}
	} else if (!strcmp(argv[1], "-r5p")) {
		R5P r5p(argv[2]);
		cout << r5p.res_nums() << ' ' << 6 << endl;
		for (auto &chain: r5p.chains) {
			for (auto &residue: chain.residues) {
				for (auto &atom: residue.atoms) {
					for (auto &chain2: r5p.chains) {
						for (auto &residue2: chain2.residues) {
							for (auto &atom2: residue2.atoms) {
								cout << Point(atom.x, atom.y, atom.z).dist(
								        Point(atom2.x, atom2.y, atom2.z)) << ' ';
							}
						}
					}
					cout << endl;
				}
			}
		}
	} else if (!strcmp(argv[1], "-r5p:mono_nuc")) {
		R5P r5p(argv[2]);
		for (auto &chain: r5p.chains) {
			for (auto &residue: chain.residues) {
				if (residue.name != string(argv[3])) continue;
				for (auto &atom1: residue.atoms) {
					if (atom1.name != string(argv[4])) continue;
					for (auto &atom2: residue.atoms) {
						if (atom2.name != string(argv[5])) continue;
						double x = atom1.x - atom2.x;
						double y = atom1.y - atom2.y;
						double z = atom1.z - atom2.z;
						double dist = sqrt(x * x + y * y + z * z);
						cout << dist << endl;
					}
				}
			}
		}
	} else if (!strcmp(argv[1], "-r5p:adj_nuc")) {
		R5P r5p(argv[2]);
		int res_num = 0;
		Residue old_res;
		for (auto &chain: r5p.chains) {
			for (auto &residue: chain.residues) {
				if (res_num != 0) {
					if (old_res.name == string(argv[3]) && residue.name == string(argv[4])) {
						for (auto &atom1: old_res.atoms) {
							if (atom1.name != string(argv[5])) continue;
							for (auto &atom2: residue.atoms) {
								if (atom2.name != string(argv[6])) continue;
								double x = atom1.x - atom2.x;
								double y = atom1.y - atom2.y;
								double z = atom1.z - atom2.z;
								double dist = sqrt(x * x + y * y + z * z);
								cout << dist << endl;
							}
						}
					}
				}
				res_num++;
				old_res = residue;
			}
		}
	} else if (!strcmp(argv[1], "-aa")) {
		RNA rna(argv[2]);
		for (auto &chain: rna.chains) {
			for (auto &residue: chain.residues) {
				if (residue.name == string(argv[3])) {
					for (auto &atom1: residue.atoms) {
						if (atom1.name != string(argv[4])) continue;
						for (auto &atom2: residue.atoms) {
							if (atom2.name != string(argv[5])) continue;
							double x = atom1.x - atom2.x;
							double y = atom1.y - atom2.y;
							double z = atom1.z - atom2.z;
							double dist = sqrt(x * x + y * y + z * z);
							cout << dist << endl;
						}
					}
				}
			}
		}
	}
	return 0;
}



