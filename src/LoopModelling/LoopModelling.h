/**
 * @file LoopModelling.cpp
 * @brief Modelling loop by using distance geometry algorithm.
 *
 * Each nucleotide is represented by 6 backbone atoms.
 *
 * @author wj_hust08@163.com
 * @version 1.0
 * @date 2015-5-25
*/

#ifndef LOOPMODELLING_H
#define LOOPMODELLING_H

#include "../Pdb.h"
#include "../Assemble_/Mol2D.h"

class LoopModelling {
public:
	LoopModelling(string seq, string ss) {
		this->seq = seq;
		this->ss = ss;
		int len1 = seq.size();
		int len2 = count_if(ss.begin(), ss.end(), [](char c) {
			return set<char>{'.', '(', ')', '[', ']'}.count(c);
		});
		assert(len1 == len2);
		resLen = len1;
		len = len1 * 6;
		type = new int[resLen];
		map<char, int> type_map{{'A', 0}, {'U', 1}, {'G', 2}, {'C', 3}};
		for (int i = 0; i < resLen; i++) {
			type[i] = type_map[seq[i]];
		}
		chir = NULL;

		string lib = getenv("RNA");

		string filename = lib + "distance.par";
		ifstream ifile(filename.c_str());
		assert(ifile);
		par = new double[144];
		string str;
		for (int i = 0; i < 72; i++) {
			ifile >> str >> par[2 * i] >> par[2 * i + 1];
		}
		ifile.close();

		// read parameters
		string helix_par_file = getenv("RNA");
		helix_par_file += "/r6p_helix.par";
		ifile.open(helix_par_file.c_str());
		int res_nums, num_atom_per_residue;
		ifile >> res_nums >> num_atom_per_residue;
		int total_atom_nums = res_nums * num_atom_per_residue;
		helix_par.resize(total_atom_nums, total_atom_nums);
		for (int i = 0; i < total_atom_nums; i++) {
			for (int j = 0; j < total_atom_nums; j++) {
				ifile >> helix_par(i, j);
			}
		}
		ifile.close();

		filename = lib + "chirality";
		ifile.open(filename.c_str());
		assert(ifile);
		for (int i = 0; i < 6; i++) {
			ifile >> chirality[i];
		}
		ifile.close();

		init();
		dg = DG(bound);
	}

	void init();
	RNA *run();
	void set_base_pairs(loop *);
	~LoopModelling();

	double at(MatrixXf &, int, int, int);
	void assign(MatrixXf &, int, int, double, int);
	Residue *buildNuc(Point *, Point *, double, int);

	string seq;
	string ss;
	int len;
	int *type;
	int resLen;
	MatrixXf bound;
	double chirality[6];
	Matr_ *chir;
	MatrixXf backbone;
	double *par;
	MatrixXf helix_par;
	DG dg;
	int atom_nums_per_nuc = 6;
	double err_radius = 0.001;
	string lib;
};

#endif

