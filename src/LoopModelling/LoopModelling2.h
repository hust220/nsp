/**
 * @file LoopModelling22.cpp
 * @brief Modelling loop by using distance geometry algorithm.
 *
 * Each nucleotide is represented by 5 atoms.
 *   A: C5* O3* C1* N6 C2
 *   U: C5* O3* C1* O2 O4
 *   G: C5* O3* C1* O6 N2
 *   C: C5* O3* C1* O2 N4
 *
 * @author wj_hust08@163.com
 * @version 1.0
 * @date 2015-5-25
*/

#ifndef LOOPMODELLING2_H
#define LOOPMODELLING2_H

#include "../Pdb.h"
#include "../Assemble_/Mol2D.h"
#include "../Assemble_/Connect.h"

class LoopModelling2 {
public:
	LoopModelling2() {
		char *lib = getenv("RNA3D");
		assert(lib);
		_lib = string() + lib;

		// read mononucleotide parameters
		string file_name = _lib + "/pars/r5p/mono_nuc_pars";
		ifstream ifile(file_name.c_str());
		assert(ifile);
		for (int ii = 0; ii < 4; ii++) {
			string c;
			ifile >> c;
			_mono_nuc_pars[c].resize(_atom_nums_per_nuc, _atom_nums_per_nuc);
			for (int i = 0; i < _atom_nums_per_nuc; i++) {
				for (int j = 0; j < _atom_nums_per_nuc; j++) {
					ifile >> _mono_nuc_pars[c](i, j);
				}
			}
		}
		ifile.close();

		// read adjacent nucleotides parameters
		file_name = _lib + "/pars/r5p/adj_nuc_pars";
		ifile.open(file_name.c_str());
		for (int ii = 0; ii < 16; ii++) {
			string str;
			ifile >> str;
			_adj_nuc_pars[str].resize(_atom_nums_per_nuc, _atom_nums_per_nuc);
			for (int i = 0; i < _atom_nums_per_nuc; i++) {
				for (int j = 0; j < _atom_nums_per_nuc; j++) {
					ifile >> _adj_nuc_pars[str](i, j);
				}
			}
		}
		ifile.close();

		// read all-atoms single nucleotide parameters
		file_name = _lib + "/pars/r5p/aa.pars";
		ifile.open(file_name.c_str());
		map<string, int> temp_map{{"A", 22}, {"U", 20}, {"G", 23}, {"C", 20}};
		for (int ii = 0; ii < 4; ii++) {
			string str;
			ifile >> str;
			int atom_nums = temp_map[str];
			_aa_pars[str].resize(atom_nums, atom_nums);
			for (int i = 0; i < atom_nums; i++) {
				for (int j = 0; j < atom_nums; j++) {
					ifile >> _aa_pars[str](i, j);
				}
			}
		}
		ifile.close();
	}
	RNA operator ()(string seq, string ss, int num = 1);
	RNA to_all_atom(const MatrixXf &);

	map<string, MatrixXf> _mono_nuc_pars;
	map<string, MatrixXf> _adj_nuc_pars;
	map<string, MatrixXf> _aa_pars;

	void init();
	void set_base_pairs(loop *);

//	double at(MatrixXf &, int, int, int);
//	void assign(MatrixXf &, int, int, double, int);
	MatrixXf get_helix_par(const R5P &);
	R5P get_helix(const helix &);
	R5P create_helix(const string &);
	Residue build_nuc(const string &, const MatrixXf &);
	Residue make_residue(const string &, const MatrixXf &);
//	MatrixXf get_par_from_file(const string &);

	string _seq;
	string _ss;
	MatrixXf _bound;
	DG dg;
	int _atom_nums_per_nuc = 5;
	double _err_radius = 0.001;
	string _lib;

	int _view = 0;
};

#endif

