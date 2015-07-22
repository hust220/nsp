#include "Connect.h"

using namespace jian;
using namespace jian::nuc3d;

Model Connect::operator ()(const Model &a, const Model &b, int a1, int a2) {
//	if (_a) delete _a;
//	if (_b) delete _b;
	_model_1 = a;
	_model_2 = b;

	int len_1 = _model_1.res_nums();
	int len_2 = _model_2.res_nums();

	MatrixXf x(12, 3), y(12, 3);
	int num_res = 0;
	int temp = 0;
	for (auto &chain : _model_1.chains) {
		for (auto &residue : chain.residues) {
			if (set<int>{a1 - 1, a1, a2, a2 + 1}.count(num_res)) {
				for (auto &atom : residue.atoms) {
					if (_superposed_atoms.count(atom.name)) {
						x(temp, 0) = atom.x;
						x(temp, 1) = atom.y;
						x(temp, 2) = atom.z;
						temp++;
					}
				}
			}
			num_res++;
		}
	}
	num_res = 0;
	temp = 0;
	for (auto &chain : _model_2.chains) {
		for (auto &residue : chain.residues) {
			if (set<int>{0, 1, len_2 - 2, len_2 - 1}.count(num_res)) {
				for (auto &atom : residue.atoms) {
					if (_superposed_atoms.count(atom.name)) {
						y(temp, 0) = atom.x;
						y(temp, 1) = atom.y;
						y(temp, 2) = atom.z;
						temp++;
					}
				}
			}
			num_res++;
		}
	}

	SupPos sp;
	sp(x, y);

	translate(_model_1, sp.c1);
	translate(_model_2, sp.c2);
	rotate(_model_1, sp.rot);

	// construct a new RNA
	RNA temp_rna;
	Chain temp_chain;
	num_res = 0;
	for (auto &chain : _model_1.chains) {
		for (auto &residue : chain.residues) {
			if (num_res < a1 - 1 || num_res > a2 + 1) {
				temp_chain.residues.push_back(residue);
			} else if (num_res == a1) {
				for (auto &chain2 : _model_2.chains) {
					for (auto &residue2 : chain2.residues) {
						temp_chain.residues.push_back(residue2);
					}
				}
			}
			num_res++;
		}
	}
	temp_rna.chains.push_back(temp_chain);
	return temp_rna;
}

void Connect::translate(Model &model, const Point &point) {
	for (auto &chain : model.chains) {
		for (auto &residue : chain.residues) {
			for (auto &atom : residue.atoms) {
				atom.x -= point.x;
				atom.y -= point.y;
				atom.z -= point.z;
			}
		}
	}
}

void Connect::rotate(Model &model, const Matrix3f &matrix) {
	for (auto &chain : model.chains) {
		for (auto &residue : chain.residues) {
			for (auto &atom : residue.atoms) {
				RowVector3f vec = RowVector3f(atom.x, atom.y, atom.z) * matrix;
				atom.x = vec(0);
				atom.y = vec(1);
				atom.z = vec(2);
			}
		}
	}
}

