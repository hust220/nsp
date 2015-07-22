#include "RNA.h"

using namespace jian;

void RNA::move(double dx, double dy, double dz) {
	int i, j, k;
	
	for (i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				chains[i].residues[j].atoms[k].x += dx;
				chains[i].residues[j].atoms[k].y += dy;
				chains[i].residues[j].atoms[k].z += dz;
			}
		}
	}
}

void RNA::rotate(Matr_ *m) {
	int i, j, k;
	double x, y, z;

	for (i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				x = chains[i].residues[j].atoms[k].x;
				y = chains[i].residues[j].atoms[k].y;
				z = chains[i].residues[j].atoms[k].z;
				chains[i].residues[j].atoms[k].x = x * m->data[0][0] + y * m->data[0][1] + z * m->data[0][2];
				chains[i].residues[j].atoms[k].y = x * m->data[1][0] + y * m->data[1][1] + z * m->data[1][2];
				chains[i].residues[j].atoms[k].z = x * m->data[2][0] + y * m->data[2][1] + z * m->data[2][2];
			}
		}
	}
}

void RNA::format() {
	int i, j, k, temp;
	string name;
	Residue *residue;

	for (i = 0; i < (int) chains.size(); i++) {
		chains[i].name = 'A' + i;
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			residue = new Residue;
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				name = chains[i].residues[j].name;
			}
			name = chains[i].residues[j].name;
			for (temp = 0; temp < (int) name.size(); j++) {
				if (name[temp] == 'A' || name[temp] == 'a') {
					chains[i].residues[j].name = 'A';
					break;
				} else if (name[temp] == 'U' || name[temp] == 'u') {
					chains[i].residues[j].name = 'U';
					break;
				} else if (name[temp] == 'G' || name[temp] == 'g') {
					chains[i].residues[j].name = 'G';
					break;
				} else if (name[temp] == 'C' || name[temp] == 'c') {
					chains[i].residues[j].name = 'C';
					break;
				}
			}
		}
	}
}

void RNA::addP() {
	int i, j, k, temp;
	vector<Atom> va;
	Point *p = NULL, *q = NULL;
	Matr_ *a, *b, *c;
	double x1, x2, x, y1, y2, y, z1, z2, z;
	string name;
	Atom *atom;
	double r, r1, r2;

	for (i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			p = new Point[5];
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				if (chains[i].residues[j].atoms[k].name == "P") {
		      p[0].x = chains[i].residues[j].atoms[k].x;
					p[0].y = chains[i].residues[j].atoms[k].y;
			    p[0].z = chains[i].residues[j].atoms[k].z;
			  } else if (chains[i].residues[j].atoms[k].name == "O3*") {
		     	p[1].x = chains[i].residues[j].atoms[k].x;
				 	p[1].y = chains[i].residues[j].atoms[k].y;
					p[1].z = chains[i].residues[j].atoms[k].z;
			  } else if (chains[i].residues[j].atoms[k].name == "O5*") {
			    p[2].x = chains[i].residues[j].atoms[k].x;
			    p[2].y = chains[i].residues[j].atoms[k].y;
			    p[2].z = chains[i].residues[j].atoms[k].z;
				} else if (chains[i].residues[j].atoms[k].name == "C3*") {
			    p[3].x = chains[i].residues[j].atoms[k].x;
		      p[3].y = chains[i].residues[j].atoms[k].y;
					p[3].z = chains[i].residues[j].atoms[k].z;
			  } else if (chains[i].residues[j].atoms[k].name == "C5*") {
			    p[4].x = chains[i].residues[j].atoms[k].x;
			    p[4].y = chains[i].residues[j].atoms[k].y;
			    p[4].z = chains[i].residues[j].atoms[k].z;
				}
			}
			for (k = 0, temp = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				name = chains[i].residues[j].atoms[k].name;
				if (name == "P" || name == "O1P" || name == "O2P") {
					temp++;
				}
			}
			
			if (temp ==	3 && j != 0) {
				x1 = p[0].x - p[2].x; y1 = p[0].y - p[2].y; z1 = p[0].z - p[2].z;
				x2 = p[0].x - q[1].x; y2 = p[0].y - q[1].y; z2 = p[0].z - q[1].z;
	
				r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
			}

			if ((temp != 3 || (temp == 3 && (r1 > 2 || r2 > 2))) && j != 0) {
				atom = new Atom[3];
				atom[0].name = "P"; atom[1].name = "O1P"; atom[2].name = "O2P";
				x1 = q[1].x - q[3].x; y1 = q[1].y - q[3].y; z1 = q[1].z - q[3].z;
				x2 = p[2].x - p[4].x; y2 = p[2].y - p[4].y; z2 = p[2].z - p[4].z;
				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				x1 /= r; y1 /= r; z1 /= r;
				r = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
				x2 /= r; y2 /= r; z2 /= r;
				x = (x1 + x2) / 2; y = (y1 + y2) / 2; z = (z1 + z2) / 2;
				x1 = p[2].x - q[1].x; y1 = p[2].y - q[1].y; z1 = p[2].z - q[1].z;
				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				x1 /= r; y1 /= r; z1 /= r;
				r = -(x1 * x + y1 * y + z1 * z);
				x1 *= r; y1 *= r; z1 *= r;
				x += x1; y += y1; z += z1;
				r = sqrt(x * x + y * y + z * z);
				x /= r; y /= r; z /= r;
				atom[0].x = 0.884 * x + (q[1].x + p[2].x) / 2; atom[0].y = 0.884 * y + (q[1].y + p[2].y) / 2; atom[0].z = 0.884 * z + (q[1].z + p[2].z) / 2;
//				atom[0].x = (q[1].x + p[2].x) / 2; atom[0].y = (q[1].y + p[2].y) / 2; atom[0].z = (q[1].z + p[2].z) / 2;
				x1 = p[2].x - q[1].x; y1 = p[2].y - q[1].y; z1 = p[2].z - q[1].z;
				x2 = x; y2 = y; z2 = z;
				a = new Matr_(2, 2);
				b = new Matr_(2, 1);
				a->data[0][0] = x1; a->data[0][1] = y1; a->data[1][0] = x2; a->data[1][1] = y2;
				b->data[0][0] = -z1; b->data[1][0] = -z2;
				c = a->inverse()->multiply(b);
				x1 = c->data[0][0]; y1 = c->data[1][0]; z1 = 1;
				delete a;
				delete b;
				delete c;
				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
				x1 /= r; y1 /= r; z1 /= r;
				x1 *= 1.25; y1 *= 1.25; z1 *= 1.25;
				if ((p[4].x - p[2].x) * x1 + (p[4].y - p[2].y) * y1 + (p[4].z - p[2].z) * z1 > 0) {
					atom[1].x = atom[0].x + x1 + 0.884 * x;
					atom[1].y = atom[0].y + y1 + 0.884 * y;
					atom[1].z = atom[0].z + z1 + 0.884 * z;
					atom[2].x = atom[0].x - x1 + 0.884 * x;
					atom[2].y = atom[0].y - y1 + 0.884 * y;
					atom[2].z = atom[0].z - z1 + 0.884 * z;
				} else {
					atom[1].x = atom[0].x - x1 + 0.884 * x;
					atom[1].y = atom[0].y - y1 + 0.884 * y;
					atom[1].z = atom[0].z - z1 + 0.884 * z;
					atom[2].x = atom[0].x + x1 + 0.884 * x;
					atom[2].y = atom[0].y + y1 + 0.884 * y;
					atom[2].z = atom[0].z + z1 + 0.884 * z;
				}
				for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
					name = chains[i].residues[j].atoms[k].name;
					if (name != "P" && name != "O1P" && name != "O2P") {
						va.push_back(chains[i].residues[j].atoms[k]);
					}
				}
				chains[i].residues[j].atoms.clear();
				chains[i].residues[j].atoms.push_back(atom[0]);
				chains[i].residues[j].atoms.push_back(atom[1]);
				chains[i].residues[j].atoms.push_back(atom[2]);
				for (k = 0; k < (int) va.size(); k++) {
					chains[i].residues[j].atoms.push_back(va[k]);
				}
			}
			if (q != NULL) {
				delete [] q;
			}
			q = p;
			va.clear();
		}
	}
}

void RNA::mutate(string seq) {
	int i, j, k, k1, k2, k3, k_1, k_2, k_3, flag, temp;
	string resName1, resName2, baseName, name;
	double x, y, z;

	for (i = 0, flag = 0, temp = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++, flag++) {
			if (flag >= (int) seq.size()) break;
			resName1 = chains[i].residues[j].name;
			if (seq[temp] != 'X') {
				resName2 = "";
				resName2 += seq[temp];
				if (chains[i].residues[j].name != resName2) {
					baseName = "";
					baseName += getenv("RNA");
					baseName += "base/";
					baseName += resName2;
					baseName += ".pdb";
					RNA *rna = new RNA(baseName);

					// move rna1 to align the C1' atom 
					for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
						if (chains[i].residues[j].atoms[k].name == "C1*") {
							k1 = k;
							x = chains[i].residues[j].atoms[k].x;
							y = chains[i].residues[j].atoms[k].y;
							z = chains[i].residues[j].atoms[k].z;
						}
						if (chains[i].residues[j].atoms[k].name == "C2") {
							k2 = k;
						}
					}
					for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
						name = chains[i].residues[j].atoms[k].name;
						if ((name == "N9" && (resName1 == "A" || resName1 == "G")) || (name == "N1" && (resName1 == "U" || resName1 == "C"))) {
							k3 = k;
							goto out1;
						}
					}
					out1:
					move(-x, -y, -z);

					// move rna2 to align the C1' atom 
					for (k = 0; k < (int) rna->chains[0].residues[0].atoms.size(); k++) {
						if (rna->chains[0].residues[0].atoms[k].name == "C1*") {
							x = rna->chains[0].residues[0].atoms[k].x;
							y = rna->chains[0].residues[0].atoms[k].y;
							z = rna->chains[0].residues[0].atoms[k].z;
							k_1 = k;
						}
						if (rna->chains[0].residues[0].atoms[k].name == "C2") {
							k_2 = k;
						}
					}
					for (k = 0; k < (int) rna->chains[0].residues[0].atoms.size(); k++) {
						name = rna->chains[0].residues[0].atoms[k].name;
						if ((name == "N9" && (resName2 == "A" || resName2 == "G")) || (name == "N1" && (resName2 == "U" || resName2 == "C"))) {
							k_3 = k;
							goto out2;
						}
					}
					out2:
					rna->move(-x, -y, -z);

					// rotate rna1 with x-axis to make N atom on the x-z plane
					Matr_ *m = new Matr_(3, 3);
					x = chains[i].residues[j].atoms[k3].x;
					y = chains[i].residues[j].atoms[k3].y;
					z = chains[i].residues[j].atoms[k3].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						rotate(m);
					}
					delete m;
				
					// rotate rna2 with x-axis to make N atom on the x-z plane
					m = new Matr_(3, 3);
					x = rna->chains[0].residues[0].atoms[k_3].x;
					y = rna->chains[0].residues[0].atoms[k_3].y;
					z = rna->chains[0].residues[0].atoms[k_3].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						rna->rotate(m);
					}
					delete m;
				
					// rotate rna1 with y-axis to make N atom on the x-axis
					m = new Matr_(3, 3);
					x = chains[i].residues[j].atoms[k3].x;
					y = chains[i].residues[j].atoms[k3].y;
					z = chains[i].residues[j].atoms[k3].z;
					if (x * x + z * z != 0) {
						m->data[0][0] = x / sqrt(x * x + z * z);  m->data[0][1] = 0; m->data[0][2] = z / sqrt(x * x + z * z);
						m->data[1][0] = 0;                        m->data[1][1] = 1; m->data[1][2] = 0;
						m->data[2][0] = -z / sqrt(x * x + z * z); m->data[2][1] = 0; m->data[2][2] = x / sqrt(x * x + z * z);
						rotate(m);
					}
					delete m;

					// rotate rna2 with y-axis to make N atom on the x-axis
					m = new Matr_(3, 3);
					x = rna->chains[0].residues[0].atoms[k_3].x;
					y = rna->chains[0].residues[0].atoms[k_3].y;
					z = rna->chains[0].residues[0].atoms[k_3].z;
					if (x * x + z * z != 0) {
						m->data[0][0] = x / sqrt(x * x + z * z);  m->data[0][1] = 0; m->data[0][2] = z / sqrt(x * x + z * z);
						m->data[1][0] = 0;                        m->data[1][1] = 1; m->data[1][2] = 0;
						m->data[2][0] = -z / sqrt(x * x + z * z); m->data[2][1] = 0; m->data[2][2] = x / sqrt(x * x + z * z);
						rna->rotate(m);
					}
					delete m;

					// rotate rna1 with x-axis to make C2 atom on the x-z plane

					m = new Matr_(3, 3);
					x = chains[i].residues[j].atoms[k2].x;
					y = chains[i].residues[j].atoms[k2].y;
					z = chains[i].residues[j].atoms[k2].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						rotate(m);
					}
					delete m;

					// rotate rna2 with x-axis to make C2 atom on the x-z plane
					m = new Matr_(3, 3);
					x = rna->chains[0].residues[0].atoms[k_2].x;
					y = rna->chains[0].residues[0].atoms[k_2].y;
					z = rna->chains[0].residues[0].atoms[k_2].z;
					if (y * y + z * z != 0) {
						m->data[0][0] = 1; m->data[0][1] = 0;                       m->data[0][2] = 0;
						m->data[1][0] = 0; m->data[1][1] = z / sqrt(y * y + z * z); m->data[1][2] = -y / sqrt(y * y + z * z);
						m->data[2][0] = 0; m->data[2][1] = y / sqrt(y * y + z * z); m->data[2][2] = z / sqrt(y * y + z * z);
						rna->rotate(m);
					}
					delete m;

					// change base
					vector<Atom> v;
					for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
						name = chains[i].residues[j].atoms[k].name;
						if (name == "P" || name == "O1P" || name == "O2P" || name[name.size() - 1] == '*') {
							v.push_back(chains[i].residues[j].atoms[k]);
						}
					}
					chains[i].residues[j].atoms.clear();
					for (k = 0; k < (int) v.size(); k++) {
						chains[i].residues[j].atoms.push_back(v[k]);
					}
					for (k = 0; k < (int) rna->chains[0].residues[0].atoms.size(); k++) {
						if (rna->chains[0].residues[0].atoms[k].name != "C1*") {
							chains[i].residues[j].atoms.push_back(rna->chains[0].residues[0].atoms[k]);
						}
					}
					chains[i].residues[j].name = resName2;
					v.clear();
					delete rna;
				}
			}
			temp++;
		}
	}
}

void RNA::rotateByX(double alpha) {
	int i, j, k;
	double c, s;
	double x, y, z;

	c = cos(alpha); s = sin(alpha);
	for (i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				x = chains[i].residues[j].atoms[k].x;
				y = chains[i].residues[j].atoms[k].y;
				z = chains[i].residues[j].atoms[k].z;
				chains[i].residues[j].atoms[k].x = x;
				chains[i].residues[j].atoms[k].y = c * y - s * z;
				chains[i].residues[j].atoms[k].z = s * y + c * z;
			}
		}
	}
}

void RNA::rotateByZ(double beta) {
	int i, j, k;
	double c, s;
	double x, y, z;

	c = cos(beta); s = sin(beta);
	for (i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				x = chains[i].residues[j].atoms[k].x;
				y = chains[i].residues[j].atoms[k].y;
				z = chains[i].residues[j].atoms[k].z;
				chains[i].residues[j].atoms[k].x = c * x - s * y;
				chains[i].residues[j].atoms[k].y = s * x + c * y;
				chains[i].residues[j].atoms[k].z = z;
			}
		}
	}
}

