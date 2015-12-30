#ifndef ADDPHOS_H
#define ADDPHOS_H

#include <pdb/util.h>

namespace jian {

namespace nuc3d {

class AddPhos {
public:
    void operator() (Model &model) {
    //	int i, j, k, temp;
    //	vector<Atom> va;
    //	Point *p = NULL, *q = NULL;
    //	Matr_ *a, *b, *c;
    //	double x1, x2, x, y1, y2, y, z1, z2, z;
    //	string name;
    //	Atom *atom;
    //	double r, r1, r2;
    //
    //	for (i = 0; i < (int) model.chains.size(); i++) {
    //		for (j = 0; j < (int) model.chains[i].residues.size(); j++) {
    //			p = new Point[5];
    //			for (k = 0; k < (int) model.chains[i].residues[j].atoms.size(); k++) {
    //				if (model.chains[i].residues[j].atoms[k].name == "P") {
    //		      p[0].x = model.chains[i].residues[j].atoms[k].x;
    //					p[0].y = model.chains[i].residues[j].atoms[k].y;
    //			    p[0].z = model.chains[i].residues[j].atoms[k].z;
    //			  } else if (model.chains[i].residues[j].atoms[k].name == "O3*") {
    //		     	p[1].x = model.chains[i].residues[j].atoms[k].x;
    //				 	p[1].y = model.chains[i].residues[j].atoms[k].y;
    //					p[1].z = model.chains[i].residues[j].atoms[k].z;
    //			  } else if (model.chains[i].residues[j].atoms[k].name == "O5*") {
    //			    p[2].x = model.chains[i].residues[j].atoms[k].x;
    //			    p[2].y = model.chains[i].residues[j].atoms[k].y;
    //			    p[2].z = model.chains[i].residues[j].atoms[k].z;
    //				} else if (model.chains[i].residues[j].atoms[k].name == "C3*") {
    //			    p[3].x = model.chains[i].residues[j].atoms[k].x;
    //		      p[3].y = model.chains[i].residues[j].atoms[k].y;
    //					p[3].z = model.chains[i].residues[j].atoms[k].z;
    //			  } else if (model.chains[i].residues[j].atoms[k].name == "C5*") {
    //			    p[4].x = model.chains[i].residues[j].atoms[k].x;
    //			    p[4].y = model.chains[i].residues[j].atoms[k].y;
    //			    p[4].z = model.chains[i].residues[j].atoms[k].z;
    //				}
    //			}
    //			for (k = 0, temp = 0; k < (int) model.chains[i].residues[j].atoms.size(); k++) {
    //				name = model.chains[i].residues[j].atoms[k].name;
    //				if (name == "P" || name == "O1P" || name == "O2P") {
    //					temp++;
    //				}
    //			}
    //			
    //			if (temp ==	3 && j != 0) {
    //				x1 = p[0].x - p[2].x; y1 = p[0].y - p[2].y; z1 = p[0].z - p[2].z;
    //				x2 = p[0].x - q[1].x; y2 = p[0].y - q[1].y; z2 = p[0].z - q[1].z;
    //	
    //				r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    //				r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    //			}
    //
    //			if ((temp != 3 || (temp == 3 && (r1 > 2 || r2 > 2))) && j != 0) {
    //				atom = new Atom[3];
    //				atom[0].name = "P"; atom[1].name = "O1P"; atom[2].name = "O2P";
    //				x1 = q[1].x - q[3].x; y1 = q[1].y - q[3].y; z1 = q[1].z - q[3].z;
    //				x2 = p[2].x - p[4].x; y2 = p[2].y - p[4].y; z2 = p[2].z - p[4].z;
    //				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    //				x1 /= r; y1 /= r; z1 /= r;
    //				r = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    //				x2 /= r; y2 /= r; z2 /= r;
    //				x = (x1 + x2) / 2; y = (y1 + y2) / 2; z = (z1 + z2) / 2;
    //				x1 = p[2].x - q[1].x; y1 = p[2].y - q[1].y; z1 = p[2].z - q[1].z;
    //				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    //				x1 /= r; y1 /= r; z1 /= r;
    //				r = -(x1 * x + y1 * y + z1 * z);
    //				x1 *= r; y1 *= r; z1 *= r;
    //				x += x1; y += y1; z += z1;
    //				r = sqrt(x * x + y * y + z * z);
    //				x /= r; y /= r; z /= r;
    //				atom[0].x = 0.884 * x + (q[1].x + p[2].x) / 2; atom[0].y = 0.884 * y + (q[1].y + p[2].y) / 2; atom[0].z = 0.884 * z + (q[1].z + p[2].z) / 2;
    ////				atom[0].x = (q[1].x + p[2].x) / 2; atom[0].y = (q[1].y + p[2].y) / 2; atom[0].z = (q[1].z + p[2].z) / 2;
    //				x1 = p[2].x - q[1].x; y1 = p[2].y - q[1].y; z1 = p[2].z - q[1].z;
    //				x2 = x; y2 = y; z2 = z;
    //				a = new Matr_(2, 2);
    //				b = new Matr_(2, 1);
    //				a->data[0][0] = x1; a->data[0][1] = y1; a->data[1][0] = x2; a->data[1][1] = y2;
    //				b->data[0][0] = -z1; b->data[1][0] = -z2;
    //				c = a->inverse()->multiply(b);
    //				x1 = c->data[0][0]; y1 = c->data[1][0]; z1 = 1;
    //				delete a;
    //				delete b;
    //				delete c;
    //				r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    //				x1 /= r; y1 /= r; z1 /= r;
    //				x1 *= 1.25; y1 *= 1.25; z1 *= 1.25;
    //				if ((p[4].x - p[2].x) * x1 + (p[4].y - p[2].y) * y1 + (p[4].z - p[2].z) * z1 > 0) {
    //					atom[1].x = atom[0].x + x1 + 0.884 * x;
    //					atom[1].y = atom[0].y + y1 + 0.884 * y;
    //					atom[1].z = atom[0].z + z1 + 0.884 * z;
    //					atom[2].x = atom[0].x - x1 + 0.884 * x;
    //					atom[2].y = atom[0].y - y1 + 0.884 * y;
    //					atom[2].z = atom[0].z - z1 + 0.884 * z;
    //				} else {
    //					atom[1].x = atom[0].x - x1 + 0.884 * x;
    //					atom[1].y = atom[0].y - y1 + 0.884 * y;
    //					atom[1].z = atom[0].z - z1 + 0.884 * z;
    //					atom[2].x = atom[0].x + x1 + 0.884 * x;
    //					atom[2].y = atom[0].y + y1 + 0.884 * y;
    //					atom[2].z = atom[0].z + z1 + 0.884 * z;
    //				}
    //				for (k = 0; k < (int) model.chains[i].residues[j].atoms.size(); k++) {
    //					name = model.chains[i].residues[j].atoms[k].name;
    //					if (name != "P" && name != "O1P" && name != "O2P") {
    //						va.push_back(model.chains[i].residues[j].atoms[k]);
    //					}
    //				}
    //				model.chains[i].residues[j].atoms.clear();
    //				model.chains[i].residues[j].atoms.push_back(atom[0]);
    //				model.chains[i].residues[j].atoms.push_back(atom[1]);
    //				model.chains[i].residues[j].atoms.push_back(atom[2]);
    //				for (k = 0; k < (int) va.size(); k++) {
    //					model.chains[i].residues[j].atoms.push_back(va[k]);
    //				}
    //			}
    //			if (q != NULL) {
    //				delete [] q;
    //			}
    //			q = p;
    //			va.clear();
    //		}
    //	}
    }};

} /// namespace nuc3d

} /// namespace jian

#endif

