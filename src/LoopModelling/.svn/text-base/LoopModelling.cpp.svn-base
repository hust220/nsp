#include "LoopModelling.h"

void LoopModelling::init() {
	/* initialize bound matrix*/
	bound.resize(len, len);
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			if (i == j) {
				bound(i, j) = 0;
			} else if (j - i <= 12) {
				bound(j, i) = par[(j - i - 1) * 12 + 2 * (i % 6)];
				bound(i, j) = par[(j - i - 1) * 12 + 2 * (i % 6) + 1];
			} else {
				bound(j, i) = 7;
				bound(i, j) = 999;
			}
		}
	}

	// base pairs
	Mol2D mol2d(ss);
	set_base_pairs(mol2d.head);

	string str = ss;
	transform(str.begin(), str.end(), str.begin(), [](const char &c) {
		map<char, char> temp_map{{'(', '['}, {')', ']'}, {'[', '('}, {']', ')'}, {'.', '.'}, {'&', '&'}};
		return temp_map[c];
	});
	Mol2D mol2d_2(str);
	set_base_pairs(mol2d_2.head);

	/* chirality */
	int flag = 0;
	for (int i = 3; i < len;) {
		if (ss[i / 6] == '(' && i % 6 == 5 && i != 5) {
			i += 4;
		} else {
			i++;
		}
		flag++;
	}
	delete chir;
	chir = new Matr_(flag, 5);
	for (int i = 3, flag = 0; i < len;) {
		chir->data[flag][0] = i - 3;
		chir->data[flag][1] = i - 2;
		chir->data[flag][2] = i - 1;
		chir->data[flag][3] = i;
		chir->data[flag][4] = chirality[(i - 2) % 6];
		if (ss[i / 6] == '(' && i % 6 == 5 && i != 5) {
			i += 4;
		} else {
			i++;
		}
		flag++;
	}

}

void LoopModelling::set_base_pairs(loop *src) {
	if (src == NULL) {
		return;
	} else {
		if (src->s.head != NULL) {
			int helix_len = src->s.getLen();
			int orig_pos_chain1 = (src->s.head->res1.num - 1) * atom_nums_per_nuc;
			int orig_pos_chain2 = (src->s.head->res2.num - helix_len) * atom_nums_per_nuc;
			int atom_nums_per_chain = atom_nums_per_nuc * helix_len;
cout << "helix_len: " << helix_len << endl;
cout << "orig_pos_chain1: " << orig_pos_chain1 << endl;
cout << "orig_pos_chain2: " << orig_pos_chain2 << endl;
cout << "atom_nums_per_chain: " << atom_nums_per_chain << endl;
cout << endl;
			for (int i: {0, 1}) {
				for (int j: {0, 1}) {
					for (int ii = 0; ii < atom_nums_per_chain; ii++) {
						for (int jj = 0; jj < atom_nums_per_chain; jj++) {
							if (i == j && ii == jj) {
								continue;
							}
							bound(orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * i + ii,
							      orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * j + jj)
							= helix_par(i * (helix_par.rows() - atom_nums_per_chain) + ii, j * (helix_par.rows() - atom_nums_per_chain) + jj) -
							  err_radius;
							bound(orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * j + jj,
							      orig_pos_chain1 + (orig_pos_chain2 - orig_pos_chain1) * i + ii)
							= helix_par(i * (helix_par.rows() - atom_nums_per_chain) + ii, j * (helix_par.rows() - atom_nums_per_chain) + jj) +
							  err_radius;
						}
					}
				}
			}
			
		}
		set_base_pairs(src->son);
		set_base_pairs(src->brother);
	}
}

LoopModelling::~LoopModelling() {
	delete chir;
	delete [] type;
	delete [] par;
//	delete [] par2;
}

RNA *LoopModelling::run() {
	backbone = dg();
	Point *p = NULL;
	RNA *rna = new RNA;
	Chain chain;
	rna->chains.push_back(chain);
	Point *o = NULL;
	for (int i = 0; i < resLen; i++) {
		delete [] p;
		p = new Point[6];
		for (int j = 0; j < 6; j++) {
			p[j].x = backbone(i * 6 + j, 0);
			p[j].y = backbone(i * 6 + j, 1);
			p[j].z = backbone(i * 6 + j, 2);
		}
		double chi = (rand() % 1000) / 1000. * 360;
		Residue *res = buildNuc(p, o, chi, type[i]);
		rna->chains[0].residues.push_back(*res);

		delete o;
		if (i != 0 && ss[i] == '(') {
			o = NULL;
		} else {
			o = new Point;
			o->x = backbone(i * 6 + 5, 0);
			o->y = backbone(i * 6 + 5, 1);
			o->z = backbone(i * 6 + 5, 2);
		}
	}
	return rna;
}

double LoopModelling::at(MatrixXf &bound, int i, int j, int len) {
    if (j >= len) {
        j -= len;
        return at(bound, j, i, len);
    }
    if (i >= len) {
        i -= len;
        return at(bound, j, i, len);
    }
    return bound(i, j);
}

void LoopModelling::assign(MatrixXf &bound, int i, int j, double d, int len) {
    if (j >= len) {
        j -= len;
        assign(bound, j, i, d, len);
        return;
    }
    if (i >= len) {
        i -= len;
        assign(bound, j, i, d, len);
        return;
    }
    bound(i, j) = d;
}

Residue *LoopModelling::buildNuc(Point *p1, Point *o, double chi, int type) {
	Point *p2 = new Point[7];
	string lib = getenv("RNA");
	string f_name = lib + "sugar";
	ifstream ifile(f_name.c_str());
	for (int i = 0; i < 7; i++) {
		ifile >> p2[i].x >> p2[i].y >> p2[i].z;
	}
	ifile.close();

	// move p2 to origin
	double x_ = p2[2].x;
	double y_ = p2[2].y;
	double z_ = p2[2].z;
	for (int i = 0; i < 7; i++) {
		p2[i].x -= x_;
		p2[i].y -= y_;
		p2[i].z -= z_;
	}

	// target vector
	Point target(p1[3].x - p1[4].x, p1[3].y - p1[4].y, p1[3].z - p1[4].z);

	Point::coincide(p2, 7, p2[0], target);

	// normal vector
	/*
	Point *n1 = Point::normalVector(p1[4], p1[3], p1[2]);
	Point *n2 = Point::normalVector(p2[2], p2[0], p2[1]);
	double a1 = Point::angle(n1, origin, n2);
	*/
	Point *origin = new Point;
	Point *a1_1 = new Point(p1[2].x - p1[4].x, p1[2].y - p1[4].y, p1[2].z - p1[4].z);
	Point *a1_2 = new Point(p1[4].x - p1[4].x, p1[4].y - p1[4].y, p1[4].z - p1[4].z);
	Point *a1_3 = new Point(p1[3].x - p1[4].x, p1[3].y - p1[4].y, p1[3].z - p1[4].z);
	Point *a1_4 = new Point(p2[1].x, p2[1].y, p2[1].z);
	double a1 = Point::dihedral(a1_1, a1_2, a1_3, a1_4);
	delete a1_1; delete a1_2; delete a1_3; delete a1_4;
	double delta_a = 120 - a1;
	Point::rotate(p2, 7, p2[2], p2[0], delta_a);
	for (int i = 0; i < 7; i++) {
		p2[i].x += p1[4].x;
		p2[i].y += p1[4].y;
		p2[i].z += p1[4].z;
	}

	Point *p, *p3;
	int lent, len3;
	if (type == 0) {
		len3 = 11;
		if (o == NULL) {
			lent = 19;
			p = new Point[19];
		} else {
			lent = 22;
			p = new Point[22];
		}
		p3 = new Point[11];
		f_name = lib + "A";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 11; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 1) {
		len3 = 9;
		if (o == NULL) {
			lent = 17;
			p = new Point[17];
		} else {
			lent = 20;
			p = new Point[20];
		}
		p3 = new Point[9];
		f_name = lib + "U";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 9; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 2) {
		len3 = 12;
		if (o == NULL) {
			lent = 20;
			p = new Point[20];
		} else {
			lent = 23;
			p = new Point[23];
		}
		p3 = new Point[12];
		f_name = lib + "G";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 12; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 3) {
		len3 = 9;
		if (o == NULL) {
			lent = 17;
			p = new Point[17];
		} else {
			lent = 20;
			p = new Point[20];
		}
		p3 = new Point[9];
		f_name = lib + "C";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 9; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else {
		cerr << "LoopModelling2::buildNuc error! The type must be one of 0, 1, 2, 3" << endl;
		exit(1);
	}

	// move p3 to origin
	x_ = p3[0].x;
	y_ = p3[0].y;
	z_ = p3[0].z;
	for (int i = 0; i < len3; i++) {
		p3[i].x -= x_;
		p3[i].y -= y_;
		p3[i].z -= z_;
	}

	// target vector
	target.x = p2[6].x - p2[5].x;
	target.y = p2[6].y - p2[5].y;
	target.z = p2[6].z - p2[5].z;

	Point::coincide(p3, len3, p3[1], target);

	// normalVector
	Point *n1 = Point::normalVector(p2[1], p2[5], p2[3]);

	// dihedral
	double dih = Point::dihedral(n1, &(p3[0]), &(p3[1]), &(p3[2]));
	double delta_dih = chi - dih;

	// rotate
	Point::rotate(p3, len3, *origin, p3[1], delta_dih);
	for (int i = 0; i < len3; i++) {
		p3[i].x += p2[5].x;
		p3[i].y += p2[5].y;
		p3[i].z += p2[5].z;
	}

	if (o == NULL) {
		p[0].x = p1[1].x; p[0].y = p1[1].y; p[0].z = p1[1].z;
		p[1].x = p1[2].x; p[1].y = p1[2].y; p[1].z = p1[2].z;
		p[2].x = p1[3].x; p[2].y = p1[3].y; p[2].z = p1[3].z;
		p[3].x = p2[1].x; p[3].y = p2[1].y; p[3].z = p2[1].z;
		p[4].x = p1[4].x; p[4].y = p1[4].y; p[4].z = p1[4].z;
		p[5].x = p1[5].x; p[5].y = p1[5].y; p[5].z = p1[5].z;
		p[6].x = p2[3].x; p[6].y = p2[3].y; p[6].z = p2[3].z;
		p[7].x = p2[4].x; p[7].y = p2[4].y; p[7].z = p2[4].z;
		for (int i = 0; i < len3; i++) {
			p[i + 8].x = p3[i].x;
			p[i + 8].y = p3[i].y;
			p[i + 8].z = p3[i].z;
		}
	} else {
		Point o3(o->x - p1[0].x, o->y - p1[0].y, o->z - p1[0].z);
		Point o5(p1[1].x - p1[0].x, p1[1].y - p1[0].y, p1[1].z - p1[0].z);
		double ratio = 1.5 / sqrt(o3.x * o3.x + o3.y * o3.y + o3.z * o3.z);
		Point o1p(ratio * o3.x, ratio * o3.y, ratio * o3.z);
		Point o2p(ratio * o3.x, ratio * o3.y, ratio * o3.z);
		Point *spindle1 = Point::normalVector(o5, *origin, o3);
		Point::rotate(&o1p, 1, *origin, *spindle1, 109);
		Point::rotate(&o2p, 1, *origin, *spindle1, 109);
		Point::rotate(&o1p, 1, *origin, o3, 120);
		Point::rotate(&o2p, 1, *origin, o3, 240);
	
		p[0].x = p1[0].x; p[0].y = p1[0].y; p[0].z = p1[0].z;
		p[1].x = o1p.x + p1[0].x; p[1].y = o1p.y + p1[0].y; p[1].z = o1p.z + p1[0].z;
		p[2].x = o2p.x + p1[0].x; p[2].y = o2p.y + p1[0].y; p[2].z = o2p.z + p1[0].z;
		p[3].x = p1[1].x; p[3].y = p1[1].y; p[3].z = p1[1].z;
		p[4].x = p1[2].x; p[4].y = p1[2].y; p[4].z = p1[2].z;
		p[5].x = p1[3].x; p[5].y = p1[3].y; p[5].z = p1[3].z;
		p[6].x = p2[1].x; p[6].y = p2[1].y; p[6].z = p2[1].z;
		p[7].x = p1[4].x; p[7].y = p1[4].y; p[7].z = p1[4].z;
		p[8].x = p1[5].x; p[8].y = p1[5].y; p[8].z = p1[5].z;
		p[9].x = p2[3].x; p[9].y = p2[3].y; p[9].z = p2[3].z;
		p[10].x = p2[4].x; p[10].y = p2[4].y; p[10].z = p2[4].z;
		for (int i = 0; i < len3; i++) {
			p[i + 11].x = p3[i].x;
			p[i + 11].y = p3[i].y;
			p[i + 11].z = p3[i].z;
		}
	}

	Residue *residue = new Residue(p, lent, type);
	return residue;

}

