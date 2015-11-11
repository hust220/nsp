#include "P2B.h"

using namespace jian;

P2B::P2B(RNA *rna) {
	len = rna->getLen();
	for (int i = 0; i < len; i++) {
		ss += '.';
	}
	b.resize(3, len);
	n.resize(3, len);
	d.resize(3, len);
	for (int i = 0, index = 0; i < (int) rna->chains.size(); i++) {
		for (int j = 0; j < (int) rna->chains[i].residues.size(); j++, index++) {
			Point *c2 = new Point;
			Point *c4 = new Point;
			Point *c6 = new Point;
			for (int k = 0; k < (int) rna->chains[i].residues[j].atoms.size(); k++) {
				if (rna->chains[i].residues[j].atoms[k].name == "C2") {
					c2->x = rna->chains[i].residues[j].atoms[k].x;
					c2->y = rna->chains[i].residues[j].atoms[k].y;
					c2->z = rna->chains[i].residues[j].atoms[k].z;
				} else if (rna->chains[i].residues[j].atoms[k].name == "C4") {
					c4->x = rna->chains[i].residues[j].atoms[k].x;
					c4->y = rna->chains[i].residues[j].atoms[k].y;
					c4->z = rna->chains[i].residues[j].atoms[k].z;
				} else if (rna->chains[i].residues[j].atoms[k].name == "C6") {
					c6->x = rna->chains[i].residues[j].atoms[k].x;
					c6->y = rna->chains[i].residues[j].atoms[k].y;
					c6->z = rna->chains[i].residues[j].atoms[k].z;
				}
			}
			b(0, index) = (c2->x + c4->x + c6->x) / 3.;
			b(1, index) = (c2->y + c4->y + c6->y) / 3.;
			b(2, index) = (c2->z + c4->z + c6->z) / 3.;
			Point *p = Point::normalVector(c2, c4, c6);
			n(0, index) = p->x;
			n(1, index) = p->y;
			n(2, index) = p->z;
			delete p;
			string name = rna->chains[i].residues[j].name;
			if (name == "A" || name == "G") {
				d(0, index) = c2->x - c6->x;
				d(1, index) = c2->y - c6->y;
				d(2, index) = c2->z - c6->z;
			} else {
				d(0, index) = c4->x - c2->x;
				d(1, index) = c4->y - c2->y;
				d(2, index) = c4->z - c2->z;
			}
			delete c2; delete c4; delete c6;
		}
	}
}

void P2B::analyze() {
	MatrixXf bb(len, len);
	MatrixXf nn(len, len);
	MatrixXf nb(len, len);
	MatrixXf dd(len, len);
	for (int i = 0; i < len; i++) {
		for (int j = i; j < len; j++) {
			if (i == j) {
				bb(i, j) = 0;
				nn(i, j) = 0;
				nb(i, j) = 0;
				dd(i, j) = 0;
			} else {
				bb(i, j) = sqrt((b(0, i) - b(0, j)) * (b(0, i) - b(0, j)) + (b(1, i) - b(1, j)) * (b(1, i) - b(1, j)) + (b(2, i) - b(2, j)) * (b(2, i) - b(2, j)));
				bb(j, i) = 0;

				Point *p1 = new Point(b(0, i) + n(0, i), b(1, i) + n(1, i), b(2, i) + n(2, i));
				Point *p2 = new Point(b(0, i), b(1, i), b(2, i));
				Point *p3 = new Point(b(0, j), b(1, j), b(2, j));
				Point *p4 = new Point(b(0, j) + n(0, j), b(1, j) + n(1, j), b(2, j) + n(2, j));
				double dih = Point::dihedral(p1, p2, p3, p4);
				nn(i, j) = dih;
				nn(j, i) = 0;
				delete p4;
				
				p1->x = (n(0, i) + n(0, j)) / 2. + b(0, i);
				p1->y = (n(1, i) + n(1, j)) / 2. + b(1, i);
				p1->z = (n(2, i) + n(2, j)) / 2. + b(2, i);
				p2->x = b(0, i);
				p2->y = b(1, i);
				p2->z = b(2, i);
				p3->x = b(0, j);
				p3->y = b(1, j);
				p3->z = b(2, j);
				double ang = Point::angle(p1, p2, p3);
				nb(i, j) = ang;
				nb(j, i) = 0;

				p1->x = d(0, i);
				p1->y = d(1, i);
				p1->z = d(2, i);
				p2->x = 0;
				p2->y = 0;
				p2->z = 0;
				p3->x = d(0, j);
				p3->y = d(1, j);
				p3->z = d(2, j);
				ang = Point::angle(p1, p2, p3);
				dd(i, j) = ang;
				dd(j, i) = 0;
				delete p1; delete p2; delete p3;
			}
		}
	}
	cout << bb << endl;
	cout << nn << endl;
	cout << nb << endl;
	cout << dd << endl;
}

string P2B::run() {
	for (int i = 0; i < len; i++) {
		for (int j = i + 1; j < len; j++) {
			double dist = sqrt((b(0, i) - b(0, j)) * (b(0, i) - b(0, j)) + (b(1, i) - b(1, j)) * (b(1, i) - b(1, j)) + (b(2, i) - b(2, j)) * (b(2, i) - b(2, j)));
			Point *p1 = new Point(b(0, i) + n(0, i), b(1, i) + n(1, i), b(2, i) + n(2, i));
			Point *p2 = new Point(b(0, i), b(1, i), b(2, i));
			Point *p3 = new Point(b(0, j), b(1, j), b(2, j));
			Point *p4 = new Point(b(0, j) + n(0, j), b(1, j) + n(1, j), b(2, j) + n(2, j));
			double dih = Point::dihedral(p1, p2, p3, p4);
			delete p4;
				
			p1->x = (n(0, i) + n(0, j)) / 2. + b(0, i);
			p1->y = (n(1, i) + n(1, j)) / 2. + b(1, i);
			p1->z = (n(2, i) + n(2, j)) / 2. + b(2, i);
			p2->x = b(0, i);
			p2->y = b(1, i);
			p2->z = b(2, i);
			p3->x = b(0, j);
			p3->y = b(1, j);
			p3->z = b(2, j);
			double ang = Point::angle(p1, p2, p3);

			p1->x = d(0, i);
			p1->y = d(1, i);
			p1->z = d(2, i);
			p2->x = 0;
			p2->y = 0;
			p2->z = 0;
			p3->x = d(0, j);
			p3->y = d(1, j);
			p3->z = d(2, j);
			double ang2 = Point::angle(p1, p2, p3);

			delete p1; delete p2; delete p3;
			if ((dist < 5.8 && dist > 5.1) && (dih > 330 || dih < 30) && (ang < 110 && ang > 70) && (ang2 > 150)) {
			 	ss[i] = '(';
				ss[j] = ')';
			}
		}
	}
	return ss;
}



