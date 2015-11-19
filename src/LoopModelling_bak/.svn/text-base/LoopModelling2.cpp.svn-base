#include "LoopModelling2.h"

LoopModelling2::LoopModelling2(string ss, string seq) {
	this->seq = seq;
	this->ss = ss;
	if (ss.size() != seq.size()) {
		cerr << "Sequence's length is not equal to second structure's length!" << endl;
		exit(1);
	}
	len = ss.size() * 6;
	d = new double[len - 1];
	for (int i = 0; i < len - 1; i++) {
		if (i % 6 == 0 || i % 6 == 5) {
			d[i] = 1.6;
		} else if (i % 6 == 1 || i % 6 == 4) {
			d[i] = 1.4;
		} else {
			d[i] = 1.5;
		}
	}
	for (int i = 0; i < ss.size(); i++) {
		if (ss[i] == '(' && i != 0) {
			d[i * 6 - 1] = 17;
		}
	}
	a = new double[len - 2];
	for (int i = 0; i < len - 2; i++) {
		a[i] = 108;
	}
	t = new double[len - 3];
	for (int i = 0; i < len - 3; i++) {
		t[i] = 0;
	}
	string f = getenv("RNA");
	string f_info = f + "dinuctor/info";
	ifstream ifile(f_info.c_str());
	string temp[16] = {"AA", "AU", "AG", "AC", "UA", "UU", "UG", "UC", "GA", "GU", "GG", "GC", "CA", "CU", "CG", "CC"};
	for (int i = 0; i < 16; i++) {
		string f_dinuc = f + "dinuctor/";
		f_dinuc += temp[i];
		ifstream if_dinuc(f_dinuc.c_str());

		ifile >> dinucQuant[i];
		dinuctor[i] = new double *[dinucQuant[i]];
		for (int j = 0; j < dinucQuant[i]; j++) {
			dinuctor[i][j] = new double[9];
			ifstream ifile2();
			for (int k = 0; k < 9; k++) {
				if_dinuc >> dinuctor[i][j][k];
			}
		}
		if_dinuc.close();
	}
	ifile.close();
}

RNA *LoopModelling2::run() {
	double dist = getDist();
	srand((unsigned int)time(0));
	int type[ss.size()];
	for (int i = 0; i < ss.size(); i++) {
		if (seq[i] == 'A') {
			type[i] = 0;
		} else if (seq[i] == 'U') {
			type[i] = 1;
		} else if (seq[i] == 'G') {
			type[i] = 2;
		} else {
			type[i] = 3;
		}
	}
	while (abs(dist - 17) > 0.5) {
		int n = int((rand() % 1000) / 1000. * (ss.size() - 1));
		if (ss[n] == '(' && n != 0) continue;
		int a = type[n] * 4 + type[n + 1];
		int temp = int((rand() % 1000) / 1000. * dinucQuant[a]);
		for (int i = 0; i < 9; i++) {
			t[6 * n + i] = dinuctor[a][temp][i];
		}
		dist = getDist();
	}
	cout << dist << endl;
	/*
	for (int i = 0; i < len; i++) {
		cout << c[i].x << '\t' << c[i].y << '\t' << c[i].z << endl;
	}
	*/
	RNA *rna = new RNA;
	Chain chain;
	rna->chains.push_back(chain);
	rna->name = "";
	rna->len = seq.size();
	Residue *residue;
	Point *o, *n;
	double *r_l, *r_a, *r_t;
	o = new Point[3];
	n = new Point[3];
	o[0].x = 0; o[0].y = 0; o[0].z = 0;
	o[1].x = 0; o[1].y = 0; o[1].z = 1.6;
	o[2].x = 0; o[2].y = 1.38564; o[2].z = 2.4;
	for (int i = 0; i < ss.size(); i++) {
		r_l = new double[6];
		r_a = new double[6];
		r_t = new double[6];
		for (int j = 0; j < 6; j++) {
			if (i == len - 1 && j == 5) {
				r_l[j] = d[0];
			} else {
				r_l[j] = d[6 * i + j];
			}
			if (i == len - 1 && j == 4) {
				r_a[j] = 120;
			} else if (i == len - 1 && j == 5) {
				r_a[j] = 120;
			} else {
				r_a[j] = a[6 * i + j];
			}
			if (i == 0 && j == 0) {
				r_t[j] = 0;
			} else if (i == len - 1 && j == 4) {
				r_t[j] = 0;
			} else if (i == len - 1 && j == 5) {
				r_t[j] = 0;
			} else {
				r_t[j] = t[6 * i + j - 1];
			}
		}
		double r_chi = (rand() % 1000) / 1000. * 360;
		residue = buildNuc(r_l, r_a, r_t, r_chi, o, n, type[i]);
		rna->chains[0].residues.push_back(*residue);
		delete residue;
		delete [] r_l;
		delete [] r_a;
		delete [] r_t;
		o[0].x = n[0].x; o[0].y = n[0].y; o[0].z = n[0].z;
		o[1].x = n[1].x; o[1].y = n[1].y; o[1].z = n[1].z;
		o[2].x = n[2].x; o[2].y = n[2].y; o[2].z = n[2].z;
	}
	return rna;
}

double LoopModelling2::getDist() {
	c = Point::t2c(d, a, t, len);
	return c[0].dist(c[len - 1]);
}

Residue *LoopModelling2::buildNuc(double *l, double *a, double *t, double chi, Point *o, Point *n, int type) {
	Point *p1 = Point::t2c(l, a, t, o, n, 6);
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

	/*
	// rotate p2 with z to let p2[0]-p2[2] and p1[3]-p1[4] overlap together
	double r = sqrt(p2[0].x * p2[0].x + p2[0].y * p2[0].y);
	if (r != 0) {
		double c = p2[0].y / r;
		double s = p2[0].x / r;
		for (int i = 0; i < 7; i++) {
			double x_ = c * p2[i].x - s * p2[i].y;
			double y_ = s * p2[i].y + c * p2[i].y;
			double z_ = p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}

	// rotate p2 with x to let p2[0]-p2[2] and p1[3]-p1[4] overlap together
	r = sqrt(p2[0].y * p2[0].y + p2[0].z * p2[0].z);
	if (r != 0) {
		double c = p2[0].z / r;
		double s = p2[0].y / r;
		for (int i = 0; i < 7; i++) {
			double x_ = p2[i].x;
			double y_ = c * p2[i].y - s * p2[i].z;
			double z_ = s * p2[i].y + c * p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}

	// rotate p2 with x
	r = sqrt(target.x * target.x + target.y * target.y + target.z * target.z);
	if (r != 0) {
		double c = target.z / r;
		double s = -sqrt(target.x * target.x + target.y * target.y) / r;
		for (int i = 0; i < 7; i++) {
			double x_ = p2[i].x;
			double y_ = c * p2[i].y - s * p2[i].z;
			double z_ = s * p2[i].y + c * p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}

	// rotate p2 with z
	r = sqrt(target.x * target.x + target.y * target.y);
	if (r != 0) {
		double c = target.y / r;
		double s = -target.x / r;
		for (int i = 0; i < 7; i++) {
			double x_ = c * p2[i].x - s * p2[i].y;
			double y_ = s * p2[i].y + c * p2[i].y;
			double z_ = p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}
	*/

	// normal vector
	Point *n1 = Point::normalVector(p1[4], p1[3], p1[2]);
	Point *n2 = Point::normalVector(p2[2], p2[0], p2[1]);
	Point *origin = new Point;
	double a1 = Point::angle(n1, origin, n2);
	double delta_a = 120 - a1;

	/* rotate p2 with p1[4]-p1[3] */
	Point::rotate(p2, 7, p2[2], p2[0], delta_a);
	/*
	// rotate p2 with z
	r = sqrt(target.x * target.x + target.y * target.y);
	if (r != 0) {
		double c = target.y / r;
		double s = target.x / r;
		for (int i = 0; i < 7; i++) {
			double x_ = c * p2[i].x - s * p2[i].y;
			double y_ = s * p2[i].y + c * p2[i].y;
			double z_ = p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}
	// rotate p2 with x
	r = sqrt(target.x * target.x + target.y * target.y + target.z * target.z);
	if (r != 0) {
		double c = target.z / r;
		double s = sqrt(target.x * target.x + target.y * target.y) / r;
		for (int i = 0; i < 7; i++) {
			double x_ = p2[i].x;
			double y_ = c * p2[i].y - s * p2[i].z;
			double z_ = s * p2[i].y + c * p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}
	// rotate p2 with z
	double c = cos(delta_a / 180 * 3.1415927);
	double s = sin(delta_a / 180 * 3.1425927);
	for (int i = 0; i < 7; i++) {
			double x_ = c * p2[i].x - s * p2[i].y;
			double y_ = s * p2[i].y + c * p2[i].y;
			double z_ = p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
	}
	// rotate p2 with x
	r = sqrt(target.x * target.x + target.y * target.y + target.z * target.z);
	if (r != 0) {
		c = target.z / r;
		s = -sqrt(target.x * target.x + target.y * target.y) / r;
		for (int i = 0; i < 7; i++) {
			double x_ = p2[i].x;
			double y_ = c * p2[i].y - s * p2[i].z;
			double z_ = s * p2[i].y + c * p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}
	// rotate p2 with z
	r = sqrt(target.x * target.x + target.y * target.y);
	if (r != 0) {
		double c = target.y / r;
		double s = -target.x / r;
		for (int i = 0; i < 7; i++) {
			double x_ = c * p2[i].x - s * p2[i].y;
			double y_ = s * p2[i].y + c * p2[i].y;
			double z_ = p2[i].z;
			p2[i].x = x_;
			p2[i].y = y_;
			p2[i].z = z_;
		}
	}
	*/

	// translate p2
	for (int i = 0; i < 7; i++) {
		p2[i].x += p1[4].x;
		p2[i].y += p1[4].y;
		p2[i].z += p1[4].z;
	}

	Point *p, *p3;
	int len3;
	if (type == 0) {
		len3 = 11;
		p = new Point[22];
		p3 = new Point[11];
		f_name = lib + "A";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 11; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 1) {
		len3 = 9;
		p = new Point[20];
		p3 = new Point[9];
		f_name = lib + "U";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 9; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 2) {
		len3 = 12;
		p = new Point[23];
		p3 = new Point[12];
		f_name = lib + "G";
		ifile.open(f_name.c_str());
		for (int i = 0; i < 12; i++) {
			ifile >> p3[i].x >> p3[i].y >> p3[i].z;
		}
		ifile.close();
	} else if (type == 3) {
		len3 = 9;
		p = new Point[20];
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

	/*
	// rotate p3 with z
	r = sqrt(p3[1].x * p3[1].x + p3[1].y * p3[1].y);
	if (r != 0) {
		c = p3[1].y / r;
		s = p3[1].x / r;
		for (int i = 0; i < len3; i++) {
			double x_ = c * p3[i].x - s * p3[i].y;
			double y_ = s * p3[i].x + c * p3[i].y;
			double z_ = p3[i].z;
			p3[i].x = x_;
			p3[i].y = y_;
			p3[i].z = z_;
		}
	}

	// rotate p3 with x
	a1 = acos(target.z / sqrt(target.x * target.x + target.y * target.y + target.z * target.z));
	r = sqrt(p3[1].y * p3[1].y + p3[1].z * p3[1].z);
	double a2 = acos(p3[1].z / r);
	delta_a = a2 - a1;
	c = cos(delta_a);
	s = sin(delta_a);
	for (int i = 0; i < len3; i++) {
		double x_ = p3[i].x;
		double y_ = c * p3[i].y - s * p3[i].z;
		double z_ = s * p3[i].y + c * p3[i].z;
		p3[i].x = x_;
		p3[i].y = y_;
		p3[i].z = z_;
	}

	// rotate p3 with z
	r = sqrt(target.x * target.x + target.y * target.y);
	if (r != 0) {
		c = target.y / r;
		s = -target.x / r;
		for (int i = 0; i < len3; i++) {
			double x_ = c * p3[i].x - s * p3[i].y;
			double y_ = s * p3[i].x + c * p3[i].y;
			double z_ = p3[i].z;
			p3[i].x = x_;
			p3[i].y = y_;
			p3[i].z = z_;
		}
	}
	*/

	// normalVector
	n1 = Point::normalVector(p2[1], p2[5], p2[3]);

	// dihedral
	double dih = Point::dihedral(n1, &(p3[0]), &(p3[1]), &(p3[2]));
	double delta_dih = chi - dih;

	// rotate
	Point::rotate(p3, len3, *origin, p3[1], delta_dih);
	/*
	Point spindle(p3[1].x, p3[1].y, p3[1].z);
	r = sqrt(p3[1].x * p3[1].x + p3[1].y * p3[1].y);
	if (r != 0) {
		c = p3[1].y / r;
		s = p3[1].x / r;
		for (int i = 0; i < len3; i++) {
			double x_ = c * p3[i].x - s * p3[i].y;
			double y_ = s * p3[i].x + c * p3[i].y;
			double z_ = p3[i].z;
			p3[i].x = x_;
			p3[i].y = y_;
			p3[i].z = z_;
		}
	}

	r = sqrt(p3[1].y * p3[1].y + p3[1].z * p3[1].z);
	if (r != 0) {
		c = p3[1].z / r;
		s = p3[1].y / r;
		for (int i = 0; i < len3; i++) {
			double x_ = p3[i].x;
			double y_ = c * p3[i].y - s * p3[i].z;
			double z_ = s * p3[i].y + c * p3[i].z;
			p3[i].x = x_;
			p3[i].y = y_;
			p3[i].z = z_;
		}
	}
	
	c = cos(delta_dih / 180. * 3.1415927);
	s = sin(delta_dih / 180. * 3.1415927);
	for (int i = 0; i < len3; i++) {
		double x_ = c * p3[i].x - s * p3[i].y;
		double y_ = s * p3[i].x + c * p3[i].y;
		double z_ = p3[i].z;
		p3[i].x = x_;
		p3[i].y = y_;
		p3[i].z = z_;
	}

	r = sqrt(spindle.x * spindle.x + spindle.y * spindle.y + spindle.z * spindle.z);
	if (r != 0) {
		c = spindle.z / r;
		s = -sqrt(spindle.x * spindle.x + spindle.y * spindle.y) / r;
		for (int i = 0; i < len3; i++) {
			double x_ = p3[i].x;
			double y_ = c * p3[i].y - s * p3[i].z;
			double z_ = s * p3[i].y + c * p3[i].z;
			p3[i].x = x_;
			p3[i].y = y_;
			p3[i].z = z_;
		}
	}

	r = sqrt(spindle.x * spindle.x + spindle.y * spindle.y);
	if (r != 0) {
		c = spindle.y / r;
		s = -spindle.x / r;
		for (int i = 0; i < len3; i++) {
			double x_ = c * p3[i].x - s * p3[i].y;
			double y_ = s * p3[i].x + c * p3[i].y;
			double z_ = p3[i].z;
			p3[i].x = x_;
			p3[i].y = y_;
			p3[i].z = z_;
		}
	}
	*/

	for (int i = 0; i < len3; i++) {
		p3[i].x += p2[5].x;
		p3[i].y += p2[5].y;
		p3[i].z += p2[5].z;
	}

	Point o3(o[0].x - o[1].x, o[0].y - o[1].y, o[0].z - o[1].z);
	Point o5(o[2].x - o[1].x, o[2].y - o[1].y, o[2].z - o[1].z);
	double ratio = 1.5 / sqrt(o3.x * o3.x + o3.y * o3.y + o3.z * o3.z);
	Point o1p(ratio * o3.x, ratio * o3.y, ratio * o3.z);
	Point o2p(ratio * o3.x, ratio * o3.y, ratio * o3.z);
	Point *spindle1 = Point::normalVector(o5, *origin, o3);
	Point::rotate(&o1p, 1, *origin, *spindle1, 109);
	Point::rotate(&o2p, 1, *origin, *spindle1, 109);
	Point::rotate(&o1p, 1, *origin, o3, 120);
	Point::rotate(&o2p, 1, *origin, o3, 240);

	p[0].x = p1[0].x; p[0].y = p1[0].y; p[0].z = p1[0].z;
	p[1].x = o1p.x + o[1].x; p[1].y = o1p.y + o[1].y; p[1].z = o1p.z + o[1].z;
	p[2].x = o2p.x + o[1].x; p[2].y = o2p.y + o[1].y; p[2].z = o2p.z + o[1].z;
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

	Residue *residue = new Residue(p, len3 + 11, type);
	return residue;
}




















