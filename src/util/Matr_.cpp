#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "Matr_.h"

namespace jian {

using namespace std;

MatrixXf mat_from_file(std::string file) {
    std::ifstream ifile(file.c_str());
    int rows, cols;
    ifile >> rows >> cols;
    MatrixXf mat(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (ifile >> mat(i, j)) {
                continue;
            } else {
                std::cerr << "jian::mat_from_file failed!" << std::endl;
                exit(1);
            }
        }
    }
    return mat;
}


Matr_::Matr_(int row, int col) {
	if (row <= 0 || col <= 0) {
		cout << "Matr_ construct error!" << endl;
		exit(1);
	}

	this->row = row;
	this->col = col;

	data = new double *[row];
	for (int i = 0; i < row; i++) {
		data[i] = new double[col];
	}
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (i == j) {
				data[i][j] = 1;
			} else {
				data[i][j] = 0;
			}
		}
	}
}

Matr_::Matr_(int row, int col, int n) {
	if (row <= 0 || col <= 0) {
		cout << "Matr_ construct error!" << endl;
		exit(1);
	}

	this->row = row;
	this->col = col;

	data = new double *[row];
	for (int i = 0; i < row; i++) {
		data[i] = new double[col];
	}
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			data[i][j] = n;
		}
	}
}

Matr_::Matr_(Matr_ *m) {
	int i, j;
	
	row = m->row;
	col = m->col;
	
	data = new double *[row];
	for (i = 0; i < row; i++) {
		data[i] = new double[col];
		for (j = 0; j < col; j++) {
			data[i][j] = m->data[i][j];
		}
	}
}

Matr_::Matr_(Obj<Matr_> m) {
	int i, j;
	
	row = m->row;
	col = m->col;
	
	data = new double *[row];
	for (i = 0; i < row; i++) {
		data[i] = new double[col];
		for (j = 0; j < col; j++) {
			data[i][j] = m->data[i][j];
		}
	}
}

Matr_::Matr_(string str) {
	row = 1;
	col = 0;
	for (int i = 0; i < str.size(); i++) {
		if (str[i] == ';') {
			row++;
		}
	}

	vector<string> line;
	tokenize(str, line, " ;");
	col = line.size() / row;

	data = new double *[row];
	for (int i = 0; i < row; i++) {
		data[i] = new double[col];
		for (int j = 0; j < col; j++) {
			data[i][j] = atof(line[i * col + j].c_str());
		}
	}
}

Matr_::~Matr_() {
	for (int i = 0; i < row; i++) {
		delete [] data[i];
	}
	delete [] data;
}

Matr_ *Matr_::copy() {
	Matr_ *temp = new Matr_(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			temp->data[i][j] = data[i][j];
		}
	}
	return temp;
}

double *Matr_::operator [](int n) {
	return data[n];
}

void Matr_::read(string name) {
	int i, j;
	ifstream ifile;

	ifile.open(name.c_str());
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			ifile >> data[i][j];
		}
	}
	ifile.close();
}

void Matr_::print(int n) {
	int i, j;
	
	cout << fixed << setprecision(3);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			cout << setw(n) << data[i][j];
		}
		cout << endl;
	}
	cout << endl;
}

void Matr_::identity() {
	int i, j;
	
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			if (i == j) {
				data[i][j] = 1;
			} else {
				data[i][j] = 0;
			}
		}
	}
}

void Matr_::changeRow(int row1, int row2) {
	int j;
	double *row3;
	
	row3 = new double[col];
	for (j = 0; j < col; j++) {
		row3[j] = data[row1][j];
		data[row1][j] = data[row2][j];
		data[row2][j] = row3[j];
	}
	delete row3;
}

Matr_ *Matr_::colVec(int j) {
	Matr_ *a = new Matr_(row, 1);
	for (int i = 0; i < row; i++) {
		a->data[i][0] = data[i][j];
	}
	return a;
}

Matr_ *Matr_::subMatr_(int x, int y) {
	if (x >= row || y >= col || x < 0 || y < 0) {
		cout << "Matr_::subMatr_ error!" << endl;
		exit(1);
	}

	Matr_ *m = new Matr_(row - 1, col - 1);
	int i, j, mi, mj;
	for (i = 0, mi = 0; i < row; i++) {
		if (i == x) continue;
		for (j = 0, mj = 0; j < col; j++) {
			if (j == y) continue;
			m->data[mi][mj] = data[i][j];
			mj++;
		}
		mi++;
	}
	return m;
}

double Matr_::det() {
	if (row != col) {
		cout << "Matr_::det error!" << endl;
		exit(1);
	}

	if (row == 0) {
		cout << "Matr_::det error!" << endl;
		exit(2);
	} else if (row == 1) {
		return data[0][0];
	} else if (row == 2) {
		return data[0][0] * data[1][1] - data[0][1] * data[1][0];
	} else if (row == 3) {
		return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) + data[0][1] * (data[1][2] * data[2][0] - data[1][0] * data[2][2]) + data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
	}

	int i, a = 1;
	double sum = 0;
	for (i = 0; i < col; i++) {
		sum += data[0][i] * a * subMatr_(0, i)->det();
		a = -a;
	}
	return sum;
}

void Matr_::assign(int row, int col) {
	if (row <= 0 || col <= 0) {
		cout << "Matr_::assign error!" << endl;
		exit(1);
	}

	for (int i = 0; i < this->row; i++) {
		delete [] data[i];
	}
	delete [] data;

	this->row = row;
	this->col = col;

	data = new double *[row];
	for (int i = 0; i < row; i++) {
		data[i] = new double[col];
	}
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (i == j) {
				data[i][j] = 1;
			} else {
				data[i][j] = 0;
			}
		}
	}
}

void Matr_::assign(Matr_ *a) {
	for (int i = 0; i < row; i++) {
		delete [] data[i];
	}
	delete [] data;
	data = a->data;
	row = a->row;
	col = a->col;
}

Matr_ *Matr_::add(double d) {
	Matr_ *a = new Matr_(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (i == j) {
				a->data[i][j] = d;
			}
		}
	}
	Matr_ *b = add(a);
	delete a;
	return b;
}

Matr_ *Matr_::add(Matr_ *a) {
	if (row != a->row || col != a->col) {
		cout << "Matr_::add error!" << endl;
		exit(1);
	}
	Matr_ *b = new Matr_(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			b->data[i][j] = data[i][j] + a->data[i][j];
		}
	}
	return b;
}

Matr_ *Matr_::minus(double d) {
	return add(-d);
}

Matr_ *Matr_::minus(Matr_ *a) {
	if (row != a->row || col != a->col) {
		cout << "Matr_::minus error!" << endl;
		exit(1);
	}
	Matr_ *b = new Matr_(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			b->data[i][j] = data[i][j] - a->data[i][j];
		}
	}
	return b;
}

Matr_ *Matr_::multiply(double d) {
	Matr_ *b = new Matr_(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			b->data[i][j] = data[i][j] * d;
		}
	}
	return b;
}

Matr_ *Matr_::multiply(Matr_ *m) {
	int i, j, k;
	Matr_ *n;
	
	if (col != m->row) {
		cout << "matrix multiply error" << endl;
		exit(1);
	}
	n = new Matr_(row, m->col);
	for (i = 0; i < row; i++) {
		for (j = 0; j < m->col; j++) {
			n->data[i][j] = 0;
			for (k = 0; k < col; k++) {
				n->data[i][j] += data[i][k] * m->data[k][j];
			}
		}
	}
	return n;
}

Obj<Matr_> Matr_::dot(Obj<Matr_> m) {
	int i, j, k;
	Obj<Matr_> n;
	
	if (col != m->row) {
		cout << "matrix multiply error" << endl;
		exit(1);
	}
	n = new Matr_(row, m->col);
	for (i = 0; i < row; i++) {
		for (j = 0; j < m->col; j++) {
			n->data[i][j] = 0;
			for (k = 0; k < col; k++) {
				n->data[i][j] += data[i][k] * m->data[k][j];
			}
		}
	}
	return n;
}

Matr_ *Matr_::sqrt() {
	/* create a new Matr_ m */
	Matr_ *m = new Matr_(row, col);

	/* initialize matrix m */
	int i, j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			if (data[i][j] == 0) {
				m->data[i][j] = 0;
			} else if (data[i][j] < 0) {
				cout << "Matr_::sqrt error!" << endl; // the data must be greater than zero
				exit(1);
			} else {
				m->data[i][j] = std::sqrt(data[i][j]);
			}
		}
	}

	/* return matrix m */
	return m;
}

Matr_* Matr_::inverse() {
	int i, j, j2, temp;
	double flag;
	Matr_ *l, *r;
	
	if (row != col) {
		cout << "Matr_::inverse error!" << endl;
		exit(1);
	}

	if (det() == 0) {
		cout << "Matr_::inverse error!" << endl;
		cout << "   determinant is zero!" << endl;
	}

	l = new Matr_(this);
	r = new Matr_(row, row);
	r->identity();

	for (j = 0; j < col; j++) {
		if (l->data[j][j] == 0) {
			for (temp = j + 1; temp < row; temp++) {
				if (l->data[temp][j] != 0) {
					l->changeRow(j, temp);
					r->changeRow(j, temp);
					break;
				}
			}
			if (temp == row) {
				cout << "Matr_::inverse error!" << endl;
				exit(1);
			}
		}
		flag = l->data[j][j];
		if (flag == 0) cout << "division error: Matr_.cpp 99\n";
		for (j2 = 0; j2 < col; j2++) {
			l->data[j][j2] /= flag;
			r->data[j][j2] /= flag;
		}

		for (i = j + 1; i < row; i++) {
			if (l->data[i][j] == 0) {
				for (temp = i + 1; temp < row; temp++) {
					if (l->data[temp][j] != 0) {
						l->changeRow(i, temp);
						r->changeRow(i, temp);
						break;
					}
				}
				if (temp == row) {
					cout << "matrix inverse error!" << endl;
					exit(1);
				}
			}
			flag = l->data[i][j];
		if (flag == 0) cout << "division error: Matr_.cpp 120\n";
			for (j2 = 0; j2 < col; j2++) {
				l->data[i][j2] /= flag;
				r->data[i][j2] /= flag;
				l->data[i][j2] -= l->data[j][j2];
				r->data[i][j2] -= r->data[j][j2];
			}
		}
	}

	for (j = col - 1; j > 0; j--) {
		for (i = j - 1; i >= 0; i--) {
			flag = l->data[i][j];
			for (j2 = col - 1; j2 >= 0; j2--) {
				l->data[i][j2] -= flag * l->data[j][j2];
				r->data[i][j2] -= flag * r->data[j][j2];
			}
		}
	}

	delete l;
	return r;
}

Matr_* Matr_::transpose() {
	Matr_ *m;
	int i, j;

	m = new Matr_(col, row);
	for (i = 0; i < col; i++) {
		for (j = 0; j < row; j++) {
			m->data[i][j] = data[j][i];
		}
	}
	return m;
}


/* qr decomposition using Householder method */
void Matr_::qr(Matr_ *q, Matr_ *r) {
	if (row != col) {
		cout << "Matr_::qr error!" << endl;
		cout << "   This Matr_ is not square!" << endl;
		exit(1);
	}
	if (det() == 0) {
		cout << "Matr_::qr error!" << endl;
		cout << "   This Matr_ is irreversible!" << endl;
		exit(1);
	}

	q->assign(row, col);
	r->assign(this);

	for (int i = 0; i < col - 1; i++) {
		Matr_ *H = new Matr_(row, col);
		Matr_ *x = r->colVec(i);
		Matr_ *x_ = x->transpose();
		x->assign(x_);
		Matr_ *y = new Matr_(1, col);
		for (int j = 0; j < col; j++) {
			if (j > i) {
				y->data[0][j] = 0;
			} else if (j == i) {
				double sum = 0;
				for (int k = 0; k < col; k++) {
					sum += x->data[0][k] * x->data[0][k];
				}
				for (int k = 0; k < i; k++) {
					sum -= x->data[0][k] * x->data[0][k];
				}
				if (x->data[0][j] > 0) {
					y->data[0][j] = -std::sqrt(sum);
				} else {
					y->data[0][j] = std::sqrt(sum);
				}
			} else {
				y->data[0][j] = x->data[0][j];
			}
		}
		Matr_ *t = x->minus(y);
		Matr_ *t_ = t->transpose();
		Matr_ *m = t_->multiply(t);
		double sum = 0;
		for (int j = 0; j < col; j++) {
			sum += t->data[0][j] * t->data[0][j];
		}
		Matr_ *m_ = m->multiply(2. / sum);
		// m->assign(m->multiply(2. / sum));
		Matr_ *h_ = H->minus(m_);
		// H->assign(H->minus(m));
		Matr_ *q_ = q->multiply(h_);
		q->assign(q_);
		// q->assign(q->multiply(H));
		Matr_ *r_ = h_->multiply(r);
		r->assign(r_);
		// r->assign(H->multiply(r));
	}
}

void Matr_::ed(Matr_ *p1, Matr_ *e, Matr_ *p2) {
	if (row != col) {
		cout << "Matr_::ed error!" << endl;
		exit(1);
	}
	Matr_ *q = new Matr_(row, col);
	Matr_ *r = new Matr_(row, col);
	e->assign(this);
	int n = 0;
	while (1) {//cout << n++ << endl; p1->print(); e->print();
		if (n >= 10) break;
		e->qr(q, r);
		//cout << "-----------------" << endl;
		e->assign(r->multiply(q));
		p1->assign(p1->multiply(q));
		int flag = 0;
		int flag2 = 0;
		for (int i = 1; i < row; i++) {
			int j;
			for (j = 0; j < i - 1; j++) {
				if (abs(e->data[i][j]) > 0.000001) {
					flag = 1;
				}
			}
			if (abs(e->data[i][j]) > 0.000001) {
				flag2++;
			} else {
				flag2 = 0;
			}
			if (flag2 == 2) {
				flag = 1;
				break;
			}
		}
		if (flag == 0) {
			break;
		}
	}
	p2->assign(p1->inverse());
}

void Matr_::svd(Matr_ *q1, Matr_ *s, Matr_ *q2) {
	/* check the matrix */
	if (row != col) {
		cout << "Matr_::svd error!" << endl;
		exit(1);
	}
	
	Matr_ *t = transpose();
	Matr_ *ts = t->multiply(this);
	Matr_ *p1, *e, *p2;
	ts->ed(p1, e, p2);
	s = e->sqrt();
	Matr_ *s_ = s->inverse();
	q2 = p1;
	q1 = multiply(q2)->multiply(s_);
}

Obj<Point> Matr_::point(int n) {
	Obj<Point> p(new Point);
	if (n >= 0 && n < row && col > 0) {
		for (int i = 0; i < col; i++) {
			(*p)[i] = data[n][i];
		}
	} else {
		cerr << "Matr_::point error!" << endl;
		exit(1);
	}
	return p;
}

} /// namespace jian


