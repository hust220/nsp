#ifndef MOL2D_H
#define MOL2D_H

#include "loop.h"

class Mol2D {
public:
	//Mol2D(string ss_, int view_ = 0) : line(del_single_pair(ss_)), view(view_) {
	Mol2D(string ss_, int view_ = 0) : line(ss_), view(view_) {
		// set ss 
		stringstream sstr;
		std::set<char> temp_set{'.', '(', ')', '[', ']', '{', '}'};
		for (int i = 0; i < line.size(); i++) {
			if (temp_set.count(line[i])) {
				sstr << line[i];
			} else {
				if (line[i] != '&') {
					cerr << "Mol2D::Mol2D() error! RNA secondary structure should only includes "
					     << "'.', '(', ')', '[', ']', '{', '}' and '&', but not includes '" 
						 << line[i] << "'." << endl;
					cerr << "But This secondary structure is: " << endl;
					cerr << line << endl;
					exit(1);
				}
			}
		}
		sstr >> ss;

		// set 2d structure tree
		res *r = new res('^', 0);
		vector<res> v;
		v.push_back(*r);
		int i = 0;
		for (i = 0; i < (int) line.size(); i++) {
			v.push_back(res(line[i], i + 1));
		}
		r = new res('&', i + 1);
		v.push_back(*r);

		// set tree
		setTree(v, view);
	}

	Mol2D(string ss_, string seq_, int view_ = 0) : Mol2D(ss_, view_) {
		readSeq(seq_);
	}

	Mol2D(Mol2D *mol2d) : line(mol2d->line), seq(mol2d->seq), ss(mol2d->ss), 
	                            head(loop::copy(mol2d->head)), pseudo_head(loop::copy(mol2d->pseudo_head)), 
								mol(new RNA(mol2d->mol)), view(mol2d->view) {}
	Mol2D(const Mol2D &mol2d) : line(mol2d.line), seq(mol2d.seq), ss(mol2d.ss), 
	                            head(loop::copy(mol2d.head)), pseudo_head(loop::copy(mol2d.pseudo_head)), 
								mol(new RNA(mol2d.mol)), view(mol2d.view) {}
	Mol2D &operator =(const Mol2D &mol2d) {
		line = mol2d.line;
		seq = mol2d.seq;
		ss = mol2d.ss;
		loop::del(head);
		loop::del(pseudo_head);
		head = loop::copy(mol2d.head);
		pseudo_head = loop::copy(mol2d.pseudo_head);
		delete mol;
		mol = new RNA(mol2d.mol);
		view = mol2d.view;
	}
	~Mol2D() {
		delete mol;
		loop::del(head);
		loop::del(pseudo_head);
	}

	string del_single_pair(string);
	void setTree(vector<res> &, int = 0);
	void resetNum(loop *, const vector<int> &);
	void print();
	void readSeq(string);
	void readMol(string);

	void getLoop(vector<res> &, vector<loop *> &, const vector<int> &);
	void printTree(loop *, int = 0);
	void delLoop(loop *);
	void setSeq(loop *, string);

	pair<vector<loop *>, int> extend_hinge(loop *, int);

	// methods for constructing pseudo-knots loop
	void setPairs();
	void setLoops(loop *);
	list<loop *> find_path(loop *, loop *);
	void constructPseudoknots();
	list<loop *> pseudo_loop(loop *, set<loop *>);
	int pseudo_tree(loop *, loop *, set<loop *>, loop *);

	int hinge_base_pair_num = 2;
	vector<loop *> loops;
	vector<int> pairs;
	loop *head = NULL;
	loop *pseudo_head = NULL;
	RNA *mol = NULL;
	string line;
	string seq;
	string ss;
	int view = 0;
};

#endif
