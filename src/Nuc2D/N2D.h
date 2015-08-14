#ifndef N2D_H
#define N2D_H

#include "loop.h"
#include "../Pdb.h"

namespace jian {

namespace nuc2d {

class N2D {
public:
    N2D();
    N2D(string ss_, int view_ = 0);
    N2D(string ss_, string seq_, int view_ = 0);
    N2D(N2D *mol2d);
    N2D(const N2D &mol2d);
    void operator ()();
    void operator ()(std::string ss_);
    void operator ()(std::string seq_, std::string ss_);
    N2D &operator =(const N2D &mol2d);
    ~N2D();

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
    Model mol;
    string line;
    string seq;
    string ss;
    int view = 0;
};

} // namespace nuc2d

} // namespace jian

#endif

