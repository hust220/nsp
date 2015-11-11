#include "RMSD.h"

namespace jian {

namespace pdb {

void RMSD::setXY() {
    int len1 = _model1.res_nums();
    int len2 = _model2.res_nums();
    assert(len1 == len2);

    std::map<std::string, Atom> atom_map;
    std::vector<Atom> atom_list1;
    std::vector<Atom> atom_list2;

    auto chain1 = _model1.begin();
    auto chain2 = _model2.begin();
    auto res1 = chain1->begin();
    auto res2 = chain2->begin();

    int n = 0;
    while (1) {
        n++;

        assert(res1->name == res2->name);
        for (auto &&atom1 : res1->atoms) {
            atom_map[atom1.name] = atom1;
        }
        for (auto &&atom2 : res2->atoms) {
            if (atom_map.count(atom2.name) == 1) {
                atom_list1.push_back(atom_map[atom2.name]);
                atom_list2.push_back(atom2);
            }
        }
        atom_map.clear();

        res1++;
        if (res1 == chain1->end()) {
            chain1++;
            if (chain1 == _model1.end()) {
                break;
            }
            res1 = chain1->begin();
        }
        res2++;
        if (res2 == chain2->end()) {
            chain2++;
            if (chain2 == _model2.end()) {
                break;
            }
            res2 = chain2->begin();
        }
    }

    _len = atom_list1.size();

    x.resize(3, _len);
    y.resize(3, _len);

    /* set X */
    double c1[3] = {0, 0, 0};
    double c2[3] = {0, 0, 0};
    for (int i = 0; i < _len; i++) {
        for (int j = 0; j < 3; j++) {
            x(j, i) = atom_list1[i][j];
            c1[j] += x(j, i);
            y(j, i) = atom_list2[i][j];
            c2[j] += y(j, i);
        }
//        x(0, i) = atom_list1[i]._x;
//        c1[0] += x(0, i);
//        x(1, i) = atom_list1[i]._y;
//        c1[1] += x(1, i);
//        x(2, i) = atom_list1[i]._z;
//        c1[2] += x(2, i);
//        y(0, i) = atom_list2[i]._x;
//        c2[0] += y(0, i);
//        y(1, i) = atom_list2[i]._y;
//        c2[1] += y(1, i);
//        y(2, i) = atom_list2[i]._z;
//        c2[2] += y(2, i);
    }
    for (int i = 0; i < 3; i++) {
        c1[i] /= _len;
        c2[i] /= _len;
    }

    /* translate the center of X and Y to origin */
    for (int i = 0; i < _len; i++) {
        for (int j = 0; j < 3; j++) {
            x(j, i) = x(j, i) - c1[j];
            y(j, i) = y(j, i) - c2[j];
        }
    }
}

double RMSD::run(const Model &model1, const Model &model2) {
    _model1 = model1;
    _model2 = model2;
    setXY();

    Matrix3f g = x * y.transpose();

    JacobiSVD<Matrix3f> svd(g, ComputeFullU|ComputeFullV);
    Matrix3f u = svd.matrixU();
    Matrix3f v = svd.matrixV();

    double det = g.determinant();
    Matrix3f I;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                I(i, j) = 1;
            } else {
                I(i, j) = 0;
            }
        }
    }
    if (det < 0) {
        I(2, 2) = -1;
    }
    Matrix3f r = v * I * u.transpose();

    MatrixXf x_ = r * x;

    rmsd = 0;
    for (int i = 0; i < _len; i++) {
        for (int j = 0; j < 3; j++) {
            rmsd += (x_(j, i) - y(j, i)) * (x_(j, i) - y(j, i));
        }
    }
    rmsd /= _len;
    rmsd = sqrt(rmsd);
    return rmsd;
}

} /// namespace pdb

} /// namespace jian







