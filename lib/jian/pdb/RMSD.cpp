#include "RMSD.hpp"

namespace jian {

double RMSD::operator ()(const Model &model1, const Model &model2) {
    return run(model1, model2);
}

void RMSD::setXY() {
    int len1 = num_residues(_model1);
    int len2 = num_residues(_model2);
    assert(len1 == len2);

    std::map<std::string, Atom> atom_map;
    std::vector<Atom> atom_list1;
    std::vector<Atom> atom_list2;

    auto chain1 = _model1.begin();
    auto chain2 = _model2.begin();
    auto res1 = chain1->begin();
    auto res2 = chain2->begin();

    int n = 0;
    while (true) {
        n++;

        assert(res1->name == res2->name);
        for (auto &&atom1 : *res1) {
            atom_map[atom1.name] = atom1;
        }
        for (auto &&atom2 : *res2) {
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

    Eigen::Matrix3f g = x * y.transpose();

    Eigen::JacobiSVD<Eigen::Matrix3f> svd(g, Eigen::ComputeFullU|Eigen::ComputeFullV);
    Eigen::Matrix3f u = svd.matrixU();
    Eigen::Matrix3f v = svd.matrixV();

    double det = g.determinant();
    Eigen::Matrix3f I;
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
    Eigen::Matrix3f r = v * I * u.transpose();

    Eigen::MatrixXf x_ = r * x;

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

} // namespace jian

