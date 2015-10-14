#include "Train.h"

namespace jian {
namespace scoring {

Train::Train(Par par) {
    if (par.count("cutoff")) _cutoff = std::stoi(par["cutoff"][0]);
    else if (par.count("bin")) _bin = std::stof(par["bin"][0]);
    else if (par.count("dist")) _par_dist = par["dist"][0];
    else if (par.count("dih")) _par_dist = par["dist"][0];
    else if (par.count("list")) _file_list = par["list"][0];
    _dist_anal = new DistAnal(_cutoff, _bin);
    if (_par_dist != "") _dist_anal->read_obs_parm(_par_dist);
}

Train::~Train() {
    delete _dist_anal;
}

void Train::operator ()() {
    std::string str;
    std::ifstream ifile(_file_list.c_str());
    int n = 0;
    while (ifile >> str) {
        n++;
        cerr << n << ". train: " << str << endl;
        _dist_anal->read_mol(RNA(str));
        _dist_anal->train();
    }
    ifile.close();
    print_par(_dist_anal->_obs_parm);
}

//void Train::dih(char *filename) {
//    string str;
//    DihAnal dihAnal;
//    if (par_dih != "") {
//        dihAnal.readParm(par_dih);
//    }
//    ifstream ifile(filename);
//    int n = 0;
//    while (ifile >> str) {
//        n++;
//        cerr << n << ". train: " << str << endl;
//        dihAnal.readRNA(RNA(str));
//        dihAnal.train();
//    }
//    ifile.close();
//    dihAnal.printParm();
//}
//
void Train::print_par(const MatrixXd &mat) {
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols() / _bins; j++) {
            for (int k = 0; k < _bins; k++) {
                std::cout << mat(i, j * _bins + k) << ' ';
            }
            std::cout << std::endl;
        }
    }
}

} /// namespace scoring
} /// namespace jian




