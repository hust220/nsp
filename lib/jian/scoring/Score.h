#ifndef SCORE_H
#define SCORE_H

#include "DistAnal.h"
#include "DihAnal.h"

namespace jian {
namespace scoring {

class Score {
public:
    DistAnal * _dist_anal;
    DihAnal * _dih_anal;

    std::string _par_dist_obs;
    std::string _par_dist_ref;
    std::string _par_dih;
    std::string _ref_state = "average";
    std::string _lib;
    std::string _file_list;
    std::string _input;

    double _dist_bin = 0.3;
    double _dih_bin = 4.5;
    int _cutoff = 20;

    double _constant = 27.1118;
    double _dist_weight = 0.433513;
    double _dih_weight = 1.59348;

    std::string _type = "RNA";

    Score(Par par) {
        std::string _lib = env("NSP");

        if (par["global"].size() == 2) _input = par["global"][1];
        else if (par.count("lib")) _lib = par["lib"][0];
        else if (par.count("list")) _file_list = par["list"][0];
        else if (par.count("type")) _type = par["type"][0];
        else if (par.count("dist")) _par_dist_obs = par["dist"][0];
        else if (par.count("dih")) _par_dih = par["dih"][0];
        else if (par.count("cutoff")) _cutoff = std::stoi(par["cutoff"][0]);

        _par_dist_obs = _lib + "/" + _type + "/dist.par";
        _par_dist_ref = "";
        _par_dih = _lib + "/" + _type + "/dih.par";

        _dist_anal = new DistAnal(_cutoff, _dist_bin, _ref_state);
        _dist_anal->read_obs_parm(_par_dist_obs);
        _dist_anal->read_ref_parm(_par_dist_ref);

        _dih_anal = new DihAnal(_dih_bin);
        _dih_anal->readParm(_par_dih);
    }

    ~Score() {
        delete _dist_anal;
        delete _dih_anal;
    }

    void operator ()() {
        if (_file_list != "") {
            std::ifstream ifile(_file_list.c_str());
            std::string line;
            while (ifile >> line) {
                (*_dist_anal)(Model(line));
                std::cout << RowVectorXf(_dist_anal->_scores) << std::endl;
            }
            ifile.close();
        } else {
            (*_dist_anal)(Model(_input));
            std::cout << RowVectorXf(_dist_anal->_scores) << std::endl;
        }
    }

};

} /// namespace scoring
} /// namespace jian

#endif

