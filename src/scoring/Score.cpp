#include "Score.h"

namespace jian {
namespace scoring {

Score::Score(Par par) {
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

Score::~Score() {
    delete _dist_anal;
    delete _dih_anal;
}

//Score::Score() {
//    std::string lib = env("NSP");
//    par_dist_obs = lib + "/" + "_type" + "/dist.par";
//    par_dist_ref = "";
//    dihPar = lib + "/" + "_type" + "/dih.par";
//
//    distAnal = new DistAnal(cutoff, _dist_bin, reference_state);
//    dihAnal = new DihAnal(_dih_bin);
//
//    distAnal->read_obs_parm(par_dist_obs);
//    distAnal->read_ref_parm(par_dist_ref);
//    dihAnal->readParm(dihPar);
//}
//
//Score::Score(char *parm) {
//    std::string lib = env("NSP");
//    par_dist_obs = lib + "/" + "_type" + "/dist.par";
//    par_dist_ref = "";
//    dihPar = lib + "/" + "_type" + "/dih.par";
//
//    std::ifstream ifile(parm);
//    while (ifile) {
//        string line;
//        getline(ifile, line);
//        vector<string> splited_line;
//        tokenize(line, splited_line, " ,:");
//        if (splited_line.size() != 2) continue;
//        if (splited_line[0] == "par_dist_obs") {
//            par_dist_obs = splited_line[1];
//        } else if (splited_line[0] == "par_dist_ref") {
//            par_dist_ref = splited_line[1]; 
//        } else if (splited_line[0] == "par_dih") {
//            dihPar = splited_line[1]; 
//        } else if (splited_line[0] == "cutoff") {
//            cutoff = atoi(splited_line[1].c_str());
//        } else if (splited_line[0] == "dist_weight") {
//            _distWeight = atof(splited_line[1].c_str());
//        } else if (splited_line[0] == "dih_weight") {
//            _dihWeight = atof(splited_line[1].c_str());
//        } else if (splited_line[0] == "reference_state") {
//            reference_state = splited_line[1];
//        }
//    }
//    ifile.close();
//
//    distAnal = new DistAnal(cutoff, _dist_bin, reference_state);
//    distAnal->read_obs_parm(par_dist_obs);
//    distAnal->read_ref_parm(par_dist_ref);
//
//    dihAnal = new DihAnal(_dih_bin);
//    dihAnal->readParm(dihPar);
//}

void Score::operator ()() {
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

//    double dist_score = (*_dist_anal)(model);
//    std::cout << _dist_anal->_scores << std::endl;
//    return dist_score;
//    dihAnal->readRNA(model);
//    double dihScore = dihAnal->scoring();
    // cerr << model->name << ' ' << distScore << ' ' << dihScore << endl;
    // double *score = distAnal->getScore();
    // cerr << model->name << ' ' << score[0] << ' ' << score[1] << ' ' << score[2] << ' ' << score[3] << ' ' << score[4] << ' ' << dihScore << endl;
    //std::cout << distScore << ' ' << dihScore << std::endl;
//    return _constant + _distWeight * distScore + _dihWeight * dihScore;
}

} /// scoring
} /// namespace jian

