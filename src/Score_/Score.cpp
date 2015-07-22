#include "Score.h"

namespace jian {

Score::Score() {
    string lib = getenv("NSP");
    par_dist_obs = lib + "/RNA/par_dist";
    par_dist_ref = "";
    dihPar = lib + "/RNA/par_dih";

    distAnal = new DistAnal(cutoff, _dist_bin, reference_state);
    dihAnal = new DihAnal(_dih_bin);

    distAnal->readObsParm(par_dist_obs);
    distAnal->readRefParm(par_dist_ref);
    dihAnal->readParm(dihPar);
}

Score::Score(char *parm) {
    string lib = getenv("NSP");
    par_dist_obs = lib + "/RNA/par_dist";
    par_dist_ref = "";
    dihPar = lib + "/RNA/par_dih";

    ifstream ifile(parm);
    while (ifile) {
        string line;
        getline(ifile, line);
        vector<string> splited_line;
        tokenize(line, splited_line, " ,:");
        if (splited_line.size() != 2) continue;
        if (splited_line[0] == "par_dist_obs") {
            par_dist_obs = splited_line[1];
        } else if (splited_line[0] == "par_dist_ref") {
            par_dist_ref = splited_line[1]; 
        } else if (splited_line[0] == "par_dih") {
            dihPar = splited_line[1]; 
        } else if (splited_line[0] == "cutoff") {
            cutoff = atoi(splited_line[1].c_str());
        } else if (splited_line[0] == "dist_weight") {
            _distWeight = atof(splited_line[1].c_str());
        } else if (splited_line[0] == "dih_weight") {
            _dihWeight = atof(splited_line[1].c_str());
        } else if (splited_line[0] == "reference_state") {
            reference_state = splited_line[1];
        }
    }
    ifile.close();

    distAnal = new DistAnal(cutoff, _dist_bin, reference_state);
    distAnal->readObsParm(par_dist_obs);
    distAnal->readRefParm(par_dist_ref);

    dihAnal = new DihAnal(_dih_bin);
    dihAnal->readParm(dihPar);
}

double Score::operator ()(const RNA &rna) {
    distAnal->readRNA(rna);
    dihAnal->readRNA(rna);
    double distScore = distAnal->scoring();
    double dihScore = dihAnal->scoring();
    // cerr << rna->name << ' ' << distScore << ' ' << dihScore << endl;
    // double *score = distAnal->getScore();
    // cerr << rna->name << ' ' << score[0] << ' ' << score[1] << ' ' << score[2] << ' ' << score[3] << ' ' << score[4] << ' ' << dihScore << endl;
    std::cout << distScore << ' ' << dihScore << std::endl;
    return _constant + _distWeight * distScore + _dihWeight * dihScore;
}

} /// namespace jian

