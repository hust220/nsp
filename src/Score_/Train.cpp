#include "Train.h"

namespace jian {

Train::Train() {
  cutoff = 20;
  bin = 0.5;
  par_dist = "";
  par_dih = "";
}

Train::Train(char *parm) {
  cutoff = 20;
  bin = 0.5;
  par_dist = "";
  par_dih = "";

  ifstream ifile(parm);
  while (ifile) {
    string line;
    getline(ifile, line);
    vector<string> splited_line;
    tokenize(line, splited_line, " ,:");
    if (splited_line.size() != 2) {
      continue;
    }
    if (splited_line[0] == "cutoff") {
      cutoff = atoi(splited_line[1].c_str());
    } else if (splited_line[0] == "bin_width") {
      bin = atof(splited_line[1].c_str());
    } else if (splited_line[0] == "par_dist") {
      par_dist = splited_line[1];
    } else if (splited_line[0] == "par_dih") {
      par_dih = splited_line[1];
    }
  }
  ifile.close();
}

void Train::dist(char *filename) {
  string str;
  DistAnal distAnal(cutoff, bin);
  if (par_dist != "") {
    distAnal.readObsParm(par_dist);
  }
  ifstream ifile(filename);
  int n = 0;
  while (ifile >> str) {
    n++;
    cerr << n << ". train: " << str << endl;
    distAnal.readRNA(RNA(str));
    distAnal.train();
  }
  ifile.close();
  distAnal.printObsParm();
}

void Train::dih(char *filename) {
  string str;
  DihAnal dihAnal;
  if (par_dih != "") {
    dihAnal.readParm(par_dih);
  }
  ifstream ifile(filename);
  int n = 0;
  while (ifile >> str) {
    n++;
    cerr << n << ". train: " << str << endl;
    dihAnal.readRNA(RNA(str));
    dihAnal.train();
  }
  ifile.close();
  dihAnal.printParm();
}

} /// namespace jian




