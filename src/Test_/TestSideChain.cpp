#include "TestSideChain.h"

TestSideChain::TestSideChain(char *filename) {
  oldRNA = new RNA(filename);
  init();
}

TestSideChain::TestSideChain(string filename) {
  oldRNA = new RNA(filename);
  init();
}

void TestSideChain::init() {
  len = oldRNA->getLen();
  seq = oldRNA->getSeq();
  type = new int[len];
  oldChi = new double[len];
  chi = new double[len];
  base = new Point *[len];
  baseLen = new int[len];
  c2 = new Point[len];
  for (int i = 0, n = 0; i < oldRNA->chains.size(); i++) {
    for (int j = 0; j < oldRNA->chains[i].residues.size(); j++, n++) {
      string resName = oldRNA->chains[i].residues[j].name;
      if (resName == "A") {
        type[n] = 0;
        baseLen[n] = 11;
        base[n] = new Point[11];
      } else if (resName == "U") {
        type[n] = 1;
        baseLen[n] = 9;
        base[n] = new Point[9];
      } else if (resName == "G") {
        type[n] = 2;
        baseLen[n] = 12;
        base[n] = new Point[12];
      } else if (resName == "C") {
        type[n] = 3;
        baseLen[n] = 9;
        base[n] = new Point[9];
      }
      int flag = 0;
      for (int k = 0, temp = 0; k < oldRNA->chains[i].residues[j].atoms.size(); k++) {
        string name = oldRNA->chains[i].residues[j].atoms[k].name;
        if (name == "C1*") {
          flag = 1;
        } else if (name == "C2*") {
          c2[n].x = oldRNA->chains[i].residues[j].atoms[k].x;
          c2[n].y = oldRNA->chains[i].residues[j].atoms[k].y;
          c2[n].z = oldRNA->chains[i].residues[j].atoms[k].z;
        }
        if (flag == 1) {
          base[n][temp].x = oldRNA->chains[i].residues[j].atoms[k].x;
          base[n][temp].y = oldRNA->chains[i].residues[j].atoms[k].y;
          base[n][temp].z = oldRNA->chains[i].residues[j].atoms[k].z;
          temp++;
        }
      }
      oldChi[n] = Point::dihedral(c2[n], base[n][0], base[n][1], base[n][2]);
      chi[n] = oldChi[n];
    }
  }
  /*
  string lib = getenv("RNA");
  string par_dist = lib + "par_dist";
  string par_dih = lib + "par_dih";
  distAnal.readParm(par_dist);
  dihAnal.readParm(par_dih);
  */
  score = new Score;
}

/*
double TestSideChain::score() {
  return score->run(newRNA);
  distAnal.readRNA(newRNA);
  double distScore = distAnal.scoring();
  dihAnal.readRNA(newRNA);
  double dihScore = dihAnal.scoring();
  return distScore + 0.005 * dihScore;
}
*/

void TestSideChain::rotateSideChain() {
  newRNA = oldRNA->copy();

  for (int i = 0; i < len; i++) {
    Point::rotate(base[i], baseLen[i], base[i][0], base[i][1], chi[i] - oldChi[i]);
  }
  for (int i = 0, n = 0; i < newRNA->chains.size(); i++) {
    for (int j = 0; j < newRNA->chains[i].residues.size(); j++, n++) {
      for (int k = 0, flag = 0, temp = 0; k < newRNA->chains[i].residues[j].atoms.size(); k++) {
        string name = newRNA->chains[i].residues[j].atoms[k].name;
        if (name == "C1*") {
          flag = 1;
        }
        if (flag == 1) {
          newRNA->chains[i].residues[j].atoms[k].x = base[n][temp].x;
          newRNA->chains[i].residues[j].atoms[k].y = base[n][temp].y;
          newRNA->chains[i].residues[j].atoms[k].z = base[n][temp].z;
          temp++;
        }
      }
    }
  }
}

void TestSideChain::run() {
  srand((unsigned int)time(0));
  for (int i = 0; i < len; i++) {
    chi[i] = rand() % 1000 / 1000. * 360;
  }
  double newScore, oldScore = 9999, minScore = 9999;
  int doneTimes = 0;
  Obj<RNA> minRNA;
//  RNA *minRNA = NULL;
  for (int i = 0; i < 1000; i++) {
    int n = int((rand() % 1000) / 1000. * len);
    double temp = chi[n];
    chi[n] += (rand() % 1000 / 1000. - 1) * 10;
    if (chi[n] >= 360) chi[n] -= 360;
    if (chi[n] < 0) chi[n] += 360;
    rotateSideChain();
    newScore = score->run(newRNA);
    double deltaScore = oldScore - newScore;
    if (deltaScore < 0) {
      double temp2 = (rand() % 1000) / 1000.0;
      if (temp2 > exp(deltaScore)) {
        chi[n] = temp;
        rotateSideChain();
      } else {
        cerr << newScore << '\t' << minScore << endl;
        oldScore = newScore;
        RMSD rmsd(oldRNA, newRNA);
      //  cerr << doneTimes << '\t' << rmsd.run() << endl;
        doneTimes++;
        if (minScore > newScore) {
          minScore = newScore;
          minRNA = newRNA->copy();
        }
      }
    } else {
      cerr << newScore << '\t' << minScore << endl;
      oldScore = newScore;
      RMSD rmsd(oldRNA, newRNA);
    //  cerr << doneTimes << '\t' << rmsd.run() << endl;
      doneTimes++;
      if (minScore > newScore) {
        minScore = newScore;
        minRNA = newRNA->copy();
      }
    }
  }
  minRNA->print();
}










