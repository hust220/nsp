#ifndef TESTSIDECHAIN_H
#define TESTSIDECHAIN_H

#include "mymol.h"

class TestSideChain {
public:
  TestSideChain(string);
  TestSideChain(char *);
  void init();
  void rotateSideChain();
  void run();
private:
  int len;
  int *type;
  double *oldChi;
  double *chi;
  int *baseLen;
  Point **base;
  Point *c2;
  DistAnal distAnal;
  DihAnal dihAnal;
  string seq;
  Obj<RNA> oldRNA;
  Obj<RNA> newRNA;
  Obj<Score> score;
//  RNA *oldRNA;
//  RNA *newRNA;
};

#endif



