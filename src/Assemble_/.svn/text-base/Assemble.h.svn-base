#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include "loop.h"
#include "../Pdb.h"
#include "../LoopModelling.h"
#include "Mol2D.h"
#include "Connect.h"

class Assemble {
public:
	Assemble(char *);
	Assemble(string, string);
	~Assemble();
	
	void log(string);
	void run();
	RNAs *findRNA(loop *, int = 1);
	RNAs *findLoop(loop *, int = 1);
	RNA *findHelix(helix *);
	RNA *createHelix(helix *);
	RNA *createHelix(string);
	RNA *createStrand(int);
	//static RNA *connect(RNA *, RNA *, int, int);
	void addP(RNA *);
	void mutate(RNA *);

	string job_name();
	string output_path();
	RNAs *getResults();

	Connect connect;
	Mol2D *mol = NULL;
	RNAs *results;
	string lib;
	string ss;
	string seq;
	string job;
	string out;
	string logfile;
	string family;
	int num;

	int is_test = 0; // Is this a test case or not?
	int _view = 0;
};



#endif




