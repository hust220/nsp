#ifndef MYTIME_H
#define MYTIME_H

#include "std.h"

class Time {
public:
	static string getTime();
	static int getYear();
	static int getMon();
	static int getDate();
	static int getHour();
	static int getMin();
	static int getSec();
};

#endif
