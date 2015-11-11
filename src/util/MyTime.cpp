#include "MyTime.h"

string Time::getTime() {
	time_t t;
	time(&t);
	string str = asctime(localtime(&t));
	return str.substr(0, str.size() - 1);
}

int Time::getYear() {
	return 0;
}

int Time::getMon() {
	return 0;
}

int Time::getDate() {
	return 0;
}

int Time::getHour() {
	int t = time(0);
	t = t % (3600 * 24);
	t = t / 3600 + 8;
	return t;
}

int Time::getMin() {
	int t = time(0);
	t = t % 3600;
	t = t / 60;
	return t;
}

int Time::getSec() {
	int t = time(0);
	t = t % 60;
	return t;
}


















