#ifndef OBJ_H
#define OBJ_H

#include "iostream"
#include "cstdlib"
using namespace std;

namespace jian {

template<class T>
class Obj {
public:
	Obj();
	Obj(T *);
	Obj(const Obj<T> &);
	~Obj();
	void init(T *);
	Obj<T> &operator =(const Obj<T> &);
	T &operator *() const;
	T *operator ->() const;
	int getCounts();
	T *get() const;
private:
	T *p;
	int *counts;
};

template<class T>
Obj<T>::Obj() {
	p = NULL;
	counts = new int(0);
}

template<class T>
Obj<T>::Obj(T *t) {
	p = t;
	counts = new int(1);
}

template<class T>
Obj<T>::Obj(const Obj<T> &_Obj) {
	p = _Obj.p;
	counts = _Obj.counts;
	(*counts)++;
}

template<class T>
void Obj<T>::init(T *t) {
	if (p != NULL) {
		cerr << "Obj init error!" << endl;
		exit(1);
	}
	p = t;
	(*counts)++;
}

template<class T>
Obj<T> &Obj<T>::operator =(const Obj<T> &_Obj) {
	if ((*_Obj.counts) == 0) {
		cerr << "Obj operator= error!" << endl;
		exit(1);
	}
	if ((*counts) == 0) {
		delete counts;
	} else {
		(*counts)--;
		if ((*counts) == 0) {
			delete p;
		}
	}
	(*_Obj.counts)++;
	p = _Obj.p;
	counts = _Obj.counts;
}

template<class T>
T &Obj<T>::operator *() const {
	return *p;
}

template<class T>
T *Obj<T>::operator ->() const {
	return p;
}

template<class T>
Obj<T>::~Obj() {
	(*counts)--;
	if ((*counts) == 0) {
		delete p;
		delete counts;
	}
}

template<class T>
int Obj<T>::getCounts() {
	return *counts;
}

template<class T>
T *Obj<T>::get() const {
	return p;
}

} /// namespace jian

#endif

