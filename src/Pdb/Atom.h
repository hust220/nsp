#ifndef ATOM_H_INCLUDED
#define ATOM_H_INCLUDED

#include "../Utils.h"

namespace jian {

class Atom {
public:
    Atom() {}
    Atom(Atom *atom) : name(atom->name), resName(atom->resName), rnaName(atom->rnaName), 
                            line(atom->line), num(atom->num), mass(atom->mass), 
                            x(atom->x), y(atom->y), z(atom->z), flag(atom->flag) {}
    Atom(const Atom &atom) : name(atom.name), resName(atom.resName), rnaName(atom.rnaName), 
                            line(atom.line), num(atom.num), mass(atom.mass), 
                            x(atom.x), y(atom.y), z(atom.z), flag(atom.flag) {}
    Atom &operator =(const Atom &atom) {
        name = atom.name;
        resName = atom.resName;
        rnaName = atom.rnaName;
        line = atom.line;
        num = atom.num;
        mass = atom.mass;
        x = atom.x;
        y = atom.y;
        z = atom.z;
        flag = atom.flag;
    }
    Atom(const Point &point, const string &name) {
        x = point.x;
        y = point.y;
        z = point.z;
        this->name = name;
    }
    Atom(const string &name, double x, double y, double z) {
        this->name = name;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Atom(string &, string);
    Atom(Point, string, string, int);
    friend ostream &operator <<(ostream &, const Atom &);
    Point *coord();

    double &operator [](int i) {
        if (i == 0) {
            return x;
        } else if (i == 1) {
            return y;
        } else if (i == 2) {
            return z;
        }
    }

    const double &operator [](int i) const {
        if (i == 0) {
            return x;
        } else if (i == 1) {
            return y;
        } else if (i == 2) {
            return z;
        }
    }

    template<class T = Point>
    T pos() const {
        T t;
        t[0] = x;
        t[1] = y;
        t[2] = z;
        return t;
    }

    double dist(Atom &);

    string name;
    string resName;
    string rnaName;
    string line;
    int num;
    double mass;
    double x, y, z;
    int flag;
};

} /// namespace jian
#endif // ATOM_H_INCLUDED
