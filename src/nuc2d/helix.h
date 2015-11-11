#ifndef HELIX_H
#define HELIX_H

#include "bp.h"

namespace jian {
namespace nuc2d {

class helix {
public:
    helix() {}
    helix(helix *h) : len(h->len), head(bp::copy(h->head)) {}
    helix(const helix &h) : len(h.len), head(bp::copy(h.head)) {}
    helix &operator =(const helix &h) {
        len = h.len;
        bp::del(head);
        head = bp::copy(h.head);
    }
    ~helix() {
        bp::del(head);
    }

    string getSS() const;
    string ss() const;
    string seq() const;

    bp *begin();
    bp *end();

    int getLen();

    int len = 0;
    bp *head = NULL;
};

class helixInfo {
public:
    string name;
    int n;
    string seq;
    string ss;
    string src;
};

class helixInfoEx {
public:
    helixInfoEx(); 
    helixInfoEx(helixInfo *hi, double score); 

    helixInfo *hi;
    helixInfoEx *next;
    double score;
};

class helixInfoList {
public:
    helixInfoList(); 
    int getLen();
    void add(helixInfo *hi, double score); 

    helixInfoEx *head;
};

} /// namespace nuc2d
} /// namespace jian

#endif

