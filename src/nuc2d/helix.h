#ifndef HELIX_H
#define HELIX_H

#include "bp.h"

namespace jian {
namespace nuc2d {

class helix {
public:
    bp *head = NULL;

    helix() { }

    helix(helix *h) : head(bp::copy(h->head)) { }

    helix(const helix &h) : head(bp::copy(h.head)) { }

    helix(helix &&h) {
        std::swap(head, h.head);
    }

    helix &operator =(const helix &h) {
        bp::del(head);
        head = bp::copy(h.head);
    }

    helix &operator =(helix &&h) {
        bp::del(head);
        std::swap(head, h.head);
    }

    ~helix() {
        bp::del(head);
    }

    template<typename Fn> void foreach(Fn &&f) {
        for (auto r = head; r != NULL; r = r->next) f(r);
    }

    template<typename Fn> void foreach(Fn &&f) const {
        for (auto r = head; r != NULL; r = r->next) f(r);
    }

    bool empty() const {
        return head == NULL;
    }

    inline int len() const {
        return getLen();
    }

    inline int size() const {
        return getLen();
    }

    inline std::string getSS() const {
        return ss();
    }

    int getLen() const {
        bp *b;
        int temp;

        for (b = head, temp = 0; b != NULL; b = b->next) {
            temp++;
        }
        return temp;
    }

    std::string ss() const {
        std::string str;
        for (bp *b = head; b != NULL; b = b->next) {
            str += '(';
        }
        for (bp *b = head; b != NULL; b = b->next) {
            str += ')';
        }
        return str;
    }

    std::string seq() const {
        string str;
        int i = 0;
        for (auto b = head; b != NULL; b = b->next) {
            str.insert(i, string() + b->res1.name + b->res2.name);
            i++;
        }
        return str;
    }

};

//class helixInfo {
//public:
//    string name;
//    int n;
//    string seq;
//    string ss;
//    string src;
//};
//
//class helixInfoEx {
//public:
//    helixInfo *hi;
//    helixInfoEx *next;
//    double score;
//
//    helixInfoEx() {
//        hi = NULL;
//        next = NULL;
//        score = 0;
//    }
//
//    helixInfoEx(helixInfo *hi, double score) {
//        this->hi = hi;
//        next = NULL;
//        this->score = score;
//    }
//};
//
//class helixInfoList {
//public:
//    helixInfoList() {
//        head = NULL;
//    }
//
//    int getLen() {
//        helixInfoEx *hie;
//        int temp = 0;
//        for (hie = head; hie != NULL; hie = hie->next) {
//            temp++;
//        }
//        return temp;
//    }
//
//    void add(helixInfo *hi, double score) {
//        helixInfoEx *hie = new helixInfoEx(hi, score);
//        if (getLen() == 0) {
//            head = hie;
//            return;
//        }
//        if (score > head->score) {
//            hie->next = head;
//            head = hie;
//            return;
//        }
//        helixInfoEx *h;
//        for (h = head; h->next != NULL; h = h->next) {
//            if (score > h->next->score) {
//                hie->next = h->next;
//                h->next = hie;
//                return;
//            }
//        }
//        if (h->next == NULL) {
//            h->next = hie;
//        }
//    }
//
//
//    helixInfoEx *head;
//};
//
} /// namespace nuc2d
} /// namespace jian

#endif

