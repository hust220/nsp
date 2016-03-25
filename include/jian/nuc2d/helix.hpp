#pragma once

#include "bp.hpp"

#define HELIX_EACH(h, c) ({\
    int N_BP = 0; TYPEOF((h).head) BP = (h).head;\
    for (; BP != NULL; BP = BP->next) {\
        c;\
        N_BP++;\
    }\
    N_BP;\
})

namespace jian {

class helix {
public:
    bp *head = NULL;

    bool empty() const { return head == NULL; }

    helix() {head = NULL;}

    helix(const helix &h) : head(bp::copy(h.head)) { }

    helix &operator =(const helix &h) {
        bp::del(head);
        head = bp::copy(h.head);
    }

    ~helix() {
        bp::del(head);
    }

    void push_back(bp *m) {
        if (head == NULL) head = m;
        else HELIX_EACH(*this, if (BP->next == NULL) {BP->next = m; m->prev = BP; break;});
    }

    void push_front(bp *m) {
        if (head == NULL) head = m;
        else {head->prev = m; m->next = head; head = m;}
    }

    int len() const {
        return HELIX_EACH(*this,);
    }

    friend std::ostream &operator <<(std::ostream &out, const helix &h) {
        out << "Helix: " << h.seq() << ' '<< h.ss() << ' ';
        std::list<int> d; 
        HELIX_EACH(h, 
            d.insert(std::next(d.begin(), N_BP), {BP->res1.num, BP->res2.num});
        );
        for (auto && i : d) {
            out << i << ' ';
        }
    }

    std::string ss() const {
        std::string str;
        HELIX_EACH(*this, str += '(');
        HELIX_EACH(*this, str += ')');
        return str;
    }

    std::string seq() const {
        std::string str;
        HELIX_EACH(*this, str.insert(N_BP, {BP->res1.name, BP->res2.name}));
        return str;
    }

};

} // namespace jian

