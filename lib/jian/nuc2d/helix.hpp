#pragma once

#include "../pp.hpp"
#include <iostream>
#include <string>
#include <list>
#include <sstream>
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

    std::list<int> nums() const {
        std::list<int> nums;
        HELIX_EACH(*this, 
            nums.insert(std::next(nums.begin(), N_BP), {BP->res1.num-1, BP->res2.num-1})
        );
        return nums;
    }

    operator std::string() const {
//        std::string str;
        std::ostringstream stream;
        stream << seq() << ' ' << ss();
        for (auto && i : nums()) {
            stream << ' ' << i;
        }
        return stream.str();
//        stream >> str;
//        std::cout << "helix: " << str << std::endl;
//        return str;
    }

};

inline std::ostream &operator <<(std::ostream &out, const helix &h) {
    out << (std::string)h;
}

} // namespace jian

