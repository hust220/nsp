#pragma once

#include "../pp.hpp"
#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include "bp.hpp"
#include "../utils/string.hpp"

#define JN_BREAK_INIT bool is_break = false
#define JN_BREAK is_break = true; break

#define BEGIN_HELIX_EACH(h) \
	do { \
		int N_BP = 0; \
		TYPEOF((h).head) BP = (h).head; \
		for (; BP != NULL; BP = BP->next, N_BP++)  

#define END_HELIX_EACH \
	} while(0)


#define HELIX_EACH(h, c) \
	BEGIN_HELIX_EACH(h) {\
		c;\
	} END_HELIX_EACH

BEGIN_JN

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
		int l = 0;
		HELIX_EACH(*this, l++);
		return l;
    }

    Str ss() const {
        Str str;
        HELIX_EACH(*this, str += '(');
        HELIX_EACH(*this, str += ')');
        return str;
    }

    Str seq() const {
        Str str;
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

    operator Str() const {
//        Str str;
        std::ostringstream stream;
        stream << seq() << ' ' << ss();
        for (auto && i : nums()) {
            stream << ' ' << i+1;
        }
        return stream.str();
//        stream >> str;
//        std::cout << "helix: " << str << std::endl;
//        return str;
    }

    friend std::ostream &operator <<(std::ostream &out, const helix &h) {
        out << "Helix: " << ' ' << (Str)h;
    }

};

END_JN

