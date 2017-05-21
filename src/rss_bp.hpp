#pragma once

#include "rss_res.hpp"

BEGIN_JN

class bp {
public:
    res res1, res2;
    bp *prev = NULL;
    bp *next = NULL;

    bp() = default;

    bp(const bp &m) {
        res1 = m.res1; res2 = m.res2;
        prev = m.prev; next = m.next;
    }

    bp(const res &res1, const res &res2) {
        this->res1 = res1; this->res2 = res2;
        prev = NULL; next = NULL;
    }

    static bp *copy(bp *b) {
        if (b == NULL) return NULL;
        bp *p = new bp(*b);
        p->next = copy(b->next);
        if (p->next != NULL) {
            p->next->prev = p;
        }
        return p;
    }

    static void del(bp *b) {
        if (b == NULL) return;
        del(b->next);
        delete b;
        return;
    }

};

END_JN

