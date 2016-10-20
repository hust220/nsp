#pragma once

namespace jian {

class res {
public:
    char type = 'X';
    char name = 'X';
    int num = 0;
    res *prev = NULL;
    res *next = NULL;

    res() {}

    res(char c, int i) : type(c), num(i) {}
    res(char c1, char c2, int i) : type(c1), name(c2), num(i) {}
    res(const res &r) : type(r.type), name(r.name), num(r.num), next(r.next) {}
    res &operator =(const res &r) {
        type = r.type;
        name = r.name;
        num = r.num;
        next = r.next;
        return *this;
    }

    static res *copy(res *head) {
        if (head == NULL) { return NULL; }
        res *p = new res(*head);
        p->next = copy(head->next);
        if (p->next != NULL) {
            p->next->prev = p;
        }
        return p;
    }

    static void del(res *head) {
        if (head == NULL) { return; }
        del(head->next);
        delete head;
        return;
    }

};

} // namespace jian

