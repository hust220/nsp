#ifndef RES_H
#define RES_H

#include "../util.h"

namespace jian {
namespace nuc2d {

class res {
public:
    char type = 'X';
    char name = 'X';
    int num = 0;
    res *next = NULL;

    res() {}
    res(char c, int i) : type(c), num(i) {}
    res(char c1, char c2, int i) : type(c1), name(c2), num(i) {}
    res(res *r) {
        if (r != NULL) {
            type = r->type;
            name = r->name;
            num = r->num;
            next = r->next;
        } else {
            cerr << "res::res(res *) error! The pointer is NULL" << endl;
            exit(1);
        }
    }
    res(const res &r) : type(r.type), name(r.name), num(r.num), next(r.next) {}
    res &operator =(const res &r) {
        type = r.type;
        name = r.name;
        num = r.num;
        next = r.next;
        return *this;
    }
    static res *copy(res *head) {
        if (head == NULL) {
            return NULL;
        }
        res *res_ = new res(head);
        res_->next = copy(head->next);
        return res_;
    }
    static void del(res *head) {
        if (head == NULL) {
            return;
        }
        del(head->next);
        delete head;
        return;
    }

};

} // namespace nuc2d
} // namespace jian
#endif

