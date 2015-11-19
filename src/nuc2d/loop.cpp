#include "loop.h"

using namespace jian;
using namespace jian::nuc2d;

int loop::getLen() {
    int temp;
    res *tempRes;

    for (temp = 0, tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
        if (tempRes->type == '(' || tempRes->type == ')' || tempRes->type == '.' || tempRes->type == '[' || tempRes->type == ']') {
            temp++;
        }
    }
    return temp;
}

int loop::size() {
    int temp;
    res *tempRes;

    for (temp = 0, tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
        temp++;
    }
    return temp;
}

int loop::getType() {
    int temp = 0;
    for (res *tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
        if (tempRes->type == ')') {
            temp++;
        }
    }
    return temp;
}

int loop::isVirtual() {
    res *tempRes;
    int flag = 0;

    if (head != NULL) {
        for (tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
            if (tempRes->type == '&') {
                flag = 1;
                break;
            }
        }
        return flag;
    } else {
        return -1;
    }
}

int loop::getLoopCounts() {
    int counts = 0;
    if (son != NULL) {
        counts += son->getLoopCounts();
    }
    if (brother != NULL) {
        counts += brother->getLoopCounts();
    }
    if (head != NULL) {
        counts++;
    }
    return counts;
}

string loop::getFlag() {
    string flag;
    res *r;
    int n = 0;
    for (r = head; r->next != NULL; r = r->next);
    if (r->type == '&') {
        int temp = 0;
        stringstream ss;
        for (res *tempRes = head; tempRes != NULL; tempRes = tempRes->next) {
            if (tempRes->type == '.' || tempRes->type == '[' || tempRes->type == ']') {
                temp++;
            } else if (tempRes->type == ')') {
                if (n % 2 == 1) {
                    ss << temp << '-';
                    temp = 0;
                }
                n++;
            }
        }
        ss << temp;
        ss >> flag;
    } else {
        int temp = 0;
        stringstream ss;
        for (res *tempRes = head->next; tempRes->next != NULL; tempRes = tempRes->next) {
            if (tempRes->type == '.' || tempRes->type == '[' || tempRes->type == ']') {
                temp++;
            } else if (tempRes->type == ')') {
                if (n % 2 == 1) {
                    ss << temp << '-';
                    temp = 0;
                }
                n++;
            }
        }
        ss << temp << '-';
        ss >> flag;
    }
    return flag;
}

string loop::ss() const {
    string ss;
    for (res *r = head; r != NULL; r = r->next) {
        ss += r->type;
    }

    return ss;
}

string loop::seq() const {
    string seq;
    for (res *r = head; r != NULL; r = r->next) {
        if (r->type != '&') {
            seq += r->name;
        }
    }
    return seq;
}
string loop::getSS() const {
    string ss;
    for (res *r = head; r != NULL; r = r->next) {
        ss += r->type;
    }

    return ss;
}

string loop::getSeq() const {
    string seq;
    for (res *r = head; r != NULL; r = r->next) {
        if (r->type != '&') {
            seq += r->name;
        }
    }
    return seq;
}

