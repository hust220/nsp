#pragma once

BEGIN_JN

class res {
public:
    char type = 'X';
    char name = 'X';
    int num = 0;

    res() {}
    res(char c, int i) : type(c), num(i) {}
    res(char c1, char c2, int i) : type(c1), name(c2), num(i) {}
    res(const res &r) : type(r.type), name(r.name), num(r.num) {}
};

END_JN

