#pragma once

#include <string>

namespace jian {

class loop;
class SSTreeImpl;

class SSTree {
public:
    SSTree();
    ~SSTree();
    loop *&head();
    bool empty() const;
    void make(const std::string &seq, const std::string &ss, int hinge = 2);
private:
    SSTreeImpl *_impl;
};

} // namespace jian

