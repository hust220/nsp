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
    // make tree with no broken tag
    void make(const std::string &seq, const std::string &ss, int hinge = 2);
    // make tree with broken tag
    void make_b(const std::string &seq, const std::string &ss, int hinge = 2);
private:
    SSTreeImpl *_impl;
};

} // namespace jian

