#include <array>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <jian/nuc2d.hpp>
#include <jian/utils/FSM.hpp>
#include "nsp.hpp"

namespace jian {

namespace ss_del_sm_hp_detail {

using ss_t = std::string;

FSM<char> fsm(
    {{0, 1, 0},
     {2, 1, 0},
     {2, 1, 3}},
    [](auto && c){
        if (c == '(') return 1;
        else if (c == ')') return 2;
        else return 0;
    }
);

int del_hp(ss_t &ss, int left, int right) {
    int len = 1;
    while (ss[left-len+1] == '(' && ss[right+len-1] == ')') len++;
    len--;
    if (len <= 2 && right - left - 1 < 6) {
        for (int i = left-len+1; i <= right+len-1; i++) {
            ss[i] = '.';
        }
        return 1;
    }
    return 0;
}

bool del_hp(ss_t &ss) {
    fsm.init();
    int old_state = fsm.state;
    int left, right;
    int i = 0;
    int flag = 0;
    for (auto && c : ss) {
        fsm.accept(c);
        if (old_state == 1 && fsm.state == 2) {
            left = i-1;
        } else if (fsm.state == 3) {
            right = i;
            flag += del_hp(ss, left, right);
            fsm.init();
        }
        old_state = fsm.state;
        i++;
    }
    if (flag == 0) return false; else return true;
}

bool check_ss(const ss_t & ss) {
    for (auto && c : ss) {
        if (c != '(' && c != ')' && c != '.') {
            return false;
        }
    }
    return true;
}

}

REGISTER_NSP_COMPONENT(ss_del_sm_hp) {
    using namespace ss_del_sm_hp_detail;
    ss_t ss = par["ss"][0];
    assert(check_ss(ss) && "Illegal 2D structure!");
    while (del_hp(ss));
    std::cout << ss << std::endl;
}

} // namespace jian

