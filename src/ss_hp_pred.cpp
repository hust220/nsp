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

namespace ss_hp_pred_detail {

using ss_t = std::string;
using seq_t = std::string;
using indices_t = std::deque<int>;

FSM<char> fsm(
    {{0, 1, 0},
     {2, 1, 3},
     {2, 1, 3}},
    [](auto && c){
        if (c == '(') return 1;
        else if (c == ')') return 2;
        else return 0;
    }
);

int hp_len(const ss_t &ss, int left, int right) {
    int len = 1;
    while (ss[left-len+1] == '(' && ss[right+len-1] == ')') len++;
    len--;
    return len;
}

ss_t ss_pred(const seq_t & seq) {
    ss_t ss = seq;
    return seq;
}

void hp_pred(const seq_t & seq, ss_t & ss, ss_t & ss_n, indices_t & indices, int left, int right) {
    seq_t seq_new;
    for (auto && i : indices) {
        if (i != -1) {
            seq_new += seq[i];
        } else {
            seq_new += 'X';
        }
    }
    ss_t ss_new = ss_pred(seq_new);
    for (int i = 0; i < indices.size(); i++) {
        if (i != -1) {
            ss[indices[i]] = ss_new[i];
        }
    }
}

bool hp_pred(const seq_t & seq, ss_t & ss, ss_t & ss_n, indices_t & indices) {
    int flag = 0;
    int left, right, state_o;
    std::deque<std::array<int, 2>> dq;

    std::cout << ss_n << std::endl;

    auto init = [&left, &right, &state_o](){
        left = -1;
        right = -1;
        state_o = 0;
        fsm.init();
    };

    init();
    int i = 0;
    for (auto && n : indices) {
        fsm.accept(ss_n[i]);
        if (state_o == 1 && fsm.state > 1) {
            left = i-1;
        }
        if (state_o < 3 && fsm.state == 3) {
            right = i;
            hp_pred(seq, ss, ss_n, indices, left, right);
            flag++;
            dq.push_back({left, right});
            init();
        }
        state_o = fsm.state;
        i++;
    }

    std::cout << "flag: " << flag << std::endl;
    if (flag == 0) {
        return false;
    } else {
        for (auto && arr : dq) {
            int len = hp_len(ss_n, arr[0], arr[1]);
            arr[0] -= len-1;
            arr[1] += len-1;
            std::cout << arr[0] << ' ' << arr[1] << std::endl;
        }

        ss_t ss_temp;
        indices_t indices_temp;
        for (int i = 0; i < ss_n.size(); i++) {
            if (std::find_if(dq.begin(), dq.end(), [&i](auto && arr){return arr[0]==i;}) != dq.end()) {
                ss_temp += "XXXX";
                for (int j = 0; j < 4; j++) indices_temp.push_back(-1);
            } else if (std::find_if(dq.begin(), dq.end(), [&i](auto && arr){return arr[0]<i&&i<=arr[1];}) != dq.end()) {
            } else {
                ss_temp += ss_n[i];
                indices_temp.push_back(indices[i]);
            }
        }
        ss_n = ss_temp;
        indices = indices_temp;

        return true;
    }
}

ss_t ss_hp_pred(seq_t seq, ss_t ss) {
    ss_t ss_new = ss;
    indices_t indices(seq.size());
    std::iota(indices.begin(), indices.end(), 0);
    while (hp_pred(seq, ss, ss_new, indices));
    return ss;
}

}

REGISTER_NSP_COMPONENT(ss_hp_pred) {
    using namespace ss_hp_pred_detail;
    seq_t seq = par["seq"][0];
    ss_t ss = par["ss"][0];
    ss = ss_hp_pred(seq, ss);
    std::cout << ss << std::endl;
}

} // namespace jian

