#include <array>
#include <list>
#include <string>
#include <iostream>
#include <numeric>
#include <vector>
#include <deque>
#include <algorithm>
#include "rss.hpp"
#include "rss_ss_pred.hpp"
#include "fsm.hpp"
#include "log.hpp"
#include "rssp_hp_pred.hpp"

BEGIN_JN
namespace lrsp {
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

        void hp_pred(const seq_t & seq, ss_t & ss, indices_t & indices, int left, int right) {
            //LOGV << "hp_pred: " << ss << std::endl;
            seq_t seq_new;
            for (int i = left+1; i < right; i++) {
                if (indices[i] != -1) {
                    seq_new += seq[indices[i]];
                } else {
                    seq_new += 'X';
                }
            }
            LOGV << "seq_new: " << seq_new << std::endl;
            ss_t ss_new = ss_pred(seq_new);
            LOGV << "ss_new: " << ss_new << std::endl;
            for (auto && i : indices) LOGV << i << ' '; LOGV << std::endl;
            for (int i = left+1; i < right; i++) {
                if (indices[i] != -1) {
                    ss[indices[i]] = ss_new[i-left-1];
                }
            }
            LOGV << "hp_pred done" << std::endl;
            //LOGV << "hp_pred: " << ss << std::endl;
        }

        void fsm_init(Int &left, Int &right, Int &state_o) {
            left = -1;
            right = -1;
            state_o = 0;
            fsm.init();
        }

        bool hp_pred(const seq_t & seq, ss_t & ss, ss_t & ss_n, indices_t & indices) {
            int flag = 0;
            int left, right, state_o;
            std::deque<std::array<int, 2>> dq;

            LOGV << ss_n << std::endl;

            fsm_init(left, right, state_o);
            int i = 0;
            for (auto && n : indices) {
                fsm.accept(ss_n[i]);
                if (state_o == 1 && fsm.state > 1) {
                    left = i-1;
                }
                if (state_o < 3 && fsm.state == 3) {
                    right = i;
                    hp_pred(seq, ss, indices, left, right);
                    flag++;
                    dq.push_back({left, right});
                    fsm_init(left, right, state_o);
                } else {
                    state_o = fsm.state;
                }
                i++;
            }

            LOGV << "flag: " << flag << std::endl;
            if (flag == 0) {
                hp_pred(seq, ss, indices, -1, indices.size());
                return false;
            } else {
                for (auto && arr : dq) {
                    int len = hp_len(ss_n, arr[0], arr[1]);
                    arr[0] -= len-1;
                    arr[1] += len-1;
                    LOGV << arr[0] << ' ' << arr[1] << std::endl;
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

    } // namespace ss_hp_pred_detail

    S ss_hp_pred(S seq, S ss) {
        using namespace ss_hp_pred_detail;

        Int l = size(ss);

        Str pks = ss;
        for (auto && c : pks) if (c == '.' || c == '(' || c == ')') c = 'X';
        LOG << "ss_hp_pred: " << ss << std::endl;
        for (Int i = 0; i < l; i++) if (pks[i] != 'X') ss[i] = 'X';
        LOG << "ss_hp_pred: " << ss << std::endl;

        ss_t ss_new = ss;
        indices_t indices(seq.size());
        std::iota(indices.begin(), indices.end(), 0);
        for (Int i = 0; i < l; i++) if (pks[i] != 'X') indices[i] = -1;
        while (hp_pred(seq, ss, ss_new, indices)) {
            LOG << "ss_hp_pred: " << ss << std::endl;
        }

        for (Int i = 0; i < l; i++) if (pks[i] != 'X') ss[i] = pks[i];
        LOG << "ss_hp_pred: " << ss << std::endl;

        return ss;
    }

} // lrsp

END_JN

