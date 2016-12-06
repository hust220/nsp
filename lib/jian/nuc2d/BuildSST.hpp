#ifndef JIAN_NUC2D_BUILDSST_H
#define JIAN_NUC2D_BUILDSST_H

#include "SST.hpp"
#include "NASS.hpp"

BEGIN_JN
namespace nuc2d {

class BuildSST {
private:
    // Finite State Machine: Row index is the old state, col index is the type id, entry is the new state
    std::vector<std::vector<int>> _fsm {{1, 0, -1}, {1, 2, 3}, {1, 2, 3}};
    std::map<char, int> _type_id;
    char _left, _right;

public:
    BuildSST() {
        set_primary_keys('(', ')');
    }

    void set_primary_keys(char left, char right) {
        for (auto &&pair: NASS::instance().paired_keys) {
            _type_id[pair.first] = 1;
            _type_id[pair.second] = 1;
        }
        for (auto &&ch: NASS::instance().unpaired_keys) _type_id[ch] = 1;
        _left = left; _right = right;
        _type_id[left] = 0; _type_id[right] = 2;
    }

    SST operator ()(const S &seq, const S &ss) {
        return build_sst(seq, ss);
    }

    SST build_sst(const S &seq, const S &ss) {
        Frag integral = Frag(seq, ss);
        std::list<SST> ssts;
        extract_hairpins(integral, ssts);
        return ssts.front();
    }
    
    void extract_hairpins(Frag &frag, std::list<SST> &ssts) {
        SSE hairpin;
        int num_sons;
        std::tie(hairpin, num_sons) = extract_first_hairpin(frag, first_hairpin_index(frag));
        if (hairpin._loop.size() == 2 && hairpin._loop.front().empty() && hairpin._loop.back().empty()) return;
        SST sst(hairpin);
        for (int i = 0; i < num_sons; i++) {
            sst._sons.push_front(std::move(ssts.back()));
            ssts.pop_back();
        }
        ssts.push_back(std::move(sst));
        extract_hairpins(frag, ssts);
    }

    std::pair<SSE, int> extract_first_hairpin(Frag &frag, const std::tuple<int, int, int> &tuple) {
        auto &left_index = std::get<0>(tuple);
        auto &right_index = std::get<1>(tuple);
        auto &len = std::get<2>(tuple);
        SSE hairpin;
        hairpin._helix[0].reserve(len);
        hairpin._helix[1].reserve(len);
        std::copy(std::next(frag.begin(), left_index - len + 1), std::next(frag.begin(), left_index + 1), std::back_inserter(hairpin._helix[0]));
        std::copy(std::next(frag.begin(), right_index), std::next(frag.begin(), right_index + len), std::back_inserter(hairpin._helix[1]));
        int index = left_index;
        int num_sons = 0;
        for (int i = left_index + 1; i < right_index; i++) {
            if (frag[i].num == -1) {
                num_sons++;
                Frag temp_frag;
                temp_frag.reserve(i - index - 1);
                std::copy(std::next(frag.begin(), index + 1), std::next(frag.begin(), i), std::back_inserter(temp_frag));
                hairpin._loop.push_back(temp_frag);
                index = i;
            }
        }
        Frag temp_frag;
        std::copy(std::next(frag.begin(), index + 1), std::next(frag.begin(), right_index), std::back_inserter(temp_frag));
        hairpin._loop.push_back(temp_frag);
        frag.erase(frag.begin() + left_index - len + 1, frag.begin() + right_index + len);
        frag.insert(frag.begin() + left_index - len + 1, Res{-1, 'X', '.'});
        return std::make_pair(hairpin, num_sons);
    }

    std::tuple<int, int, int> first_hairpin_index(const Frag &frag) {
        int state = 0, init_state = 0, final_state = 3, error_state = -1;
        int index = 0, left_index = -1, right_index = frag.size();
        for (auto &&res: frag) {
            auto new_state = _fsm[state][_type_id[res.ss]];
            if (state == 1 && new_state > state) {
                left_index = index - 1;
            } 
            if (new_state == final_state) {
                right_index = index;
                break;
            } else if (new_state == error_state) {
                throw "jian::BuildSST::extract_hairpins(Frag &, std::vector<SSE> &) error!";
            }
            state = new_state;
            index++;
        }
        if (left_index == -1) {
            return std::make_tuple(left_index, right_index, 0);
        } else {
            int len = 1;
            while (left_index >= len && frag[left_index - len].ss == _left &&  frag[right_index + len].ss == _right) len++;
            return std::make_tuple(left_index, right_index, len);
        }
    }
};

} // namespace nuc2d
END_JN

#endif

