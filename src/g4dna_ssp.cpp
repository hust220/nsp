#include "g4dna_ssp.hpp"

BEGIN_JN

using Frag = Deque<Int>;

static Deque<Frag> split_frags(Str seq, Int len) {
    Deque<Frag> frags;
    Frag f;
    for (Int i = 0; i < size(seq); i++) {
        if (seq[i] == 'G') {
            f.push_back(i);
        }
        else {
            if (i > 0 && seq[i-1] == 'G') {
                if (size(f) >= len) {
                    frags.push_back(f);
                }
                f.clear();
            }
        }
    }
    if (size(f) >= len) frags.push_back(f);
    return frags;
}

static g4dna_ss_info_t merge_ss_info(const Frag &f, g4dna_ss_info_t info) {
    info.push_front(f);
    return info;
}

static Deque<g4dna_ss_info_t> get_ss_infos(const Deque<Frag> &frags, Int last_i, Int last_j, Int len, Int num) {
    Deque<g4dna_ss_info_t> ls;
    if (num <= 0) return ls;
    Int l = size(frags);
    for (Int i = 0; i < l; i++) {
        Int l2 = size(frags[i]);
        //std::cout << l2 << std::endl;
        for (Int j = 0; j+len-1 < l2; j++) {
            if (i > last_i || (i == last_i && j > last_j+1)) {
                Frag f;
                for (Int n = 0; n < len; n++) {
                    f.push_back(frags[i][j]+n);
                }
                if (num > 1) {
                    auto ls2 = get_ss_infos(frags, i, j+len-1, len, num-1);
                    for (auto && it : ls2) {
                        ls.push_back(merge_ss_info(f, it));
                    }
                }
                else {
                    ls.push_back({f});
                }
            }
        }
    }
    return ls;
}

template<typename _Frags>
static void print_frags(_Frags &&frags) {
    for (auto && frag : frags) {
        Int i = 0;
        for (auto && n : frag) {
            std::cout << n;
            if (i+1 < size(frag)) std::cout << '-';
            i++;
        }
        std::cout << ' ';
    }
    std::cout << std::endl;
}

static void add_ss_infos(Deque<g4dna_ss_info_t> &ls, const Deque<Frag> &frags, Int len) {
    for (auto && info : get_ss_infos(frags, 0, -99, len, 4)) {
        ls.push_back(info);
    }
}

Deque<g4dna_ss_info_t> g4dna_ssp(Str seq) {
    Deque<g4dna_ss_info_t> ls;
    Int len = 1;
    while (true) {
        auto frags = split_frags(seq, len);
        //print_frags(frags);
        //std::cout << size(frags) << std::endl;
        if (frags.empty()) break;
        add_ss_infos(ls, frags, len);
        len++;
    }
    //std::cout << size(ls) << std::endl;
    return ls;
}

END_JN

