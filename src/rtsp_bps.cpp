#include "rtsp_bps.hpp"
#include "rss.hpp"

BEGIN_JN

Bps ss_to_bps(const Str &ss) {
    const NASS &nass = NASS::instance();
    const auto & pks = nass.paired_keys;
    const auto & bks = nass.break_keys;

    std::deque<std::deque<int>> dq;
    Bps bps;
    int i, j, k;

    dq.resize(pks.size());
    i = 0;
    for (auto && c : ss) {
        if (std::find(bks.begin(), bks.end(), c)!= bks.end()) continue;
        auto it1 = std::find_if(pks.begin(), pks.end(), [c](auto && p){ return p.first == c; });
        auto it2 = std::find_if(pks.begin(), pks.end(), [c](auto && p) { return p.second == c; });
        if (it1 != pks.end()) {
            j = std::distance(pks.begin(), it1);
            dq[j].push_back(i);
        }
        else if (it2 != pks.end()) {
            j = std::distance(pks.begin(), it2);
            if (size(dq[j]) > 0) {
                k = dq[j].back();
                bps.push_back({k, i});
                dq[j].pop_back();
            }
        }
        i++;
    }
    return bps;
}

static Int bp_level(const Str &ss, const Bp &bp) {
    const auto & pks = NASS::instance().paired_keys;
    Int left = bp[0], right = bp[1];
    if (ss[left] != '.' || ss[right] != '.') return -1;
    for (Int i = 0; i < size(pks); i++) {
        Int score = 0, flag = 1;
        for (Int j = left + 1; j < right; j++) {
            if (ss[j] == pks[i].first) score++;
            else if (ss[j] == pks[i].second) score--;
            if (score < 0) break;
        }
        if (score != 0) flag = 0;
        if (flag == 1) return i;
    }
    return -1;
}

Str bps_to_ss(const Bps &bps, int len) {
    Str ss(len, '.');
    const auto & pks = NASS::instance().paired_keys;
    for (auto && bp : bps) {
        Int l = bp_level(ss, bp);
        if (l >= 0) {
            ss[bp[0]] = pks[l].first;
            ss[bp[1]] = pks[l].second;
        }
    }
    return ss;
}

END_JN

