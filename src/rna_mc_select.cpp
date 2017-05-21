#include "rna_mc_select.hpp"

BEGIN_JN

static void set_selected(DHMC &m) {
    Int l = size(m._seq);
    Int type = m.m_selected_mvel.t;

    m.m_selected.resize(l);
    for (Int i = 0; i < l; i++) m.m_selected[i] = false;

    if (type == DHMC_MVEL_HELIX) {
        const SSE & sse = *(SSE *) m.m_selected_mvel.p;
        auto & bp1 = sse.helix.front();
        for (Int i = bp1.res1.num - 1; i <= bp1.res2.num - 1; i++) m.m_selected[i] = true;
    }
    else if (type == DHMC_MVEL_LOOP) {
        const SSE & sse = *(SSE *) m.m_selected_mvel.p;
        Int a = sse.loop.front().num - 1;
        Int b = sse.loop.back().num - 1;
        for (Int i = a; i <= b; i++) m.m_selected[i] = true;
    }
    else if (type == DHMC_MVEL_FRAG) {
        const Array<Int, 2> & arr = *(Array<Int, 2> *) m.m_selected_mvel.p;
        for (Int i = arr[0]; i <= arr[1]; i++) m.m_selected[i] = true;
    }
    else if (type == DHMC_MVEL_FRAG3) {
        const Array<Int, 2> & arr = *(Array<Int, 2> *) m.m_selected_mvel.p;
        for (Int i = arr[0]; i <= arr[1]; i++) m.m_selected[i] = true;
    }
    else {
        std::cout << to_str("Unknown mvel type: ", type, "!") << std::endl;
        throw "error!";
    }
}

void dhmc_select(DHMC &m) {
    m.m_selected_mvel = m.m_mvels[int(rand() * size(m.m_mvels))];
    set_selected(m);
}

END_JN

