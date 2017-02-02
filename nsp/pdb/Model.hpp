#pragma once

#include "Chain.hpp"
#include <string>
#include "jian/utils/traits.hpp"

#define EACH_INDEX_RES_HELPER(n, m, c) ({\
    int N_CHAIN##n = 0; \
    int N_RES##n = 0; \
    for (auto && CHAIN##n : m) {\
        for (auto && RES##n : CHAIN##n) {\
            c;\
            N_RES##n++;\
        }\
        N_CHAIN##n++;\
    }\
    N_RES##n;\
})
#define EACH_RES1(m, c) EACH_INDEX_RES_HELPER(1, m, c)
#define EACH_RES2(m, c) EACH_INDEX_RES_HELPER(2, m, c)
#define EACH_RES3(m, c) EACH_INDEX_RES_HELPER(3, m, c)
#define EACH_RES4(m, c) EACH_INDEX_RES_HELPER(4, m, c)
#define EACH_RES(m, c)  EACH_INDEX_RES_HELPER( , m, c)

#define EACH_ATOM_HELPER(n, m, c) ({\
    int N_ATOM##n = 0;\
    int N_RES##n = 0; \
    int N_CHAIN##n = 0; \
    for (auto && CHAIN##n : m) {\
        for (auto && RES##n : CHAIN##n) {\
            for (auto && ATOM##n : RES##n) {\
                c;\
                N_ATOM##n++;\
            }\
            N_RES##n++;\
        }\
        N_CHAIN##n++;\
    }\
    N_ATOM##n;\
})
#define EACH_ATOM1(m, c) EACH_ATOM_HELPER(1, m, c)
#define EACH_ATOM2(m, c) EACH_ATOM_HELPER(2, m, c)
#define EACH_ATOM3(m, c) EACH_ATOM_HELPER(3, m, c)
#define EACH_ATOM4(m, c) EACH_ATOM_HELPER(4, m, c)
#define EACH_ATOM(m, c)  EACH_ATOM_HELPER( , m, c)

BEGIN_JN

	class Model : public std::deque<Chain> {
	public:
		S name = "unknown";
		S type = "unknown";
		int num = 1;
		S m_cg = "aa";

		JN_DEF_ATOMS;
		JN_DEF_RESIDUES;
		JN_DEF_CHAINS;
	};

#define JN_DEF_MODELS \
	refs<Model> models() { return refs<Model>().append(*this); }\
	refs<const Model> models() const { return refs<const Model>().append(*this); }

END_JN

