#pragma once

#include "par.hpp"
#include "pdb.hpp"

namespace jian {

class DHMC;

class SampleLoop
{
public:
	SP<DHMC> m_mc;

	void init(const Chain &chain, S ss);

	Chain operator ()();
};

}
