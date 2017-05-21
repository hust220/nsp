#pragma once

#include "par.hpp"
#include "pdb.hpp"

BEGIN_JN

class DHMC;

class SampleLoop
{
public:
	SP<DHMC> m_mc;

	void init(const Chain &chain, S ss);

	Chain operator ()();
};

END_JN