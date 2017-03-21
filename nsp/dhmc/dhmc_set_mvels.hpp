#pragma once

BEGIN_JN

class DHMC;

void dhmc_set_mvels(DHMC &dhmc);

void mvels_set_fixed_els(Deque<MvEl *> &mvels, const Deque<Frags> &els);

END_JN

