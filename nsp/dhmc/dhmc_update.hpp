#pragma once

#include "DHMC.hpp"

BEGIN_JN

void dhmc_sample_res(DHMC &m);

void update_fragment(DHMC &m);

void translate_mvel(DHMC &m);

void rotate_about_center(DHMC &m);

Bool is_fixed(DHMC &m, Int i);

END_JN

