#pragma once

#include "jian.hpp"

namespace jian {

using Bp = Array<Int, 2>;
using Bps = Deque<Bp>;

Bps ss_to_bps(const Str &ss);
Str bps_to_ss(const Bps &bps, int len);

}

