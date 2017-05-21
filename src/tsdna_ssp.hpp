#include "jian.hpp"

BEGIN_JN

struct tsdna_ss_info_t {
    Str ss;
    Num score;
};

List<tsdna_ss_info_t> tsdna_ssp(Str seq);

END_JN

