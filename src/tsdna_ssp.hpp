#include "jian.hpp"

namespace jian {

struct tsdna_ss_info_t {
    Str ss;
    Num score;
};

List<tsdna_ss_info_t> tsdna_ssp(Str seq);

}

