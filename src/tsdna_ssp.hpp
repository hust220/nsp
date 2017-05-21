#include "jian.hpp"

BEGIN_JN

namespace tsdna {

    struct ss_info_t {
        Str ss;
        Num score;
    };

    List<ss_info_t> seq_ss(Str seq);

}

END_JN

