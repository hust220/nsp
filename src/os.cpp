#include "os.hpp"

namespace jian {

namespace os {

#if defined(JN_OS_LINUX)

# include "sys/stat.h"

#elif defined(JN_OS_WIN)
#elif defined(JN_OS_MACX)
#endif

    void rm(Str filename) {
    }

}

}

