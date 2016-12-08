#include <memory>
#include "../utils/traits.hpp"
#include "MPI.hpp"

BEGIN_JN

#ifdef JN_PARA
SP<MPI> g_mpi;
#endif

END_JN