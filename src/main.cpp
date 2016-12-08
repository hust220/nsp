#include "nsp.hpp"
#include <jian/utils/string.hpp>
#include <jian/utils/exception.hpp>
#include <jian/utils/Env.hpp>

#ifdef JN_PARA
#  include <jian/mpi.hpp>
#  define JN_INIT_MPI jian::g_mpi.reset(new jian::MPI)
#else
#  define JN_INIT_MPI
#endif

int main(int argc, char **argv) {
	JN_ g_argc = argc;
	JN_ g_argv = argv;
	JN_INIT_MPI;
#ifdef NDEBUG
	try {
		JN_ NSP::run(argc, argv);
	}
	catch (const jian::Error &inf) {
		Out << inf.what() << Endl;
	}
	catch (const char * inf) {
		Err << inf << Endl;
	}
	catch (const jian::Str &s) {
		Err << s << Endl;
	}
#else
	JN_ NSP::run(argc, argv);
#endif
}

