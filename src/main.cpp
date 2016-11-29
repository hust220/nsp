#include "nsp.hpp"
#include <jian/utils/string.hpp>
#include <jian/utils/exception.hpp>
#include <jian/utils/Env.hpp>

#ifdef JN_PARA

#include <jian/mpi.hpp>
#define JN_DECLARE_MPI std::shared_ptr<MPI> g_mpi
#define JN_INIT_MPI jian::g_mpi.reset(new jian::MPI)

#else

#define JN_DECLARE_MPI 
#define JN_INIT_MPI

#endif

namespace jian {
	JN_DECLARE_MPI;
	int g_argc;
	char **g_argv;
}

int main(int argc, char **argv) {
	jian::g_argc = argc;
	jian::g_argv = argv;
	JN_INIT_MPI;
	try {
		jian::NSP::run(argc, argv);
	}
	catch (const jian::Error &inf) {
		std::cout << inf.what() << std::endl;
	}
	catch (const char * inf) {
		std::cout << inf << std::endl;
	}
	catch (const jian::str_t &s) {
		std::cout << s << std::endl;
	}
}

