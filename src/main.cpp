#include "nsp.hpp"
#include <jian/utils/exception.hpp>
#include <jian/utils/Env.hpp>

#ifdef JN_PARA

#include <jian/mpi.hpp>
namespace jian {
	std::shared_ptr<MPI> g_mpi;
	int g_argc;
	char **g_argv;
}

#endif

int main(int argc, char **argv) {
	jian::g_argc = argc;
	jian::g_argv = argv;
#ifdef JN_PARA
	jian::g_mpi.reset(new jian::MPI);
#endif
	try {
        jian::NSP::run(argc, argv);
    } catch (const jian::Error &inf) {
        std::cout << inf.what() << std::endl;
    } catch (const char * inf) {
        std::cout << inf << std::endl;
    } catch (const std::string &s) {
        std::cout << s << std::endl;
    }
}

