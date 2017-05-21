#include "nsp.hpp"
#include "string.hpp"
#include "exception.hpp"
#include "env.hpp"
#include "mpi.hpp"

int main(int argc, char **argv) {
	JN_MPI_INIT(argc, argv);
#ifdef NDEBUG
	try {
		JN_ NSP::run(argc, argv);
	}
	catch (const jian::Error &inf) {
		Err << inf.what() << Endl;
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
    JN_MPI_FREE;
}

