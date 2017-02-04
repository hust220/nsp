#include <jian/test.hpp>
#include <jian/mpi.hpp>

int main(int argc, char **argv) {
    JN_MPI_INIT(argc, argv);
    jian::run_test()
    JN_MPI_FREE;
}

