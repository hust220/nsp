#include <jntest/UnitTest.hpp>
#include <jian/mpi.hpp>

int main(int argc, char **argv) {
    JN_MPI_INIT(argc, argv);
    run_test()
    JN_MPI_FREE;
}

