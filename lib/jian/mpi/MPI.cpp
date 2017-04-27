#ifdef JN_PARA

#include <memory>
#include "MPI.hpp"

BEGIN_JN

void mpi_init(int *argc, char ***argv) {
    MPI_Init(argc, argv);
}

int mpi_rank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int mpi_size() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void mpi_send(Str s, int rank) {
    MPI_Send(s.c_str(), s.size(), MPI_CHAR, rank, 0, MPI_COMM_WORLD);
}

Str mpi_recv(int rank) {
    char hi[1024];
    MPI_Status status;
    //int n = MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
    MPI_Recv(hi, 1024, MPI_CHAR, rank, 0, MPI_COMM_WORLD, &status);
    int n = ((int*)&status)[0];
    return Str(hi, n);
}

void mpi_free() {
    MPI_Finalize();
}

END_JN

#endif


