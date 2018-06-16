#pragma once

#ifdef JN_PARA

#include <string>
#include <memory>
#include "mpi.h"
#include "env.hpp"

namespace jian {

void mpi_init(int *argc, char ***argv);

int mpi_rank();

int mpi_size();

void mpi_send(Str s, int rank);

Str mpi_recv(int rank);

void mpi_free();

}

#define JN_MPI_INIT(argc, argv) jian::mpi_init(&argc, &argv)
#define JN_MPI_FREE

#else

#define JN_MPI_INIT(argc, argv) 
#define JN_MPI_FREE

#endif
