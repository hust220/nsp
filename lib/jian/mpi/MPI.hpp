#pragma once

#ifdef JN_PARA

#include <string>
#include <memory>
#include "mpi.h"
#include "../utils/Env.hpp"

BEGIN_JN

void mpi_init(int *argc, char ***argv);

int mpi_rank();

int mpi_size();

void mpi_send(Str s, int rank);

Str mpi_recv(int rank);

void mpi_free();

END_JN

#define JN_MPI_INIT(argc, argv) JN_ mpi_init(&argc, &argv)
#define JN_MPI_FREE

#else

#define JN_MPI_INIT(argc, argv) 
#define JN_MPI_FREE

#endif
