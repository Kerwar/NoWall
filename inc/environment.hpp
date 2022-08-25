#ifndef _ENVIRONMENT_HPP_
#define _ENVIRONMENT_HPP_

#include "mpi.h"

class Environment {
 public:
  Environment(int argc, char** argv) : mpi_error_(MPI_Init(&argc, &argv)){};
  ~Environment() { MPI_Finalize(); };

 private:
  int mpi_error_;
};

#endif
