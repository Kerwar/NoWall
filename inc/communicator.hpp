#ifndef _COMMUNICATOR_HPP_
#define _COMMUNICATOR_HPP_

#include <vector>

#include "mpi.h"

using std::vector;

class Communicator {
 public:
  Communicator() : Communicator(MPI_COMM_WORLD){};
  explicit Communicator(MPI_Comm comm)
      : comm_(comm),
        mpi_error_(MPI_Comm_size(comm_, &size_)),
        mpi_error2_(MPI_Comm_rank(comm_, &rank_)){};
  ~Communicator(){};

  inline int rank() const { return rank_; };
  inline int size() const { return size_; };
  inline MPI_Comm comm() const { return comm_; };
  MPI_Comm *comm_ref() { return &comm_; };

  void wait() { MPI_Barrier(comm_); };
  
  void send(const vector<double> &info, int destination);
  void recv(vector<double> &info, int origin);

  void send_recv(vector<double> info, int other_proc);
  
  bool is_main_proc() const { return rank_ == 0 ? true : false; };
  
  void update() {
    MPI_Comm_size(comm_, &size_);
    MPI_Comm_rank(comm_, &rank_);
  }

 public:
  MPI_Comm comm_;
  int rank_, size_;
  int mpi_error_, mpi_error2_;
  MPI_Request request_;
};
#endif