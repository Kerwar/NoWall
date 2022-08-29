#ifndef PARALEL_H
#define PARALEL_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "communicator.hpp"
#include "field.hpp"
#include "mpi.h"
#include "parameters.hpp"

using std::to_string;
using std::vector;

using namespace parameters;

enum Loc {
  up,
  down,
};

struct Paralel {
 public:
  Paralel();
  virtual ~Paralel();

  int myRowId, myColId;
  int myLeft, myRight, myBot, myTop;
  Loc loc;
  int UCMainProc, DCMainProc;

  int locIStr, locIEnd, locJStr, locJEnd;
  int iStr, iEnd;
  int jStr, jEnd;

  vector<int> iStrAllProc, iEndAllProc;
  vector<int> jStrAllProc, jEndAllProc;

  int fixPointProc{0};
  Communicator channel;
  MPI_Datatype column_type;

  int id(int, int);
  void SendInfoToNeighbours(Field &vec);
  void SendInfoToCommMainProc(Field &vec, Field &sol);
  vector<double> ExchangeWallTemperature(const Field &TWall) const;
  void GatherWallTemperature(Field &TWall, Field &T) const;
  void ShareWallTemperatureInfo(Field &TWall, Field &T) const;
  void SendWallTInTheChannel(Field &TWall) const;
  void MainsProcs(Loc location, int &proc);
  void distributeToProcs(Field &sol, Field &vec);

  void setProcWithFixPoint(const bool &fixPointIntThisProc);
  void scatter(double &from, int size, double &to, int originProc);

  bool is_left2right() const;
  bool is_right2left() const;

  bool is_main_proc_up() const;
  bool is_main_proc_down() const;

  inline bool isProcNull(int &proc) const {
    return proc == MPI_PROC_NULL ? true : false;
  };
  
  inline bool isFirstProcOfRow() const { return myRowId == 0 ? true : false; };
  inline bool isLastProcOfRow() const {
    return myRowId == NPROCSINROW - 1 ? true : false;
  };
  
  inline bool isFirstProcOfCol() const { return myColId == 0 ? true : false; };
  inline bool isLastProcOfCol() const {
    return myColId == NPROCSINROW - 1 ? true : false;
  };
  string PrintMyRank();

 private:
  int worldMyProc;
  void setProcMesh();
  void GatherTemperature(Field &T, Field &TWall, int J) const;
};

#endif
