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

using std::cout;
using std::endl;
using std::to_string;

using std::vector;

enum Loc {
  up,
  down,
};

template <int NTOTAL, int MINPUT, int NPROCS>
struct Paralel {
 public:
  Paralel();
  virtual ~Paralel();

  static constexpr int nProcsInRow = NPROCS / 2;
  static constexpr int nProcsInCol = 1;
  static constexpr int NI = NTOTAL / nProcsInRow + 2;
  static constexpr int NJ = MINPUT / 2 + 2;

  int worldNProcs, worldMyProc;

  int nProcs, myProc;
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
  Communicator myComm;
  MPI_Datatype column_type;

  int id(int, int);
  inline int worldid(int row, int col) { return col + row * nProcsInCol; };
  void SendInfoToNeighbours(Field &vec);
  void SendInfoToCommMainProc(Field &vec, Field &sol);
  vector<double> ExchangeWallTemperature(const Field &TWall) const;
  void GatherWallTemperature(Field &TWall, Field &T) const;
  void ShareWallTemperatureInfo(Field &TWall, Field &T) const;
  void SendWallTInTheChannel(Field &TWall) const;
  void MainsProcs(Loc location, int &proc);
  void PrintProgress(const string &message, bool showProgress);
  void distributeToProcs(Field &sol, Field &vec);

  void setProcWithFixPoint(const bool &fixPointIntThisProc);
  void scatter(double &from, int size, double &to, int originProc);

  inline bool isLeftToRight() const { return loc == up ? true : false; };
  inline bool isRightToLeft() const { return loc == down ? true : false; };

  inline bool isProcNull(int &proc) const {
    return proc == MPI_PROC_NULL ? true : false;
  };
  inline bool isFirstProcOfRow() const { return myRowId == 0 ? true : false; };
  inline bool isLastProcOfRow() const {
    return myRowId == nProcsInRow - 1 ? true : false;
  };
  inline bool isFirstProcOfCol() const { return myColId == 0 ? true : false; };
  inline bool isLastProcOfCol() const {
    return myColId == nProcsInCol - 1 ? true : false;
  };
  string PrintMyRank();

 private:
  void setProcMesh();
  void GatherTemperature(Field &T, Field &TWall, int J) const;
  // Field::vec1dfield TIntoDerivative(const Field &TWall,const Field &TW2R,
  // double exCte, double ex2);
};

template <int NTOTAL, int MINPUT, int NPROCS>
Paralel<NTOTAL, MINPUT, NPROCS>::Paralel() {
  PROFILE_FUNCTION();
  MPI_Comm_rank(MPI_COMM_WORLD, &worldMyProc);
  MPI_Comm_size(MPI_COMM_WORLD, &worldNProcs);

  if (worldMyProc < NPROCS / 2)
    loc = up;
  else
    loc = down;

  MPI_Comm helpComm;
  MPI_Comm_split(MPI_COMM_WORLD, loc, worldMyProc, &helpComm);
  MPI_Comm_size(helpComm, &nProcs);
  int dims[2] = {nProcsInRow, nProcsInCol};

  int period[2] = {false, false};
  int reorder = true;

  MPI_Cart_create(helpComm, 2, dims, period, reorder, myComm.comm_ref());
  MPI_Comm_free(&helpComm);

  MPI_Comm_rank(myComm.comm(), &myProc);

  int myCoords[2];
  MPI_Cart_coords(myComm.comm(), myProc, 2, myCoords);

  myRowId = myCoords[0];
  myColId = myCoords[1];

  MPI_Cart_shift(myComm.comm(), 0, 1, &myLeft, &myRight);
  MPI_Cart_shift(myComm.comm(), 1, 1, &myBot, &myTop);

  MainsProcs(up, UCMainProc);
  MainsProcs(down, DCMainProc);

  switch (loc) {
    case up:
      setProcMesh();

      locIStr = 0;
      locIEnd = NTOTAL;
      locJStr = NJ - 2;
      locJEnd = MINPUT - 1;

      for (int i = 0; i < NPROCS / 2; i++) {
        jStrAllProc[i] += locJStr;
        jEndAllProc[i] += locJStr;
      }
      break;

    case down:
      setProcMesh();

      locIStr = 0;
      locIEnd = NTOTAL;
      locJStr = 0;
      locJEnd = NJ - 2;
      break;
  }

  iStr = iStrAllProc[myProc];
  iEnd = iEndAllProc[myProc];
  jStr = jStrAllProc[myProc];
  jEnd = jEndAllProc[myProc];

  MPI_Type_vector(NJ, 1, NI, MPI_DOUBLE_PRECISION, &column_type);
  MPI_Type_commit(&column_type);
}

template <int NTOTAL, int MINPUT, int NPROCS>
Paralel<NTOTAL, MINPUT, NPROCS>::~Paralel() {
  // MPI_Type_free(&column_type);
}

template <int NTOTAL, int MINPUT, int NPROCS>
int Paralel<NTOTAL, MINPUT, NPROCS>::id(int x, int y) {
  int myCoords[2] = {x, y};
  int result;

  MPI_Cart_rank(myComm.comm(), myCoords, &result);

  return result;
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::SendInfoToNeighbours(Field &vec) {
  int leftTag = 10, rightTag = 11, topTag = 12, botTag = 13;

  leftTag *= (1 + loc);
  rightTag *= (1 + loc);
  topTag *= (1 + loc);
  botTag *= (1 + loc);

  myComm.wait();
  // SENDING INFO TO THE LEFT
  MPI_Sendrecv(&vec.value[1], 1, column_type, myLeft, leftTag,
               &vec.value[NI - 1], 1, column_type, myRight, leftTag,
               myComm.comm(), MPI_STATUS_IGNORE);

  if (!isProcNull(myRight))
    for (int j = 0; j < NJ; j++)
      vec.value[NI - 1 + j * NI] =
          (vec.value[NI - 1 + j * NI] + vec.value[NI - 2 + j * NI]) * 0.5;

  myComm.wait();
  // SENDING INFO TO THE RIGHT

  MPI_Sendrecv(&vec.value[NI - 1], 1, column_type, myRight, rightTag,
               &vec.value[0], 1, column_type, myLeft, rightTag, myComm.comm(),
               MPI_STATUS_IGNORE);

  myComm.wait();

  // SENDING INFO TO THE BOT
  MPI_Sendrecv(&vec.value[NI], NI, MPI_DOUBLE_PRECISION, myBot, botTag,
               &vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop,
               botTag, myComm.comm(), MPI_STATUS_IGNORE);

  if (!isProcNull(myTop))
    for (int i = 0; i < NI; i++)
      vec.value[i + (NJ - 1) * NI] =
          (vec.value[i + (NJ - 1) * NI] + vec.value[i + (NJ - 2) * NI]) * 0.5;
  myComm.wait();

  // SENDING INFO TO THE TOP
  MPI_Sendrecv(&vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop,
               topTag, &vec.value[0], NI, MPI_DOUBLE_PRECISION, myBot, topTag,
               myComm.comm(), MPI_STATUS_IGNORE);
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::MainsProcs(Loc location, int &proc) {
  PROFILE_FUNCTION();
  int buffer;

  if (location == loc) {
    if (myProc == 0) buffer = worldMyProc;
    PMPI_Bcast(&buffer, 1, MPI_INT, 0, myComm.comm());
  }

  int sendProc;
  switch (location) {
    case up:
      sendProc = 0;
      break;

    case down:
      sendProc = NPROCS / 2;
      break;
  }
  PMPI_Bcast(&buffer, 1, MPI_INT, sendProc, MPI_COMM_WORLD);

  proc = buffer;
}

template <int NTOTAL, int MINPUT, int NPROCS>
vector<double> Paralel<NTOTAL, MINPUT, NPROCS>::ExchangeWallTemperature(
    const Field &TWall) const {
  PROFILE_FUNCTION();

  int NW = NTOTAL + 2;

  vector<double> recv(NW, 0);
  int wallID = 22;

  if (isLeftToRight())
    PMPI_Sendrecv(&TWall.value[0], NW, MPI_DOUBLE_PRECISION, DCMainProc, wallID,
                  &recv[0], NW, MPI_DOUBLE_PRECISION, DCMainProc, wallID,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  else if (isRightToLeft())
    PMPI_Sendrecv(&TWall.value[0], NW, MPI_DOUBLE_PRECISION, UCMainProc, wallID,
                  &recv[0], NW, MPI_DOUBLE_PRECISION, UCMainProc, wallID,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return recv;
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::ShareWallTemperatureInfo(Field &TWall,
                                                               Field &T) const {
  PROFILE_FUNCTION();

  int n = NI - 1;

  switch (loc) {
    case up:
      SendWallTInTheChannel(TWall);
      for (int i = 1; i < n; i++) T.value[i] = TWall.value[i + iStr];
      break;

    case down:
      SendWallTInTheChannel(TWall);
      for (int i = 1; i < n; i++)
        T.value[i + (NJ - 1) * NI] = TWall.value[i + iStr];
      break;
  }
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::SendWallTInTheChannel(
    Field &TWall) const {
  PROFILE_FUNCTION();

  PMPI_Bcast(&TWall.value[0], NTOTAL + 2, MPI_DOUBLE, 0, myComm.comm());
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::setProcMesh() {
  PROFILE_FUNCTION();

  iStrAllProc = vector<int>(nProcs, 0);
  iEndAllProc = vector<int>(nProcs, 0);
  jStrAllProc = vector<int>(nProcs, 0);
  jEndAllProc = vector<int>(nProcs, 0);

  for (int RowId = 0; RowId < nProcsInRow; RowId++)
    for (int ColId = 0; ColId < nProcsInCol; ColId++) {
      int idN = id(RowId, ColId);
      iStrAllProc[idN] = RowId * (NI - 2);
      iEndAllProc[idN] = (RowId + 1) * (NI - 2);

      jStrAllProc[idN] = ColId * (NJ - 2);
      jEndAllProc[idN] = (ColId + 1) * (NJ - 2);
    }
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::SendInfoToCommMainProc(Field &vec,
                                                             Field &sol) {
  PROFILE_FUNCTION();
  const int iIS = NI - 2;
  const int jIS = NJ - 2;
  const int sendCount = iIS * jIS;

  double sendBuffer[sendCount];

  forAllInterior(NI, NJ) {
    sendBuffer[(i - 1) + (j - 1) * iIS] = vec.value[i + j * NI];
  }

  int recvCount = 0;
  int counts[NPROCS / 2] = {};
  int displacements[NPROCS / 2] = {};

  for (int k = 0; k < nProcs; k++) {
    displacements[k] = recvCount;
    counts[k] = iIS * jIS;
    recvCount += counts[k];
  }
  std::vector<double> recvBuffer(recvCount);
  if (myProc == 0) {
    MPI_Gatherv(&sendBuffer[0], sendCount, MPI_DOUBLE, &recvBuffer[0], counts,
                displacements, MPI_DOUBLE, 0, myComm.comm());
  } else {
    MPI_Gatherv(&sendBuffer, sendCount, MPI_DOUBLE, NULL, NULL, NULL,
                MPI_DOUBLE, 0, myComm.comm());
  }

  if (myProc == 0)
    for (int k = 0; k < nProcs; k++) {
      int iSize = iEndAllProc[k] - iStrAllProc[k];
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++)
        for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++) {
          int iRelative = i - iStrAllProc[k];
          int jRelative = j - jStrAllProc[k];

          int ij = displacements[k] + iRelative + jRelative * iSize;
          sol.value[i + j * NTOTAL] = recvBuffer[ij];
        }
    }

  myComm.wait();
}

template <int NTOTAL, int MINPUT, int NPROCS>
string Paralel<NTOTAL, MINPUT, NPROCS>::PrintMyRank() {
  // string result = to_string();
  return to_string(myProc) + "/" + to_string(nProcs - 1) + "-(" +
         to_string(myRowId) + "/" + to_string(nProcsInRow - 1) + "," +
         to_string(myColId) + "/" + to_string(nProcsInCol - 1) + ")" + "-" +
         to_string(worldMyProc) + "/" + to_string(worldNProcs - 1);
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::PrintProgress(const string &message,
                                                    bool showProgress) {
  if (showProgress == true) {
    MPI_Barrier(MPI_COMM_WORLD);
    cout << message << PrintMyRank() << ".\n";
    if (worldMyProc == 0)
      cout << "------------------------------------------------------------\n";
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::distributeToProcs(Field &sol,
                                                        Field &vec) {
  PROFILE_FUNCTION();

  int counts[NPROCS / 2];
  int displacements[NPROCS / 2];
  int sendCount = 0;

  double buffer[NTOTAL * MINPUT];

  for (int k = 0; k < nProcs; k++) {
    displacements[k] = sendCount;

    for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++)
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++) {
        int relative_i = i - iStrAllProc[k];
        int relative_j = j - jStrAllProc[k];
        int index = displacements[k] +
                    relative_j * (iEndAllProc[k] - iStrAllProc[k]) + relative_i;
        buffer[index] = sol.value[i + j * NTOTAL];
        // if(isWall()&&myProc == 0)
        // std::cout << index << " " << sol[i][j].value << std::endl;
        sendCount++;
      }
    counts[k] = sendCount - displacements[std::max(k, 0)];
    // std::cout << k << "·" << counts[k] << std::endl;
  }
  vector<double> recvValue(counts[myProc]);
  if (myProc == 0)
    MPI_Scatterv(buffer, counts, displacements, MPI_DOUBLE, &recvValue[0],
                 counts[myProc], MPI_DOUBLE, 0, myComm.comm());
  else
    MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, &recvValue[0], counts[myProc],
                 MPI_DOUBLE, 0, myComm.comm());

  for (int i = 0; i < counts[myProc]; i++) {
    int vec_i = i % (NI - 2);
    int vec_j = floor(i / (NI - 2));
    // std::cout << i << " " << vec_i << " " << vec_j << " " << counts[myProc]
    // << " " << recvValue[i]<< endl;
    vec.value[vec_i + 1 + (vec_j + 1) * NI] = recvValue[i];
  }
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::GatherWallTemperature(Field &TWall,
                                                            Field &T) const {
  PROFILE_FUNCTION();
  if (isLeftToRight())
    GatherTemperature(T, TWall, 1);
  else if (isRightToLeft())
    GatherTemperature(T, TWall, T.NJ - 2);
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::GatherTemperature(Field &T, Field &TWall,
                                                        int J) const {
  MPI_Gather(&T.value[1 + J * NI], NI - 2, MPI_DOUBLE, &TWall.value[1], NI - 2,
             MPI_DOUBLE, 0, myComm.comm());
  // for (int i = 0; i < exI1; i++) TWall.value[i] = 0;

  // for (int i = exI2; i < N + 2; i++) TWall.value[i] = 0;
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::scatter(double &from, int size,
                                              double &to, int originProc) {
  MPI_Scatter(&from, size, MPI_DOUBLE, &to, size, MPI_DOUBLE, originProc,
              myComm.comm());
}

template <int NTOTAL, int MINPUT, int NPROCS>
void Paralel<NTOTAL, MINPUT, NPROCS>::setProcWithFixPoint(
    const bool &fixPointInThisProc) {
  if (fixPointInThisProc) fixPointProc += worldMyProc;

  MPI_Allreduce(&fixPointProc, &fixPointProc, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
}
#endif
