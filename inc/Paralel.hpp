#ifndef PARALEL_H
#define PARALEL_H

#include <cmath>
#include <vector>
#include <string>
#include "Field.hpp"
#include "mpi.h"

using std::cout;
using std::endl;
using std::to_string;

using std::vector;

enum Loc
{
  up,
  down,
};

template <int N, int M, int NPROCS>
struct Paralel
{
public:
  Paralel();
  virtual ~Paralel();

  static constexpr int nProcsInRow = NPROCS / 2;
  static constexpr int nProcsInCol = 1;
  static constexpr int NI = N / nProcsInRow + 2;
  static constexpr int NJ = M / 2 + 2;

  int worldNProcs, worldMyProc;

  int channelURank;
  int channelDRank;

  int nProcs, myProc;
  int myRowId, myColId;
  int myLeft, myRight, myBot, myTop;
  Loc loc;
  int UCMainProc, DCMainProc;

  int NYChannel, NXWall;
  int locIStr, locIEnd, locJStr, locJEnd;
  int exI1, exI2;
  int iStr, iEnd;
  int jStr, jEnd;

  int nProcsExtraXPoint, nProcsExtraYPoint;

  vector<int> iStrAllProc, iEndAllProc;
  vector<int> jStrAllProc, jEndAllProc;

  int fixPointProc;
  MPI_Comm myComm;
  MPI_Datatype column_type;

  int id(int, int);
  inline int worldid(int row, int col) { return col + row * nProcsInCol; };
  void setUpMesh(int &exi1, int &exi2);
  void SendInfoToNeighbours(Field &vec);
  void SendInfoToCommMainProc(Field &vec, Field &sol);
  void ExchangeWallTemperature(Field &TWall, Field &TNextToWall, double &exCte, int &solExI1, int &solExI2);
  void GatherWallTemperature(Field &TWall, Field &TNextToWall, Field &T);
  void ShareWallTemperatureInfo(Field &TWall, Field &T);
  void SendWallTInTheChannel(Field &TWall);
  void MainsProcs(Loc location, int &proc);
  void freeComm();
  void PrintProgress(string message, bool showProgress);
  void distributeToProcs(Field &sol, Field &vec);

  inline bool isLeftToRight() { return loc == up ? true : false; };
  inline bool isRightToLeft() { return loc == down ? true : false; };

  inline bool isProcNull(int &proc) { return proc == MPI_PROC_NULL ? true : false; };
  inline bool isFirstProcOfRow() { return myRowId == 0 ? true : false; };
  inline bool isLastProcOfRow() { return myRowId == nProcsInRow - 1 ? true : false; };
  inline bool isFirstProcOfCol() { return myColId == 0 ? true : false; };
  inline bool isLastProcOfCol() { return myColId == nProcsInCol - 1 ? true : false; };
  string PrintMyRank();

private:
  void setProcMesh();
  void GatherTemperature(Field &T, Field &TWall, int J);
  // Field::vec1dfield TIntoDerivative(const Field &TWall,const Field &TW2R, double exCte, double ex2);
};

template <int N, int M, int NPROCS>
Paralel<N, M, NPROCS>::Paralel()
{
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

  MPI_Cart_create(helpComm, 2, dims, period, reorder, &myComm);
  MPI_Comm_free(&helpComm);

  MPI_Comm_rank(myComm, &myProc);

  int myCoords[2];
  MPI_Cart_coords(myComm, myProc, 2, myCoords);

  myRowId = myCoords[0];
  myColId = myCoords[1];

  MPI_Cart_shift(myComm, 0, 1, &myLeft, &myRight);
  MPI_Cart_shift(myComm, 1, 1, &myBot, &myTop);

  MainsProcs(up, UCMainProc);
  MainsProcs(down, DCMainProc);
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::setUpMesh(int &exi1, int &exi2)
{
  PROFILE_FUNCTION();
  exI1 = exi1;
  exI2 = exi2;

  NXWall = exI2 - exI1;

  switch (loc)
  {
  case up:
    setProcMesh();

    locIStr = 0;
    locIEnd = N;
    locJStr = NJ - 2;
    locJEnd = M - 1;

    for (int i = 0; i < NPROCS / 2; i++)
    {
      jStrAllProc[i] += locJStr;
      jEndAllProc[i] += locJStr;
    }
    break;

  case down:
    setProcMesh();

    locIStr = 0;
    locIEnd = N;
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

template <int N, int M, int NPROCS>
Paralel<N, M, NPROCS>::~Paralel()
{
  MPI_Type_free(&column_type);
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::freeComm()
{
  if (myComm != MPI_COMM_NULL)
    MPI_Comm_free(&myComm);
}

template <int N, int M, int NPROCS>
int Paralel<N, M, NPROCS>::id(int x, int y)
{
  int myCoords[2] = {x, y};
  int result;

  MPI_Cart_rank(myComm, myCoords, &result);

  return result;
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::SendInfoToNeighbours(Field &vec)
{
  int leftTag = 10, rightTag = 11, topTag = 12, botTag = 13;

  leftTag *= (1 + loc);
  rightTag *= (1 + loc);
  topTag *= (1 + loc);
  botTag *= (1 + loc);

  MPI_Barrier(myComm);
  // SENDING INFO TO THE LEFT
  MPI_Sendrecv(&vec.value[1], 1, column_type, myLeft, leftTag,
               &vec.value[NI - 1], 1, column_type, myRight, leftTag, myComm, MPI_STATUS_IGNORE);

  if (!isProcNull(myRight))
    for (int j = 0; j < NJ; j++)
      vec.value[NI - 1 + j * NI] = (vec.value[NI - 1 + j * NI] + vec.value[NI - 2 + j * NI]) * 0.5;

  MPI_Barrier(myComm);
  // SENDING INFO TO THE RIGHT

  MPI_Sendrecv(&vec.value[NI - 1], 1, column_type, myRight, rightTag,
               &vec.value[0], 1, column_type, myLeft, rightTag, myComm, MPI_STATUS_IGNORE);

  MPI_Barrier(myComm);

  // SENDING INFO TO THE BOT
  MPI_Sendrecv(&vec.value[NI], NI, MPI_DOUBLE_PRECISION, myBot, botTag,
               &vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop, botTag, myComm, MPI_STATUS_IGNORE);

  if (!isProcNull(myTop))
    for (int i = 0; i < NI; i++)
      vec.value[i + (NJ - 1) * NI] = (vec.value[i + (NJ - 1) * NI] + vec.value[i + (NJ - 2) * NI]) * 0.5;
  MPI_Barrier(myComm);

  // SENDING INFO TO THE TOP
  MPI_Sendrecv(&vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop, topTag,
               &vec.value[0], NI, MPI_DOUBLE_PRECISION, myBot, topTag, myComm, MPI_STATUS_IGNORE);
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::MainsProcs(Loc location, int &proc)
{
  PROFILE_FUNCTION();
  int buffer;

  if (location == loc)
  {
    if (myProc == 0)
      buffer = worldMyProc;
    PMPI_Bcast(&buffer, 1, MPI_INT, 0, myComm);
  }

  int sendProc;
  switch (location)
  {
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

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::ExchangeWallTemperature(Field &TWall, Field &TNextToWall, double &exCte, int &solExI1, int &solExI2)
{
  PROFILE_FUNCTION();

  int NW = N + 2;
  double wallRecv[NW] = {};
  double nextToWallRecv[NW] = {};

  int wallID = 22;
  int nextToWallID = 23;

  if (isLeftToRight())
  {
    PMPI_Sendrecv(&TWall.value[0], NW, MPI_DOUBLE_PRECISION, DCMainProc, wallID, &wallRecv, NW, MPI_DOUBLE_PRECISION, DCMainProc, wallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    PMPI_Sendrecv(&TNextToWall.value[0], NW, MPI_DOUBLE_PRECISION, DCMainProc, nextToWallID, &nextToWallRecv, NW, MPI_DOUBLE_PRECISION, DCMainProc, nextToWallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = solExI1; i < solExI2; i++)
      TWall.value[i] = ((8.0 + exCte) * (9.0 * TWall.value[i] - TNextToWall.value[i]) + 9.0 * exCte * wallRecv[i] - exCte * nextToWallRecv[i]) / (64.0 + 16 * exCte);
  }
  else if (isRightToLeft())
  {
    PMPI_Sendrecv(&TWall.value[0], NW, MPI_DOUBLE_PRECISION, UCMainProc, wallID, &wallRecv, NW, MPI_DOUBLE_PRECISION, UCMainProc, wallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    PMPI_Sendrecv(&TNextToWall.value[0], NW, MPI_DOUBLE_PRECISION, UCMainProc, nextToWallID, &nextToWallRecv, NW, MPI_DOUBLE_PRECISION, UCMainProc, nextToWallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = solExI1; i < solExI2; i++)
      TWall.value[i] = ((8.0 + exCte) * (9.0 * TWall.value[i] - TNextToWall.value[i]) + 9.0 * exCte * wallRecv[i] - exCte * nextToWallRecv[i]) / (64.0 + 16 * exCte);
  }
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::ShareWallTemperatureInfo(Field &TWall, Field &T)
{
  PROFILE_FUNCTION();

  int n = NI - 1;

  switch (loc)
  {
  case up:
    SendWallTInTheChannel(TWall);
    for (int i = 1; i < n; i++)
      T.value[i] = TWall.value[i + iStr];
    break;

  case down:
    SendWallTInTheChannel(TWall);
    for (int i = 1; i < n; i++)
      T.value[i + (NJ - 1) * NI] = TWall.value[i + iStr];
    break;
  }
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::SendWallTInTheChannel(Field &TWall)
{
  PROFILE_FUNCTION();

  PMPI_Bcast(&TWall.value[0], N + 2, MPI_DOUBLE, 0, myComm);
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::setProcMesh()
{
  PROFILE_FUNCTION();

  iStrAllProc = vector<int>(nProcs, 0);
  iEndAllProc = vector<int>(nProcs, 0);
  jStrAllProc = vector<int>(nProcs, 0);
  jEndAllProc = vector<int>(nProcs, 0);

  for (int RowId = 0; RowId < nProcsInRow; RowId++)
    for (int ColId = 0; ColId < nProcsInCol; ColId++)
    {
      int idN = id(RowId, ColId);
      iStrAllProc[idN] = RowId * (NI - 2);
      iEndAllProc[idN] = (RowId + 1) * (NI - 2);

      jStrAllProc[idN] = ColId * (NJ - 2);
      jEndAllProc[idN] = (ColId + 1) * (NJ - 2);
    }
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::SendInfoToCommMainProc(Field &vec, Field &sol)
{
  PROFILE_FUNCTION();
  int iIS = NI - 2;
  int jIS = NJ - 2;
  int sendCount = iIS * jIS;

  double sendBuffer[sendCount] = {};

  forAllInterior(NI, NJ)
  {
    sendBuffer[(i - 1) + (j - 1) * iIS] = vec.value[i + j * NI];
  }

  int recvCount = 0;
  int counts[nProcs] = {};
  int displacements[nProcs] = {};

  for (int k = 0; k < nProcs; k++)
  {
    displacements[k] = recvCount;
    counts[k] = iIS * jIS;
    recvCount += counts[k];
  }
  double recvBuffer[recvCount];
  if (myProc == 0)
  {
    MPI_Gatherv(&sendBuffer[0], sendCount, MPI_DOUBLE, recvBuffer, counts, displacements,
                MPI_DOUBLE, 0, myComm);
  }
  else
  {
    MPI_Gatherv(&sendBuffer, sendCount, MPI_DOUBLE, NULL, NULL, NULL,
                MPI_DOUBLE, 0, myComm);
  }

  if (myProc == 0)
    for (int k = 0; k < nProcs; k++)
    {
      int iSize = iEndAllProc[k] - iStrAllProc[k];
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++)
        for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++)
        {
          int iRelative = i - iStrAllProc[k];
          int jRelative = j - jStrAllProc[k];

          int ij = displacements[k] + iRelative + jRelative * iSize;
          sol.value[i + j * N] = recvBuffer[ij];
        }
    }

  MPI_Barrier(myComm);
}

template <int N, int M, int NPROCS>
string Paralel<N, M, NPROCS>::PrintMyRank()
{
  // string result = to_string();
  return to_string(myProc) + "/" + to_string(nProcs - 1) +
         "-(" + to_string(myRowId) + "/" + to_string(nProcsInRow - 1) +
         "," + to_string(myColId) + "/" + to_string(nProcsInCol - 1) + ")" +
         "-" + to_string(worldMyProc) + "/" + to_string(worldNProcs - 1);
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::PrintProgress(string message, bool showProgress)
{
  if (showProgress == true)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    cout << message << PrintMyRank() << ".\n";
    if (worldMyProc == 0)
      cout << "------------------------------------------------------------\n";
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::distributeToProcs(Field &sol, Field &vec)
{
  PROFILE_FUNCTION();

  int counts[nProcs];
  int sendCount = 0;
  int displacements[nProcs];

  double buffer[N * M];

  for (int k = 0; k < nProcs; k++)
  {
    displacements[k] = sendCount;

    for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++)
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++)
      {
        int relative_i = i - iStrAllProc[k];
        int relative_j = j - jStrAllProc[k];
        int index = displacements[k] + relative_j * (iEndAllProc[k] - iStrAllProc[k]) + relative_i;
        buffer[index] = sol.value[i + j * N];
        // if(isWall()&&myProc == 0)
        // std::cout << index << " " << sol[i][j].value << std::endl;
        sendCount++;
      }
    counts[k] = sendCount - displacements[std::max(k, 0)];
    // std::cout << k << "·" << counts[k] << std::endl;
  }
  double recvValue[counts[myProc]];
  if (myProc == 0)
    MPI_Scatterv(buffer, counts, displacements, MPI_DOUBLE, &recvValue, counts[myProc], MPI_DOUBLE, 0, myComm);
  else
    MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvValue, counts[myProc], MPI_DOUBLE, 0, myComm);

  for (int i = 0; i < counts[myProc]; i++)
  {
    int vec_i = i % (NI - 2);
    int vec_j = floor(i / (NI - 2));
    // std::cout << i << " " << vec_i << " " << vec_j << " " << counts[myProc] << " " << recvValue[i]<< endl;
    vec.value[vec_i + 1 + (vec_j + 1) * NI] = recvValue[i];
  }
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::GatherWallTemperature(Field &TWall, Field &TNextToWall, Field &T)
{
  PROFILE_FUNCTION();
  if (isLeftToRight())
  {
    GatherTemperature(T, TNextToWall, 2);
    GatherTemperature(T, TWall, 1);
  }
  else if (isRightToLeft())
  {
    GatherTemperature(T, TNextToWall, T.NJ - 3);
    GatherTemperature(T, TWall, T.NJ - 2);
  }
}

template <int N, int M, int NPROCS>
void Paralel<N, M, NPROCS>::GatherTemperature(Field &T, Field &TWall, int J)
{
  int count[nProcs];
  int displacement[nProcs];

  for (int i = 0; i < nProcs; i++)
  {
    count[i] = iEndAllProc[i] - iStrAllProc[i];
    displacement[i] = iStrAllProc[i];
  }

  if (myProc == 0)
    MPI_Gatherv(&T.value[1 + J * NI], NI - 2, MPI_DOUBLE, &TWall.value[1], count, displacement, MPI_DOUBLE, 0, myComm);
  else
    MPI_Gatherv(&T.value[1 + J * NI], NI - 2, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, myComm);

  for (int i = 0; i < exI1; i++)
    TWall.value[i] = 0;

  for (int i = exI2; i < N + 2; i++)
    TWall.value[i] = 0;
}

#endif
