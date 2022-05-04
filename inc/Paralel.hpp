#ifndef PARALEL_H
#define PARALEL_H

#include <cmath>
#include <vector>
#include <string>
#include "Field.hpp"
#include "mpi.h"

using std::vector;

class Paralel
{
public:
  Paralel();
  virtual ~Paralel();

  enum Loc
  {
    up,
    down,
  };

  int worldNProcs, worldMyProc;

  int channelURank;
  int channelDRank;
  int cUSize;
  int cDSize;

  int nProcsInRow, nProcsInCol;
  int nProcs, myProc;
  int myRowId, myColId;
  int myLeft, myRight, myBot, myTop;
  Loc loc;
  int UCMainProc, DCMainProc;

  int NYChannel, NXWall;
  int locIStr, locIEnd, locJStr, locJEnd;
  int exI1, exI2;
  int myNx, myNy;
  int iStr, iEnd;
  int jStr, jEnd;

  int nProcsExtraXPoint, nProcsExtraYPoint;

  vector<int> iStrAllProc, iEndAllProc;
  vector<int> jStrAllProc, jEndAllProc;

  int fixPointProc;
  MPI_Comm myComm;
  MPI_Datatype column_type;

  int id(int, int);
  int worldid(int, int);
  void setUpComm(const int &NX, const int &NY);
  void setUpMesh(const int &NX, const int &NY, int &exi1, int &exi2);
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
  void setProcMesh(const int &NX, const int &NY);
  double TWUinBorder(double (&T)[4], double &exCte);
  double TWDinBorder(double (&T)[4], double &exCte);
  void GatherTemperature(Field &T, Field &TWall, int J);
  int SetXProcs(double &nmRelation);
  // Field::vec1dfield TIntoDerivative(const Field &TWall,const Field &TW2R, double exCte, double ex2);
};

#endif
