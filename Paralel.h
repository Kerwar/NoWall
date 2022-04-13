#ifndef PARALEL_H
#define PARALEL_H

#include <cmath>
#include <vector>
#include <string>
#include "Field.h"
#include "mpi.h"

using std::vector;

class Paralel
{
public:
  Paralel(int worldsize);
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

  int id(int, int);
  int worldid(int, int);
  void setUpComm(int &NX, int &NY, double &xMax, double &xExMin, double &xExMax);
  void setUpMesh(int &NX, int &NY, double &channelWidth, int &exi1, int &exi2);
  void SendInfoToNeighbours(Field::vectorField &vec);
  void SendInfoToCommMainProc(Field::vectorField &vec, Field::vectorField &sol);
  void ExchangeWallTemperature(Field::vec1dfield &TWall, Field::vec1dfield &TNextToWall, double &exCte);
  void GatherWallTemperature(Field::vec1dfield &TWall, Field::vec1dfield &TNextToWall, Field::vectorField &T, double &exCte);
  void ShareWallTemperatureInfo(Field::vec1dfield &TWall, Field::vectorField &T);
  void SendWallTInTheChannel(Field::vec1dfield &TWall, vector<double> &TW);
  void MainsProcs(Loc location, int &proc);
  void freeComm();
  void PrintProgress(string message, bool showProgress);
  void distributeToProcs(Field::vectorField &sol, Field::vectorField &vec);

  inline bool isLeftToRight() { return loc == up ? true : false; };
  inline bool isRightToLeft() { return loc == down ? true : false; };

  inline bool isProcNull(int &proc) { return proc == MPI_PROC_NULL ? true : false; };
  inline bool isFirstProcOfRow() { return myRowId == 0 ? true : false; };
  inline bool isLastProcOfRow() { return myRowId == nProcsInRow - 1 ? true : false; };
  inline bool isFirstProcOfCol() { return myColId == 0 ? true : false; };
  inline bool isLastProcOfCol() { return myColId == nProcsInCol - 1 ? true : false; };
  string PrintMyRank();

private:
  void setProcMesh(int &NX, int &NY);
  double TWUinBorder(double (&T)[4], double &exCte);
  double TWDinBorder(double (&T)[4], double &exCte);
  void GatherTemperature(Field::vectorField &T, Field::vec1dfield &TWall, int J);
  int SetXProcs(double &nmRelation);
  Field::vec1dfield TIntoDerivative(const Field::vec1dfield &TWall,const Field::vec1dfield &TW2R, double exCte, double ex2);
};

#endif
