#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <mpi.h>

class Processor
{
public:
  Processor(MPI_Comm myComm);
  virtual ~Processor();

  int id;
  int size;

  int worldId;
  int worldSize;

  int myRowId, myColId;
  int nProcsInRow, nProcsInCol;
  int myLeft, myRight, myBot, myTop;

  int myNx, myNy;
  int iStr, iEnd;
  int jStr, jEnd;
};

#endif