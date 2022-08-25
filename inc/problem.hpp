#ifndef PROBLEM_H
#define PROBLEM_H

#include <memory>

#include "parameters.hpp"
#include "Equation.hpp"
#include "FileWriter.hpp"
#include "Grid.hpp"
#include "Paralel.hpp"
#include "Variable.hpp"

using namespace parameters;

class Problem {

 public:
  Problem(const double &prevm);
  virtual ~Problem();

  void setUpProblem();
  void initializeVariables();
  void writeSolution(string &filename, int i);
  void freeComm();
  double mainIter(int i);
  void writefilename(string &filename);
  void retrieveNandM(int &nxOut, int &nyOut, double &mOut);
  bool isSolutionNotGoodEnough() {
    return !variables.isTOutEqualToQ(q, paralel);
  };

  int maxit, iShow;
  double m,q;
  double xfix;
  Paralel<NTOTAL, MINPUT, NPROCS> paralel;

  inline bool isErrorSmallEnough(double &error) {
    return error < TOL ? true : false;
  };
 private:
  double alpha;

  double yfix;
  double yWall, yChannel, yMin, yMax;
  double xHS_U, xHS_D, z0hs;

  double alphaWall, DT;

  int itersol;

  double bK, exCte;

  double viscX, viscY;

  bool readFromFile;
  Grid mainGrid, myGrid;
  Variable<NTOTAL, MINPUT, NTOTAL+2, M+2, NPROCS> variables;

  int myProc, nProcs;
  double myXMin, myXMax;
  double myYMin, myYMax;

  Equation *Teqn;
  Equation *Feqn;
  Equation *Zeqn;

  string sufix;
  
  bool fixPointInThisProc();
  void setExchangeConstant();
  bool isFixPointXCoordinateInThisProc();
  bool isFixPointYCoordinateInThisProc();
};

#endif
