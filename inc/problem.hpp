#ifndef PROBLEM_H
#define PROBLEM_H

#include <memory>

#include "filewriter.hpp"
#include "variable.hpp"

using namespace parameters;

class Problem {
 public:
  Problem(const double &prevm);
  virtual ~Problem();

  void setUpProblem();
  void initializeVariables();
  void writeSolution(string &filename, int i);
  double mainIter(int i);
  void writefilename(string &filename);
  void retrieveNandM(int &nxOut, int &nyOut, double &mOut);
  bool isSolutionNotGoodEnough() {
    return !variables.isTOutEqualToQ(q, paralel);
  };

  int maxit, iShow;
  double m, q;
  double xfix;
  Paralel paralel;

  inline bool isErrorSmallEnough(double &error) {
    return error < TOL ? true : false;
  };

 private:
  double alpha;

  double yfix;
  double xHS_U, xHS_D, z0hs;

  double alphaWall, DT;

  int itersol;

  double bK, exCte;

  double viscX, viscY;

  bool readFromFile;
  Grid mainGrid, myGrid;
  VariableManager variables;

  int myProc, nProcs;
  double myXMin, myXMax;
  double myYMin, myYMax;

  SystemOfEquations sys;

  string sufix;

  bool fixPointInThisProc();
  void setExchangeConstant();
  bool isFixPointXCoordinateInThisProc();
  bool isFixPointYCoordinateInThisProc();
};

#endif
