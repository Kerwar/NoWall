#ifndef PROBLEM_H
#define PROBLEM_H

#include "Paralel.h" 
#include "Grid.h"
#include "Variable.h"
#include "FileWriter.h"
#include "Equation.h"

class Problem
{
public:
  Problem(int worldsize);
  Problem(int worldsize, int nx, int my, double prevm);
  virtual ~Problem();

  void setUpProblem();
  void initializeVariables();
  void writeSolution(string &filename, int i);
  void freeComm();
  double mainIter(int i);
  void setExchangeConstant(double &deltaY);
  void writefilename(string &filename);
  bool isErrorSmallEnough(double error);
  bool fixPointInThisProc();
  void retrieveNandM(int &nxOut, int &nyOut, double &mOut);
  bool isSolutionNotGoodEnough();

  int maxit, iShow;
  double m;
private:

  int N, M;
  double tolerance;
  Paralel paralel;
  double alpha;

  static double beta, gamma;
  static double xMin, xMax, xExMin, xExMax;
  static double xfix, yfix;
  static double yWall, yChannel, yMin, yMax;
  static double q, xHS_U, xHS_D, r0hs, z0hs;

  static double alphaWall, DT;

  int itersol;
  static double LeF, LeZ;

  static double a, a2;
  static double bK, exCte;

  static double viscX, viscY;

  Grid mainGrid, myGrid;
  Variable *variables;

  int myProc, nProcs;
  double myXMin, myXMax;
  double myYMin, myYMax;

  int solN, solM;
  bool readFromFile;

  string sufix;
};

#endif