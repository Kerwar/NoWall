#ifndef PROBLEM_H
#define PROBLEM_H

#include "Paralel.hpp"
#include "Grid.hpp"
#include "Variable.hpp"
#include "FileWriter.hpp"
#include "Equation.hpp"

class Problem
{
public:
  Problem();
  Problem(const int &nx, const int &ny,
   const int &ni, const int &nj, const double &prevm);
  virtual ~Problem();

  void setUpProblem();
  void initializeVariables();
  void writeSolution(string &filename, int i);
  void freeComm();
  double mainIter(int i);
  void setExchangeConstant(double &deltaY);
  void writefilename(string &filename);
  void retrieveNandM(int &nxOut, int &nyOut, double &mOut);
  bool fixPointInThisProc();
  inline bool isSolutionNotGoodEnough(){return !variables.isTOutEqualToQ(q, paralel);};

  int maxit, iShow;
  double m;
  const int N, M;
  Paralel paralel;

private:
  static double tolerance;
  static double alpha;

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

  static bool readFromFile;
  Grid mainGrid, myGrid;
  Variable variables;

  int myProc, nProcs;
  double myXMin, myXMax;
  double myYMin, myYMax;

  int solN, solM;
  
  Equation *Teqn;
  Equation *Feqn;
  Equation *Zeqn;
  
  string sufix;
public:
  inline bool isErrorSmallEnough(double &error){return error < tolerance ? true : false;};
};

#endif