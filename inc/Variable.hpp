#ifndef VARIABLE_H
#define VARIABLE_H

#include "Field.hpp"
#include "Grid.hpp"
#include "Paralel.hpp"
#include "Equation.hpp"
#include "FiniteVolumeOperations.hpp"
#include "FileReader.hpp"

class Variable
{
public:
  // Variable();
  Variable(const int &n,const int &m,const int &soln,const int &solm);
  virtual ~Variable();

  void passInfoGridToAll(const Grid &solGrid, const Grid &myGrid, double &viscX, double &viscY, double &LeF, double &LeZ);
  void initializeLeftToRight(double m, double yMin, double yMax, double q, double xHS, double r0hs, double z0hs, double xMin, double xMax);
  void initializeRightToLeft(double m, double yMin, double yMax, double q, double xHS, double r0hs, double z0hs, double xMin, double xMax);
  void initializeWall();
  void readFile(Paralel &paralel, int block, double &m, double yMin, double yMax);
  void setMassFluxes(const Grid &myGrid);
  void setInletBoundaryConditionLeftToRight();
  void setInletBoundaryConditionRightToLeft();
  void sendInfoToCommMainProc(Paralel &paralel);
  void sendInfoToNeighbours(Paralel &paralel);
  void exchangeTemperature(Paralel &paralel, double &exCte, int &solExI1, int &solExI2);
  void setWallEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double &alphaWall, double &DT);
  void setChannelEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double &m, double &q, double &beta, double &gamma, double &DT);
  void setWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side);
  void setExchangeWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side, int iStr, int iEnd, int mainexI1, int mainexI2, int exI1, int exI2);
  void setDirichlet(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side);
  void assembleEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn);
  double solveEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double alpha, int &niter, int &itersol, int changeIter);
  double calculateNewM(double alpha, double m, double q);
  void setFixIndex(double xfix, double yfix);
  void updateBoundFactors();
  bool isTOutEqualToQ(double q, Paralel &paralel);
  void writeTInWall(Paralel &paralel, const Grid &mainGrid, const Grid &myGrid, int iter);

  const int NI, NJ;
  const int solN, solM;
  Field solT;
  Field solF, solZ;
  Field solU, solV;
  bool manyIter = false;

private:

  int ifix, jfix;
  Field U, V;

  Field T;
  Field F, Z;

  Field TWall, TNextToWall;

  Field massFluxE, massFluxN;
  double lowerBoundFactorm = 0.9, upperBoundFactorm = 1.1;

  string initialsol = "Sol-1.f";
};
#endif
