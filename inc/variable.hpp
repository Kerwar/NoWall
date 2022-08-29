#ifndef VARIABLE_H
#define VARIABLE_H

#include "field.hpp"
#include "filereader.hpp"
#include "finitevolumeoperations.hpp"
#include "grid.hpp"
#include "paralel.hpp"
#include "systemofequations.hpp"

using namespace parameters;

class Variable {
 public:
  Variable();
  virtual ~Variable();

  void passInfoSolutionGrid(const Grid &solGrid, double viscX, double viscY);
  void passInfoMyGrid(const Grid &myGrid, double viscX, double viscY);

  void initializeWall();
  void readFile(Paralel &paralel, int block, double &m, double yMin,
                double yMax);

  inline void setMassFluxes(const Grid &myGrid) {
    PROFILE_FUNCTION();

    massFluxE.getGridInfoPassed(myGrid, U.viscX[0], U.viscY[0]);
    massFluxN.getGridInfoPassed(myGrid, V.viscX[0], V.viscY[0]);

    massFluxE.computeEastMassFluxes(U);
    massFluxN.computeNorthMassFluxes(V);
  }

  inline void setInletBoundaryConditionLeftToRight() {
    PROFILE_FUNCTION();
    T.inletBoundaryCondition(west, 0.0);
    F.inletBoundaryCondition(west, 1.0);
    Z.inletBoundaryCondition(west, 0.0);
  }
  inline void setInletBoundaryConditionRightToLeft() {
    PROFILE_FUNCTION();
    T.inletBoundaryCondition(east, 0.0);
    F.inletBoundaryCondition(east, 1.0);
    Z.inletBoundaryCondition(east, 0.0);
  }

  inline void sendInfoToCommMainProc(Paralel &paralel) {
    PROFILE_FUNCTION();
    paralel.SendInfoToCommMainProc(T, solT);
    paralel.SendInfoToCommMainProc(F, solF);
    paralel.SendInfoToCommMainProc(Z, solZ);
    paralel.SendInfoToCommMainProc(U, solU);
    paralel.SendInfoToCommMainProc(V, solV);
  }
  inline void sendInfoToNeighbours(Paralel &paralel) {
    PROFILE_FUNCTION();
    paralel.SendInfoToNeighbours(T);
    paralel.SendInfoToNeighbours(F);
    paralel.SendInfoToNeighbours(Z);
  }

  void exchangeTemperature(const Paralel &paralel, const Grid &mainGrid);

  void setChannelEquations(SystemOfEquations &sys, const Paralel &paralel,
                           const double &m, const double &q, double &DT,
                           const double &exCte, const int &iter);
  void setWallShear(SystemOfEquations &sys, Direction side);
  void setExchangeWallShear(SystemOfEquations &sys, Direction side, int iStr,
                            int iEnd, int mainexI1, int mainexI2, int exI1,
                            int exI2);
  inline void setDirichlet(SystemOfEquations &sys, Direction side) {
    PROFILE_FUNCTION();
    sys.T->SetDirichlet(T, side);
    sys.F->SetDirichlet(F, side);
    sys.Z->SetDirichlet(Z, side);
  }

  inline void assembleEquations(SystemOfEquations &sys) {
    PROFILE_FUNCTION();

    sys.T->assembleEquation();
    sys.F->assembleEquation();
    sys.Z->assembleEquation();

    sys.T->relax(T);
    sys.F->relax(F);
    sys.Z->relax(Z);
  }
  double solveEquations(SystemOfEquations &sys, double alpha, int &niter,
                        int &itersol, int changeIter);
  double calculateNewM(const double &alpha, const double &m, const double &q);
  double calculateNewQ(const double &alpha, const double &m, const double &q);
  void setFixIndex(double xfix, double yfix);
  inline void updateBoundFactors() {
    PROFILE_FUNCTION();

    lowerBoundFactorm *= 0.9;
    upperBoundFactorm *= 1.1;
  }
  bool isTOutEqualToQ(double q, const Paralel &paralel);

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

  vector<double> DTinWall;
  Field massFluxE, massFluxN;
  double lowerBoundFactorm = 0.9, upperBoundFactorm = 1.1;
  double DY;

  string initialsol = "Sol-1.f";

 public:
  inline void initializeLeftToRight(const double &m, const double &yMin,
                                    const double &yMax, const double &q,
                                    const double &xHS, const double &r0hs,
                                    const double &z0hs, const double &xMin,
                                    const double &xMax) {
    PROFILE_FUNCTION();
    U.laminarFlow(m, yMin, yMax);
    V.initializeInternalField(0);
    T.InitializeT(q, xHS, xMin, xMax);
    F.InitializeF(xHS, xMin, xMax);
    Z.InitializeZ(z0hs, r0hs, xHS, (yMin + yMax) / 2.0);
  }
  inline void initializeRightToLeft(const double &m, const double &yMin,
                                    const double &yMax, const double &q,
                                    const double &xHS, const double &r0hs,
                                    const double &z0hs, const double &xMin,
                                    const double &xMax) {
    PROFILE_FUNCTION();
    U.laminarFlow(-m, yMin, yMax);
    V.initializeInternalField(0);
    T.InitializeT(q, xHS, xMax, xMin);
    F.InitializeF(xHS, xMax, xMin);
    Z.InitializeZ(z0hs, r0hs, xHS, (yMin + yMax) / 2.0);
  }
};

#endif
