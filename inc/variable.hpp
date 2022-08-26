#ifndef VARIABLE_H
#define VARIABLE_H

#include "equation.hpp"
#include "field.hpp"
#include "filereader.hpp"
#include "finitevolumeoperations.hpp"
#include "grid.hpp"
#include "paralel.hpp"

using namespace parameters;

class Variable {
 public:
  Variable();
  virtual ~Variable();

  void passInfoSolutionGrid(const Grid &solGrid, double viscX, double viscY);
  void passInfoMyGrid(const Grid &myGrid, double viscX, double viscY);

  void initializeWall();
  void readFile(Paralel<NTOTAL, MINPUT, NPROCS> &paralel, int block, double &m,
                double yMin, double yMax);

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

  inline void sendInfoToCommMainProc(Paralel<NTOTAL, MINPUT, NPROCS> &paralel) {
    PROFILE_FUNCTION();
    paralel.SendInfoToCommMainProc(T, solT);
    paralel.SendInfoToCommMainProc(F, solF);
    paralel.SendInfoToCommMainProc(Z, solZ);
    paralel.SendInfoToCommMainProc(U, solU);
    paralel.SendInfoToCommMainProc(V, solV);
  }
  inline void sendInfoToNeighbours(Paralel<NTOTAL, MINPUT, NPROCS> &paralel) {
    PROFILE_FUNCTION();
    paralel.SendInfoToNeighbours(T);
    paralel.SendInfoToNeighbours(F);
    paralel.SendInfoToNeighbours(Z);
  }

  inline void exchangeTemperature(
      const Paralel<NTOTAL, MINPUT, NPROCS> &paralel, const Grid &mainGrid) {
    PROFILE_FUNCTION();

    paralel.GatherWallTemperature(TWall, T);

    for (int i = 0; i < mainGrid.exI1; i++) TWall.value[i] = 0;
    for (int i = mainGrid.exI2; i < NTOTAL + 2; i++) TWall.value[i] = 0;

    if (paralel.myProc == 0)
      TWall.value = paralel.ExchangeWallTemperature(TWall);
    MPI_Barrier(MPI_COMM_WORLD);

    paralel.ShareWallTemperatureInfo(TWall, T);

    // paralel.scatter(TWall.value[1], NI - 2, DTinWall[1], 0);
    MPI_Scatter(&TWall.value[1], NI - 2, MPI_DOUBLE, &DTinWall[1], NI - 2,
                MPI_DOUBLE, 0, paralel.myComm.comm());
  }

  void setChannelEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn,
                           const Paralel<NTOTAL, MINPUT, NPROCS> &paralel,
                           const double &m, const double &q, double &DT,
                           const double &exCte, const int &iter);
  void setWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn,
                    Direction side);
  void setExchangeWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn,
                            Direction side, int iStr, int iEnd, int mainexI1,
                            int mainexI2, int exI1, int exI2);
  inline void setDirichlet(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn,
                           Direction side) {
    PROFILE_FUNCTION();
    Teqn->SetDirichlet(T, side);
    Feqn->SetDirichlet(F, side);
    Zeqn->SetDirichlet(Z, side);
  }

  inline void assembleEquations(Equation *&Teqn, Equation *&Feqn,
                                Equation *&Zeqn) {
    PROFILE_FUNCTION();

    Teqn->assembleEquation();
    Feqn->assembleEquation();
    Zeqn->assembleEquation();

    Teqn->relax(T);
    Feqn->relax(F);
    Zeqn->relax(Z);
  }
  double solveEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn,
                        double alpha, int &niter, int &itersol, int changeIter);
  double calculateNewM(const double &alpha, const double &m, const double &q);
  double calculateNewQ(const double &alpha, const double &m, const double &q);
  void setFixIndex(double xfix, double yfix);
  inline void updateBoundFactors() {
    PROFILE_FUNCTION();

    lowerBoundFactorm *= 0.9;
    upperBoundFactorm *= 1.1;
  }
  bool isTOutEqualToQ(double q, const Paralel<NTOTAL, MINPUT, NPROCS> &paralel);
  void writeTInWall(Paralel<NTOTAL, MINPUT, NPROCS> &paralel,
                    const Grid &mainGrid, const Grid &myGrid, int iter);

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
