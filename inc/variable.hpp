#ifndef _VARIABLE_MANAGER_HPP
#define _VARIABLE_MANAGER_HPP

#include "field.hpp"
#include "filereader.hpp"
#include "finitevolumeoperations.hpp"
#include "grid.hpp"
#include "paralel.hpp"
#include "systemofequations.hpp"

using namespace parameters;

class VariableManager {
 public:
  VariableManager();
  virtual ~VariableManager();

  void passInfoSolutionGrid(const Grid &solGrid, double viscX, double viscY);
  void passInfoMyGrid(const Grid &myGrid, double viscX, double viscY);

  void initializeWall();
  void readFile(Paralel &paralel, int block, double &m, double yMin,
                double yMax);

  inline void setMassFluxes(const Grid &myGrid) {
    PROFILE_FUNCTION();

    massFluxE.getGridInfoPassed(myGrid, var.U.viscX[0], var.U.viscY[0]);
    massFluxN.getGridInfoPassed(myGrid, var.V.viscX[0], var.V.viscY[0]);

    massFluxE.computeEastMassFluxes(var.U);
    massFluxN.computeNorthMassFluxes(var.V);
  }

  inline void setInletBoundaryConditionLeftToRight() {
    PROFILE_FUNCTION();
    var.T.inletBoundaryCondition(west, 0.0);
    var.F.inletBoundaryCondition(west, 1.0);
    var.Z.inletBoundaryCondition(west, 0.0);
  }
  inline void setInletBoundaryConditionRightToLeft() {
    PROFILE_FUNCTION();
    var.T.inletBoundaryCondition(east, 0.0);
    var.F.inletBoundaryCondition(east, 1.0);
    var.Z.inletBoundaryCondition(east, 0.0);
  }

  inline void sendInfoToCommMainProc(Paralel &paralel) {
    PROFILE_FUNCTION();
    paralel.SendInfoToCommMainProc(var.T, sol.T);
    paralel.SendInfoToCommMainProc(var.F, sol.F);
    paralel.SendInfoToCommMainProc(var.Z, sol.Z);
    paralel.SendInfoToCommMainProc(var.U, sol.U);
    paralel.SendInfoToCommMainProc(var.V, sol.V);
  }
  inline void sendInfoToNeighbours(Paralel &paralel) {
    PROFILE_FUNCTION();
    paralel.SendInfoToNeighbours(var.T);
    paralel.SendInfoToNeighbours(var.F);
    paralel.SendInfoToNeighbours(var.Z);
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
    sys.T->SetDirichlet(var.T, side);
    sys.F->SetDirichlet(var.F, side);
    sys.Z->SetDirichlet(var.Z, side);
  }

  inline void assembleEquations(SystemOfEquations &sys) {
    PROFILE_FUNCTION();

    sys.T->assembleEquation();
    sys.F->assembleEquation();
    sys.Z->assembleEquation();

    // sys.T->relax(T);
    // sys.F->relax(F);
    // sys.Z->relax(Z);
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

  Variables sol;

  bool manyIter = false;

 private:
  int ifix, jfix;

  Variables var;

  Field TWall, TNextToWall;

  vector<double> DTinWall;
  Field massFluxE, massFluxN;
  double lowerBoundFactorm = 0.9, upperBoundFactorm = 1.1;
  double DY;

  string initialsol = "Sol-1.f";

 public:
  inline void initialize_top_channel(const double &m, const double &q,
                                     const double &xHS, const double &z0hs) {
    PROFILE_FUNCTION();
    var.U.laminarFlow(m, y_top_min, y_top_max);
    var.V.initializeInternalField(0);
    var.T.InitializeT(q, xHS, channel_xmin, channel_xmax);
    var.F.InitializeF(xHS, channel_xmin, channel_xmax);
    var.Z.InitializeZ(z0hs, r0hs, xHS, (y_top_min + y_top_max) / 2.0);
  }
  inline void initialize_bot_channel(const double &m, const double &q,
                                     const double &xHS, const double &z0hs) {
    PROFILE_FUNCTION();
    var.U.laminarFlow(-m, y_bot_min, y_bot_max);
    var.V.initializeInternalField(0);
    var.T.InitializeT(q, xHS, channel_xmax, channel_xmin);
    var.F.InitializeF(xHS, channel_xmax, channel_xmin);
    var.Z.InitializeZ(z0hs, r0hs, xHS, (y_bot_min + y_bot_max) / 2.0);
  }
};

#endif
