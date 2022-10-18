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

  void set_mass_fluxes(const Grid &myGrid);

  void set_inlet_up_channel();
  void set_inlet_bot_channel();

  void sendInfoToCommMainProc(Paralel &paralel);
  void sendInfoToNeighbours(Paralel &paralel);

  void exchangeTemperature(const Paralel &paralel, const Grid &mainGrid);

  void setChannelEquations(SystemOfEquations &sys, const Paralel &paralel,
                           const double &m, const double &q, double &DT,
                           const int &iter);
  void setWallShear(SystemOfEquations &sys, Direction side);
  void setExchangeWallShear(SystemOfEquations &sys, Direction side);
  void setDirichlet(SystemOfEquations &sys, Direction side);

  void assembleEquations(SystemOfEquations &sys);
  double solveEquations(SystemOfEquations &sys, double alpha, int &niter,
                        int &itersol, int changeIter);

  double calculateNewM(const double &alpha, const double &m, const double &q);
  double calculateNewQ(const double &alpha, const double &m, const double &q);

  void setFixIndex(double xfix, double yfix);

  bool isTOutEqualToQ(double q, const Paralel &paralel);

  Variables sol;

  bool manyIter = false;

  void initialize_top_channel(const double &m, const double &q,
                              const double &xHS);
  void initialize_bot_channel(const double &m, const double &q,
                              const double &xHS);

  void setExchangeIndexes(const Grid &mainGrid, const Grid &myGrid, const Paralel &paralel);

 private:
  int ifix, jfix;
  int end_ex_zone_1, begin_ex_zone_2;

  Variables var;

  Field TWall, TNextToWall;

  vector<double> DTinWall;
  Field massFluxE, massFluxN;
  double lowerBoundFactorm = 0.9, upperBoundFactorm = 1.1;
  double DY;

  string initialsol = "Sol-1.f";

};

#endif
