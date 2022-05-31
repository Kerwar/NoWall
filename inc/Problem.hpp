#ifndef PROBLEM_H
#define PROBLEM_H

#include <memory>

#include "Equation.hpp"
#include "FileWriter.hpp"
#include "Grid.hpp"
#include "Paralel.hpp"
#include "Variable.hpp"

template <int N, int M, int NPROCS>
class Problem {
  static constexpr int NI = N / (NPROCS / 2) + 2;
  static constexpr int NJ = M / 2 + 2;

 public:
  Problem(const double &prevm);
  virtual ~Problem();

  void setUpProblem();
  void initializeVariables();
  void writeSolution(string &filename, int i);
  void freeComm();
  double mainIter(int i);
  void inline setExchangeConstant() {
    exCte = bK * a * a / 2.0;
    exCte /= (1.0 + mainGrid.DY * exCte);
  };
  void writefilename(string &filename);
  void retrieveNandM(int &nxOut, int &nyOut, double &mOut);
  inline bool fixPointInThisProc();
  inline bool isSolutionNotGoodEnough() {
    return !variables.isTOutEqualToQ(q, paralel);
  };

  int maxit, iShow;
  double m;
  Paralel<N, M, NPROCS> paralel;

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
  Variable<N, M, NI, NJ, NPROCS> variables;

  int myProc, nProcs;
  double myXMin, myXMax;
  double myYMin, myYMax;

  int solN, solM;

  Equation *Teqn;
  Equation *Feqn;
  Equation *Zeqn;

  string sufix;

 public:
  inline bool isErrorSmallEnough(double &error) {
    return error < tolerance ? true : false;
  };
};

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::beta;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::gamma;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::xMin;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::xMax;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::xExMin;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::xExMax;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::yWall;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::yChannel;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::yMin;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::yMax;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::xfix;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::yfix;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::q;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::xHS_U;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::xHS_D;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::r0hs;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::z0hs;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::alphaWall;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::DT;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::LeF;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::LeZ;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::a;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::a2;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::bK;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::exCte;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::viscX;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::viscY;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::tolerance;

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::alpha;

template <int N, int M, int NPROCS>
bool Problem<N, M, NPROCS>::readFromFile;

template <int N, int M, int NPROCS>
Problem<N, M, NPROCS>::Problem(const double &prevm)
    : maxit(10E2), iShow(10E1), m(prevm) {
  PROFILE_FUNCTION();
  readFromFile = false;
  tolerance = 10e-12;
  alpha = 0.95;

  beta = 10;
  gamma = 0.7;
  xMin = 0;
  xMax = 120;
  xExMin = 40;
  xExMax = 80;
  yWall = 0.0;
  yChannel = 0.5;
  yMin = 0;
  yMax = 1;
  xfix = 43.95;
  yfix = 0.75;
  q = 1.2;
  xHS_U = 44;
  xHS_D = 76;
  r0hs = 1;
  z0hs = 0.5;
  alphaWall = 1;
  DT = 10e-3;
  itersol = 4;
  LeF = 1.0;
  LeZ = 0.3;
  a = 0.1;
  a2 = 1 / (a * a);
  bK = 0.15;
  viscX = 1;
  viscY = a2;  // if (nx > 20000)
  // {
  // tolerance = 10e-6;
  // maxit = 1000000;}
}

template <int N, int M, int NPROCS>
Problem<N, M, NPROCS>::~Problem() {}

template <int N, int M, int NPROCS>
void Problem<N, M, NPROCS>::setUpProblem() {
  PROFILE_FUNCTION();
  mainGrid = Grid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);

  myProc = paralel.myProc;
  nProcs = paralel.nProcs;

  myXMin = mainGrid.X[paralel.iStr];
  myXMax = mainGrid.X[paralel.iEnd];
  myYMin = 0.0;
  myYMax = yChannel;

  if (paralel.isLeftToRight()) {
    myYMin += yChannel + yWall;
    myYMax += yChannel + yWall;
  }

  myGrid = Grid(NI - 2, NJ - 2, myXMin, myXMax, myYMin, myYMax);
  myGrid.SetIEx(xExMin, xExMax);

  solN = N;
  solM = M;
}

template <int N, int M, int NPROCS>
void Problem<N, M, NPROCS>::initializeVariables() {
  PROFILE_FUNCTION();
  variables.passInfoGridToAll(mainGrid, myGrid, viscX, viscY, LeF, LeZ);
  setExchangeConstant();
  int fixPointProcBuffer = 0;

  if (fixPointInThisProc()) {
    fixPointProcBuffer += paralel.worldMyProc;
    variables.setFixIndex(xfix, yfix);
  }

  MPI_Allreduce(&fixPointProcBuffer, &paralel.fixPointProc, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

  if (!readFromFile) {
    if (paralel.isLeftToRight())
      variables.initializeLeftToRight(m, myYMin, myYMax + yChannel, q, xHS_U,
                                      r0hs, z0hs, xMin, xMax);
    else if (paralel.isRightToLeft())
      variables.initializeRightToLeft(m, myYMin - yChannel, myYMax, q, xHS_D,
                                      r0hs, z0hs, xMin, xMax);

    readFromFile = true;
  } else {
    if (paralel.isLeftToRight())
      variables.readFile(paralel, 2, m, myYMin, myYMax + yChannel);

    MPI_Barrier(MPI_COMM_WORLD);
    if (paralel.isRightToLeft())
      variables.readFile(paralel, 1, m, myYMin - yChannel, myYMax);

    variables.sendInfoToNeighbours(paralel);
  }
  //  If density non constant this interpolate the Velocities and the mass
  //  fluxes have to be calculated with them
  // Field::vectorField UE(fieldOper.interpolatedFieldEast(U, myGrid));
  // Field::vectorField UN(fieldOper.interpolatedFieldNorth(V, myGrid));

  // fieldOper.getGridInfoPassed(UE, myGrid, viscX, viscY);
  // fieldOper.getGridInfoPassed(UN, myGrid, viscX, viscY);

  variables.setMassFluxes(myGrid);
  if (paralel.isLeftToRight() && paralel.isProcNull(paralel.myLeft)) {
    variables.setInletBoundaryConditionLeftToRight();
  } else if (paralel.isRightToLeft() && paralel.isProcNull(paralel.myRight)) {
    variables.setInletBoundaryConditionRightToLeft();
  }

  variables.sendInfoToCommMainProc(paralel);

  variables.exchangeTemperature(paralel, mainGrid);
}

template <int N, int M, int NPROCS>
void Problem<N, M, NPROCS>::writeSolution(string &prefix, int i) {
  PROFILE_FUNCTION();

  FileWriter fileWriter;

  if (i == 0) {
    sufix = "";
    writefilename(sufix);
  } else if (i == -1) {
    sufix = "";
    variables.writeTInWall(paralel, mainGrid, myGrid, 0);
  }

  variables.sendInfoToCommMainProc(paralel);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myProc == 0 && paralel.isRightToLeft())
    fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.solU,
                          variables.solV, variables.solT, variables.solF,
                          variables.solZ, paralel.locIStr, paralel.locIEnd,
                          paralel.locJStr, paralel.loc);
  // fileWriter.WriteBin(prefix, sufix, i, mainGrid, variables.solU,
  // variables.solV, variables.solT, variables.solF, variables.solZ,
  // paralel.locIStr, paralel.locIEnd, paralel.locJStr, paralel.locJEnd,
  // paralel.loc);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myProc == 0 && paralel.isLeftToRight())
    fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.solU,
                          variables.solV, variables.solT, variables.solF,
                          variables.solZ, paralel.locIStr, paralel.locIEnd,
                          paralel.locJStr, paralel.loc);
  // fileWriter.WriteBin(prefix, sufix, i, mainGrid, variables.solU,
  // variables.solV, variables.solT, variables.solF, variables.solZ,
  // paralel.locIStr, paralel.locIEnd, paralel.locJStr, paralel.locJEnd,
  // paralel.loc);
  MPI_Barrier(MPI_COMM_WORLD);
}

template <int N, int M, int NPROCS>
double Problem<N, M, NPROCS>::mainIter(int i) {
  PROFILE_FUNCTION();

  double error = 1;

  variables.setChannelEquations(Teqn, Feqn, Zeqn, paralel, m, q, beta, gamma,
                                DT, exCte, i);
  Direction side = west;

  if (paralel.isFirstProcOfRow() && paralel.isRightToLeft())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);
  else
    variables.setDirichlet(Teqn, Feqn, Zeqn, side);

  side = east;

  if (paralel.isLastProcOfRow() && paralel.isLeftToRight())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);
  else
    variables.setDirichlet(Teqn, Feqn, Zeqn, side);

  side = south;

  if (paralel.isRightToLeft())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);
  else if (paralel.isLeftToRight())
    variables.setExchangeWallShear(Teqn, Feqn, Zeqn, side, paralel.iStr,
                                   paralel.iEnd, mainGrid.exI1, mainGrid.exI2,
                                   myGrid.exI1, myGrid.exI2);

  side = north;

  if (paralel.isRightToLeft())
    variables.setExchangeWallShear(Teqn, Feqn, Zeqn, side, paralel.iStr,
                                   paralel.iEnd, mainGrid.exI1, mainGrid.exI2,
                                   myGrid.exI1, myGrid.exI2);
  else if (paralel.isLeftToRight())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);

  variables.assembleEquations(Teqn, Feqn, Zeqn);

  Teqn->EqnName = "T-Eqn";
  Feqn->EqnName = "F-Eqn";
  Zeqn->EqnName = "Z-Eqn";

  error =
      variables.solveEquations(Teqn, Feqn, Zeqn, alpha, i, itersol, 0 * iShow);

  double error_result = 0;

  if (fixPointInThisProc())
    m = variables.calculateNewM(0.75 * alpha, m, q);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&m, 1, MPI_DOUBLE, paralel.fixPointProc, MPI_COMM_WORLD);

  MPI_Allreduce(&error, &error_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  variables.sendInfoToNeighbours(paralel);

  variables.exchangeTemperature(paralel, mainGrid);

  MPI_Barrier(MPI_COMM_WORLD);
  return error_result;
}

template <int N, int M, int NPROCS>
void Problem<N, M, NPROCS>::freeComm() {
  paralel.freeComm();
}

template <int N, int M, int NPROCS>
inline bool Problem<N, M, NPROCS>::fixPointInThisProc() {
  if (myXMin <= xfix && xfix <= myXMax)
    if (myYMin <= yfix && yfix <= myYMax) return true;

  return false;
}

template <int N, int M, int NPROCS>
void Problem<N, M, NPROCS>::writefilename(string &filename) {
  PROFILE_FUNCTION();
  std::ostringstream temp;

  filename.append("_NxM-");
  temp.str("");
  temp.clear();
  temp << N;
  filename.append(temp.str());

  filename.append("x");
  temp.str("");
  temp.clear();
  temp << M;
  filename.append(temp.str());

  filename.append("_Lxa-");
  temp.str("");
  temp.clear();
  temp << xMax - xMin;
  filename.append(temp.str());

  filename.append("x");
  temp.str("");
  temp.clear();
  temp << a;
  filename.append(temp.str());

  filename.append("_Ex-");
  temp.str("");
  temp.clear();
  temp << xExMax - xExMin;
  filename.append(temp.str());

  filename.append("_q-");
  temp.str("");
  temp.clear();
  temp << q;
  filename.append(temp.str());

  filename.append("_m-");
  temp.str("");
  temp.clear();
  temp << m;
  filename.append(temp.str());

  filename.append("_beta-");
  temp.str("");
  temp.clear();
  temp << beta;
  filename.append(temp.str());

  filename.append("_LeFxZ-");
  temp.str("");
  temp.clear();
  temp << LeF;
  filename.append(temp.str());

  filename.append("x");
  temp.str("");
  temp.clear();
  temp << LeZ;
  filename.append(temp.str());

  filename.append("_");
}

template <int N, int M, int NPROCS>
void Problem<N, M, NPROCS>::retrieveNandM(int &nxOut, int &nyOut,
                                          double &mOut) {
  PROFILE_FUNCTION();
  // double xcellSize = (xMax - xMin) / N;
  // double ycellSize = (yMax - yMin) / M;

  nxOut = N * 5;
  nyOut = M * 2;
  // if (xcellSize >= ycellSize * a)
  //  {
  //  nxOut = N * 2;
  // nyOut = M;
  //  }
  mOut = m;
}

#endif