#ifndef VARIABLE_H
#define VARIABLE_H

#include "Field.hpp"
#include "Grid.hpp"
#include "Paralel.hpp"
#include "Equation.hpp"
#include "FiniteVolumeOperations.hpp"
#include "FileReader.hpp"

template <int N, int M, int NI, int NJ, int NPROCS>
class Variable
{
public:
  // Variable();
  Variable();
  virtual ~Variable();

  void passInfoGridToAll(const Grid &solGrid, const Grid &myGrid, double &viscX, double &viscY, double &LeF, double &LeZ);
  void initializeWall();
  void readFile(Paralel<N, M, NPROCS> &paralel, int block, double &m, double yMin, double yMax);
  void setMassFluxes(const Grid &myGrid);
  void setInletBoundaryConditionLeftToRight();
  void setInletBoundaryConditionRightToLeft();
  void sendInfoToCommMainProc(Paralel<N, M, NPROCS> &paralel);
  void sendInfoToNeighbours(Paralel<N, M, NPROCS> &paralel);
  void exchangeTemperature(Paralel<N, M, NPROCS> &paralel, double &exCte, int &solExI1, int &solExI2);
  void setChannelEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, const double &m, const double &q,
                           const double &beta, double &gamma, double &DT,
                           const int &iter);
  void setWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side);
  void setExchangeWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side, int iStr, int iEnd, int mainexI1, int mainexI2, int exI1, int exI2);
  void setDirichlet(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side);
  void assembleEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn);
  double solveEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double alpha, int &niter, int &itersol, int changeIter);
  double calculateNewM(const double &alpha, const double &m, const double &q);
  void setFixIndex(double xfix, double yfix);
  void updateBoundFactors();
  bool isTOutEqualToQ(double q, Paralel<N, M, NPROCS> &paralel);
  void writeTInWall(Paralel<N, M, NPROCS> &paralel, const Grid &mainGrid, const Grid &myGrid, int iter);

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

public:
  inline void initializeLeftToRight(const double &m, const double &yMin, const double &yMax,
                                    const double &q, const double &xHS, const double &r0hs, const double &z0hs, const double &xMin, const double &xMax)
  {
    PROFILE_FUNCTION();
    U.laminarFlow(m, yMin, yMax);
    V.initializeInternalField(0);
    T.InitializeT(q, xHS, xMin, xMax);
    F.InitializeF(xHS, xMin, xMax);
    Z.InitializeZ(z0hs, r0hs, xHS, (yMin + yMax) / 2.0);
  }
  inline void initializeRightToLeft(const double &m, const double &yMin, const double &yMax,
                                    const double &q, const double &xHS, const double &r0hs, const double &z0hs, const double &xMin, const double &xMax)
  {
    PROFILE_FUNCTION();
    U.laminarFlow(-m, yMin, yMax);
    V.initializeInternalField(0);
    T.InitializeT(q, xHS, xMax, xMin);
    F.InitializeF(xHS, xMax, xMin);
    Z.InitializeZ(z0hs, r0hs, xHS, (yMin + yMax) / 2.0);
  }
};

template <int N, int M, int NI, int NJ, int NPROCS>
Variable<N, M, NI, NJ, NPROCS>:: Variable() : solT(N, M), solF(N, M),
                                            solZ(N, M), solU(N, M),
                                            solV(N, M), U(NI, NJ),
                                            V(NI, NJ), T(NI, NJ),
                                            F(NI, NJ), Z(NI, NJ),
                                            TWall(N + 2, 1), TNextToWall(N + 2, 1),
                                            massFluxE(NI, NJ), massFluxN(NI, NJ)
{
}

template <int N, int M, int NI, int NJ, int NPROCS>
Variable<N, M, NI, NJ, NPROCS>:: ~Variable() {}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: passInfoGridToAll(const Grid &solGrid, const Grid &myGrid, double &viscX, double &viscY, double &LeF, double &LeZ)
{
  PROFILE_FUNCTION();

  double viscTx = viscX;
  double viscTy = viscY;

  double viscFx = viscX / LeF;
  double viscFy = viscY / LeF;

  double viscZx = viscX / LeZ;
  double viscZy = viscY / LeZ;

  solU.getGridInfoPassed(solGrid, viscTx, viscTy);
  solV.getGridInfoPassed(solGrid, viscTx, viscTy);
  solT.getGridInfoPassed(solGrid, viscTx, viscTy);
  solF.getGridInfoPassed(solGrid, viscFx, viscFy);
  solZ.getGridInfoPassed(solGrid, viscZy, viscZy);

  U.getGridInfoPassed(myGrid, viscX, viscY);
  V.getGridInfoPassed(myGrid, viscX, viscY);

  T.getGridInfoPassed(myGrid, viscTx, viscTy);
  F.getGridInfoPassed(myGrid, viscFx, viscFy);
  Z.getGridInfoPassed(myGrid, viscZx, viscZy);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: initializeWall()
{
  PROFILE_FUNCTION();
  T.initializeInternalField(1.0);
  T.linearExtrapolateCondition(Field::south);
  T.linearExtrapolateCondition(Field::north);
  T.linearExtrapolateCondition(Field::west);
  T.linearExtrapolateCondition(Field::east);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setMassFluxes(const Grid &myGrid)
{
  PROFILE_FUNCTION();

  massFluxE.getGridInfoPassed(myGrid, U.viscX[0], U.viscY[0]);
  massFluxN.getGridInfoPassed(myGrid, V.viscX[0], V.viscY[0]);

  massFluxE.computeEastMassFluxes(U);
  massFluxN.computeNorthMassFluxes(V);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setInletBoundaryConditionLeftToRight()
{
  PROFILE_FUNCTION();
  T.inletBoundaryCondition(Field::west, 0.0);
  F.inletBoundaryCondition(Field::west, 1.0);
  Z.inletBoundaryCondition(Field::west, 0.0);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setInletBoundaryConditionRightToLeft()
{
  PROFILE_FUNCTION();
  T.inletBoundaryCondition(Field::east, 0.0);
  F.inletBoundaryCondition(Field::east, 1.0);
  Z.inletBoundaryCondition(Field::east, 0.0);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>::sendInfoToCommMainProc(Paralel<N, M, NPROCS> &paralel)
{
  PROFILE_FUNCTION();
  paralel.SendInfoToCommMainProc(T, solT);
  paralel.SendInfoToCommMainProc(F, solF);
  paralel.SendInfoToCommMainProc(Z, solZ);
  paralel.SendInfoToCommMainProc(U, solU);
  paralel.SendInfoToCommMainProc(V, solV);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: sendInfoToNeighbours(Paralel<N, M, NPROCS> &paralel)
{
  PROFILE_FUNCTION();
  paralel.SendInfoToNeighbours(T);
  paralel.SendInfoToNeighbours(F);
  paralel.SendInfoToNeighbours(Z);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setChannelEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, const double &m, const double &q,
                                                        const double &beta, double &gamma, double &DT,
                                                        const int &iter)
{
  PROFILE_FUNCTION();
  if (iter == 1)
  {
    Teqn = new Equation(fvm::diffusiveTerm(T) +
                        fvm::convectiveTerm(T, massFluxE, massFluxN, m) +
                        fvm::heatProduction(T, Z, q));

    Feqn = new Equation(fvm::diffusiveTerm(F) +
                        fvm::convectiveTerm(F, massFluxE, massFluxN, m) -
                        fvm::intermidiateReaction(F, Z, T, beta, gamma));

    Zeqn = new Equation(fvm::diffusiveTerm(Z) +
                        fvm::convectiveTerm(Z, massFluxE, massFluxN, m) +
                        fvm::intermidiateReaction(F, Z, T, beta, gamma) - fvm::zComsumptium(Z));

    Teqn->DT = DT;
    Feqn->DT = DT;
    Zeqn->DT = DT;
  }
  else
  {
    Teqn->updateEquation(fvm::diffusiveTerm(T) +
                         fvm::convectiveTerm(T, massFluxE, massFluxN, m) +
                         fvm::heatProduction(T, Z, q));

    Feqn->updateEquation(fvm::diffusiveTerm(F) +
                         fvm::convectiveTerm(F, massFluxE, massFluxN, m) -
                         fvm::intermidiateReaction(F, Z, T, beta, gamma));

    Zeqn->updateEquation(fvm::diffusiveTerm(Z) +
                         fvm::convectiveTerm(Z, massFluxE, massFluxN, m) +
                         fvm::intermidiateReaction(F, Z, T, beta, gamma) - fvm::zComsumptium(Z));
  }
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side)
{
  PROFILE_FUNCTION();
  if (side == Field::west || side == Field::east)
  {
    Teqn->SetWallShearY(T, side);
    Feqn->SetWallShearY(F, side);
    Zeqn->SetWallShearY(Z, side);
  }
  else
  {
    Teqn->SetWallShearX(T, side);
    Feqn->SetWallShearX(F, side);
    Zeqn->SetWallShearX(Z, side);
  }
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setDirichlet(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side)
{
  PROFILE_FUNCTION();
  Teqn->SetDirichlet(T, side);
  Feqn->SetDirichlet(F, side);
  Zeqn->SetDirichlet(Z, side);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setExchangeWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side, int iStr, int iEnd, int mainexI1, int mainexI2, int exI1, int exI2)
{
  PROFILE_FUNCTION();

  Teqn->SetWallShearTX(T, iStr, iEnd, mainexI1, mainexI2,
                       exI1, exI2, side);
  Feqn->SetWallShearX(F, side);
  Zeqn->SetWallShearX(Z, side);
  
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: assembleEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn)
{
  PROFILE_FUNCTION();

  Teqn->assembleEquation();
  Feqn->assembleEquation();
  Zeqn->assembleEquation();

  Teqn->relax(T);
  Feqn->relax(F);
  Zeqn->relax(Z);
}

template <int N, int M, int NI, int NJ, int NPROCS>
double Variable<N, M, NI, NJ, NPROCS>:: solveEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double alpha, int &niter, int &itersol, int changeIter)
{
  PROFILE_FUNCTION();
  double error = 10;
  double newerror = 10;

  error = Teqn->solve(T, alpha, niter, itersol, changeIter);

  newerror = Feqn->solve(F, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  newerror = Zeqn->solve(Z, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  forAllInterior(NI, NJ)
  {
    F.value[i + j * NI] = std::max(F.value[i + j * NI], 0.0);
  }
  return error;
}

template <int N, int M, int NI, int NJ, int NPROCS>
double Variable<N, M, NI, NJ, NPROCS>:: calculateNewM(const double &alpha, const double &m, const double &q)
{
  PROFILE_FUNCTION();
  int index = ifix + jfix * NI;
  T.value[index] = 0.7;

  double diffusiveTermX = T.viscX[index] * T.Se[jfix] *
                              (T.value[index + 1] - T.value[index]) / T.DXPtoE[ifix] -
                          T.viscX[index - 1] * T.Se[jfix - 1] *
                              (T.value[index] - T.value[index - 1]) / T.DXPtoE[ifix - 1];

  double diffusiveTermY = T.viscY[index] * T.Sn[ifix] *
                              (T.value[index + NI] - T.value[index]) / T.DYPtoN[jfix] -
                          T.viscY[index - NI] * T.Sn[ifix - 1] *
                              (T.value[index] - T.value[index - NI]) / T.DYPtoN[jfix];

  double convectiveTermX = (T.value[index + 1] * massFluxE.value[index] * T.FXE[ifix] -
                            T.value[index - 1] * massFluxE.value[index] * T.FXP[ifix]);

  double convectiveTermY = (T.value[index + NI] * massFluxN.value[index] * T.FYN[jfix] -
                            T.value[index - NI] * massFluxN.value[index] * T.FYP[jfix]);

  double reactionTerm = q * Z.value[index] * T.Se[jfix] * T.Sn[ifix];

  double newM = alpha * (reactionTerm + diffusiveTermX + diffusiveTermY) / (convectiveTermX + convectiveTermY) + (1 - alpha) * m;

  newM = std::min(m * upperBoundFactorm, newM);
  newM = std::max(m * lowerBoundFactorm, newM);
  return newM;
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: setFixIndex(double xfix, double yfix)
{
  PROFILE_FUNCTION();

  ifix = 0;
  jfix = 0;

  double error = 1000;
  double newerror = abs(T.XC[ifix] - xfix);

  while (newerror < error)
  {
    error = abs(T.XC[ifix] - xfix);
    ifix++;
    newerror = abs(T.XC[ifix] - xfix);
  }
  ifix--;

  error = 1000;
  newerror = abs(T.YC[jfix] - yfix);

  while (newerror < error)
  {
    error = abs(T.YC[jfix] - yfix);
    jfix++;
    newerror = abs(T.YC[jfix] - yfix);
  }
  jfix--;
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: updateBoundFactors()
{
  PROFILE_FUNCTION();

  lowerBoundFactorm *= 0.9;
  upperBoundFactorm *= 1.1;
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: readFile(Paralel<N, M, NPROCS> &paralel, int block, double &m, double yMin, double yMax)
{
  PROFILE_FUNCTION();

  FileReader fileread;

  double laminarm = m;
  if (paralel.isRightToLeft())
    laminarm = -m;

  U.laminarFlow(laminarm, yMin, yMax);
  V.initializeInternalField(0);

  if (paralel.myProc == 0)
  {
    fileread.readField(initialsol, block, 1, solT, paralel.locIStr, paralel.locIEnd);
    fileread.readField(initialsol, block, 2, solF, paralel.locIStr, paralel.locIEnd);
    fileread.readField(initialsol, block, 3, solZ, paralel.locIStr, paralel.locIEnd);
  }

  MPI_Barrier(paralel.myComm);

  paralel.distributeToProcs(solT, T);
  paralel.distributeToProcs(solF, F);
  paralel.distributeToProcs(solZ, Z);
}

template <int N, int M, int NI, int NJ, int NPROCS>
bool Variable<N, M, NI, NJ, NPROCS>:: isTOutEqualToQ(double q, Paralel<N, M, NPROCS> &paralel)
{
  PROFILE_FUNCTION();

  bool result = false;
  double expectedTOut = q * M / 2.0;
  if (paralel.myProc == 0 && paralel.isRightToLeft())
  {
    double TOut = 0;
    for (int j = 0; j < M; j++)
      TOut += solT.value[1 + j * N];

    if (abs(TOut - expectedTOut) < 10e-3)
      result = true;
  }
  else if (paralel.myProc == 0 && paralel.isLeftToRight())
  {
    double TOut = 0;
    for (int j = 0; j < M; j++)
      TOut += solT.value[N - 2 + j * N];

    if (abs(TOut - expectedTOut) < 10e-3)
      result = true;
  }

  bool buffer = result;

  MPI_Bcast(&buffer, 1, MPI_CXX_BOOL, paralel.UCMainProc, MPI_COMM_WORLD);

  result = buffer || result;

  MPI_Bcast(&buffer, 1, MPI_CXX_BOOL, paralel.DCMainProc, MPI_COMM_WORLD);

  result = buffer || result;

  return result;
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: exchangeTemperature(Paralel<N, M, NPROCS> &paralel, double &exCte, int &solExI1, int &solExI2)
{
  PROFILE_FUNCTION();

  paralel.GatherWallTemperature(TWall, TNextToWall, T);

  if (paralel.myProc == 0)
    paralel.ExchangeWallTemperature(TWall, TNextToWall, exCte, solExI1, solExI2);
  MPI_Barrier(MPI_COMM_WORLD);

  paralel.ShareWallTemperatureInfo(TWall, T);
}

template <int N, int M, int NI, int NJ, int NPROCS>
void Variable<N, M, NI, NJ, NPROCS>:: writeTInWall(Paralel<N, M, NPROCS> &paralel, const Grid &mainGrid, const Grid &myGrid, int iter)
{
  PROFILE_FUNCTION();

  string file = "TW";
  std::ostringstream filename;
  filename << file << "--" << paralel.worldMyProc << "--" << iter << ".txt";
  char cstr[filename.str().size() + 2];
  strcpy(cstr, filename.str().c_str());

  std::ofstream outfile;
  outfile.open(cstr, std::ios::out);
  outfile << myGrid.exI1 << " " << myGrid.exI2 << std::endl;
  for (int i = paralel.iStr; i < paralel.iEnd; i++)
  {
    if (paralel.isLeftToRight())
      outfile << i << " " << i - paralel.iStr + 1 << " " << mainGrid.XC[i + 1] << " " << myGrid.XC[i - paralel.iStr + 1] << " " << TWall.value[i + 1] << " " << T.value[i - paralel.iStr + 1] << " " << std::endl;
    else
      outfile << i << " " << i - paralel.iStr + 1 << " " << mainGrid.XC[i + 1] << " " << myGrid.XC[i - paralel.iStr + 1] << " " << TWall.value[i + 1] << " " << T.value[i - paralel.iStr + 1 + (NJ - 1) * NI] << " " << std::endl;
  }

  outfile.close();
}

#endif
