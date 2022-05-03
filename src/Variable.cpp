#include "Variable.hpp"

// Variable::Variable()
// {
// }

Variable::Variable(const int &ni, const int &nj, const int &soln, const int &solm) : NI(ni), NJ(nj), solN(soln), solM(solm),
                                                                                     solT(solN, solM), solF(solN, solM),
                                                                                     solZ(solN, solM), solU(solN, solM),
                                                                                     solV(solN, solM), U(NI, NJ),
                                                                                     V(NI, NJ), T(NI, NJ),
                                                                                     F(NI, NJ), Z(NI, NJ),
                                                                                     TWall(solN + 2, 1), TNextToWall(solN + 2, 1),
                                                                                     massFluxE(NI, NJ), massFluxN(NI, NJ)
{
}

Variable::~Variable() {}

void Variable::passInfoGridToAll(const Grid &solGrid, const Grid &myGrid, double &viscX, double &viscY, double &LeF, double &LeZ)
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

void Variable::initializeLeftToRight(double m, double yMin, double yMax, double q, double xHS, double r0hs, double z0hs, double xMin, double xMax)
{
  PROFILE_FUNCTION();
  U.laminarFlow(m, yMin, yMax);
  V.initializeInternalField(0);
  T.InitializeT(q, xHS, (yMin + yMax) / 2.0, xMin, xMax);
  F.InitializeF(xHS, xMin, xMax);
  Z.InitializeZ(z0hs, r0hs, xHS, (yMin + yMax) / 2.0, xMax);
}

void Variable::initializeRightToLeft(double m, double yMin, double yMax, double q, double xHS, double r0hs, double z0hs, double xMin, double xMax)
{
  PROFILE_FUNCTION();
  U.laminarFlow(-m, yMin, yMax);
  V.initializeInternalField(0);
  T.InitializeT(q, xHS, (yMin + yMax) / 2.0, xMax, xMin);
  F.InitializeF(xHS, xMax, xMin);
  Z.InitializeZ(z0hs, r0hs, xHS, (yMin + yMax) / 2.0, xMax);
}

void Variable::initializeWall()
{
  PROFILE_FUNCTION();
  T.initializeInternalField(1.0);
  T.linearExtrapolateCondition(Field::south);
  T.linearExtrapolateCondition(Field::north);
  T.linearExtrapolateCondition(Field::west);
  T.linearExtrapolateCondition(Field::east);
}

void Variable::setMassFluxes(const Grid &myGrid)
{
  PROFILE_FUNCTION();

  massFluxE.getGridInfoPassed(myGrid, U.viscX[0], U.viscY[0]);
  massFluxN.getGridInfoPassed(myGrid, V.viscX[0], V.viscY[0]);

  massFluxE.computeEastMassFluxes(U);
  massFluxN.computeNorthMassFluxes(V);
}

void Variable::setInletBoundaryConditionLeftToRight()
{
  PROFILE_FUNCTION();
  T.inletBoundaryCondition(Field::west, 0.0);
  F.inletBoundaryCondition(Field::west, 1.0);
  Z.inletBoundaryCondition(Field::west, 0.0);
}

void Variable::setInletBoundaryConditionRightToLeft()
{
  PROFILE_FUNCTION();
  T.inletBoundaryCondition(Field::east, 0.0);
  F.inletBoundaryCondition(Field::east, 1.0);
  Z.inletBoundaryCondition(Field::east, 0.0);
}

void Variable::sendInfoToCommMainProc(Paralel &paralel)
{
  PROFILE_FUNCTION();
  paralel.SendInfoToCommMainProc(T, solT);
  paralel.SendInfoToCommMainProc(F, solF);
  paralel.SendInfoToCommMainProc(Z, solZ);
  paralel.SendInfoToCommMainProc(U, solU);
  paralel.SendInfoToCommMainProc(V, solV);
}

// void Variable::exchangeTemperature(Paralel &paralel, double &exCte)
// {
//   if (paralel.myProc == 0)
//     paralel.ExchangeWallTemperature(solT, TWall, exCte);
//   paralel.ShareWallTemperatureInfo(TWall, T);
// }

void Variable::sendInfoToNeighbours(Paralel &paralel)
{
  PROFILE_FUNCTION();
  paralel.SendInfoToNeighbours(T);
  paralel.SendInfoToNeighbours(F);
  paralel.SendInfoToNeighbours(Z);
}

void Variable::setWallEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double &alphaWall, double &DT)
{
  PROFILE_FUNCTION();

  Teqn = new Equation(alphaWall * fvm::diffusiveTerm(T));
  Feqn = new Equation(fvm::diffusiveTerm(F));
  Zeqn = new Equation(fvm::diffusiveTerm(Z));
  Teqn->DT = DT;
}

void Variable::setChannelEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double &m, double &q, double &beta, double &gamma, double &DT)
{
  PROFILE_FUNCTION();
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

void Variable::setWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side)
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

void Variable::setDirichlet(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side)
{
  PROFILE_FUNCTION();
  Teqn->SetDirichlet(T, side, massFluxE, massFluxN);
  Feqn->SetDirichlet(F, side, massFluxE, massFluxN);
  Zeqn->SetDirichlet(Z, side, massFluxE, massFluxN);
}

void Variable::setExchangeWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, Field::Direction side, int iStr, int iEnd, int mainexI1, int mainexI2, int exI1, int exI2)
{
  PROFILE_FUNCTION();

  Teqn->SetWallShearTX(T, iStr, iEnd, mainexI1, mainexI2,
                       exI1, exI2, massFluxE, massFluxN, side);
  Feqn->SetWallShearX(F, side);
  Zeqn->SetWallShearX(Z, side);
}

void Variable::assembleEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn)
{
  PROFILE_FUNCTION();

  Teqn->assembleEquation();
  Feqn->assembleEquation();
  Zeqn->assembleEquation();

  Teqn->relax(T);
  Feqn->relax(F);
  Zeqn->relax(Z);
}

double Variable::solveEquations(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn, double alpha, int &niter, int &itersol, int changeIter)
{
  PROFILE_FUNCTION();
  double error = 10;
  double newerror = 10;

  FiniteMatrix::finiteMat ST(Teqn->sourceRelaxed);
  error = Teqn->solve(T, ST, alpha, niter, itersol, changeIter);

  FiniteMatrix::finiteMat SF(Feqn->sourceRelaxed);
  newerror = Feqn->solve(F, SF, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  FiniteMatrix::finiteMat SZ(Zeqn->sourceRelaxed);
  newerror = Zeqn->solve(Z, SZ, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  forAllInterior(NI, NJ)
  {
    F.value[i + j * NI] = std::max(F.value[i + j * NI], 0.0);
  }
  return error;
}

double Variable::calculateNewM(double alpha, double m, double q)
{
  PROFILE_FUNCTION();

  int index = ifix + jfix * NI;
  T.value[index] = 0.7;

  double diffusiveTermX = T.viscX[index] * T.Se[index] *
                              (T.value[index + 1] - T.value[index]) / T.DXPtoE[index] -
                          T.viscX[index - 1] * T.Se[index - 1] *
                              (T.value[index] - T.value[index - 1]) / T.DXPtoE[index - 1];

  double diffusiveTermY = T.viscY[index] * T.Sn[index] *
                              (T.value[index + NI] - T.value[index]) / T.DYPtoN[index] -
                          T.viscY[index - NI] * T.Sn[index - NI] *
                              (T.value[index] - T.value[index - NI]) / T.DYPtoN[index - NI];

  double convectiveTermX = (T.value[index + 1] * massFluxE.value[index] * T.FXE[index] -
                            T.value[index - 1] * massFluxE.value[index] * T.FXP[index]);

  double convectiveTermY = (T.value[index + NI] * massFluxN.value[index] * T.FYN[index] -
                            T.value[index - NI] * massFluxN.value[index] * T.FYP[index]);

  double reactionTerm = q * Z.value[index] * T.Se[index] * T.Sn[index];

  double newM = alpha * (reactionTerm + diffusiveTermX + diffusiveTermY) / (convectiveTermX + convectiveTermY) + (1 - alpha) * m;

  newM = std::min(m * upperBoundFactorm, newM);
  newM = std::max(m * lowerBoundFactorm, newM);
  return newM;
}

void Variable::setFixIndex(double xfix, double yfix)
{
  PROFILE_FUNCTION();

  ifix = 0;
  jfix = 0;

  double error = 1000;
  double newerror = abs(T.XC[ifix + jfix * NI] - xfix);

  while (newerror < error)
  {
    error = abs(T.XC[ifix + jfix * NI] - xfix);
    ifix++;
    newerror = abs(T.XC[ifix + jfix * NI] - xfix);
  }
  ifix--;

  error = 1000;
  newerror = abs(T.YC[ifix + jfix * NI] - yfix);

  while (newerror < error)
  {
    error = abs(T.YC[ifix + jfix * NI] - yfix);
    jfix++;
    newerror = abs(T.YC[ifix + jfix * NI] - yfix);
  }
  jfix--;
}

void Variable::updateBoundFactors()
{
  PROFILE_FUNCTION();

  lowerBoundFactorm *= 0.9;
  upperBoundFactorm *= 1.1;
}

void Variable::readFile(Paralel &paralel, int block, double &m, double yMin, double yMax)
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

bool Variable::isTOutEqualToQ(double q, Paralel &paralel)
{
  PROFILE_FUNCTION();

  bool result = false;
  double expectedTOut = q * solM / 2.0;
  if (paralel.myProc == 0 && paralel.isRightToLeft())
  {
    double TOut = 0;
    for (int j = 0; j < solM; j++)
      TOut += solT.value[1 + j * solN];

    if (abs(TOut - expectedTOut) < 10e-3)
      result = true;
  }
  else if (paralel.myProc == 0 && paralel.isLeftToRight())
  {
    double TOut = 0;
    for (int j = 0; j < solM; j++)
      TOut += solT.value[solN - 2 + j * solN];

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

void Variable::exchangeTemperature(Paralel &paralel, double &exCte, int &solExI1, int &solExI2)
{
  PROFILE_FUNCTION();

  paralel.GatherWallTemperature(TWall, TNextToWall, T, exCte);

  if (paralel.myProc == 0)
    paralel.ExchangeWallTemperature(TWall, TNextToWall, exCte, solExI1, solExI2);
  MPI_Barrier(MPI_COMM_WORLD);

  paralel.ShareWallTemperatureInfo(TWall, T);
}

void Variable::writeTInWall(Paralel &paralel, const Grid &mainGrid, const Grid &myGrid, int iter)
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
