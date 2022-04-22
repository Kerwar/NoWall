#include "Variable.h"

// Variable::Variable()
// {
// }

Variable::Variable(int &n, int &m, int &soln, int &solm) : NI(n), NJ(m), solN(soln), solM(solm),
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

void Variable::passInfoGridToAll(Grid &solGrid, Grid &myGrid, double &viscX, double &viscY, double &LeF, double &LeZ)
{
  PROFILE_FUNCTION();

  double viscTx = viscX;
  double viscTy = viscY;

  double viscFx = viscX / LeF;
  double viscFy = viscY / LeF;

  double viscZx = viscX / LeZ;
  double viscZy = viscY / LeZ;

  fieldOper.getGridInfoPassed(solU, solGrid, viscTx, viscTy);
  fieldOper.getGridInfoPassed(solV, solGrid, viscTx, viscTy);

  fieldOper.getGridInfoPassed(solT, solGrid, viscTx, viscTy);
  fieldOper.getGridInfoPassed(solF, solGrid, viscFx, viscFy);
  fieldOper.getGridInfoPassed(solZ, solGrid, viscZy, viscZy);

  fieldOper.getGridInfoPassed(U, myGrid, viscX, viscY);
  fieldOper.getGridInfoPassed(V, myGrid, viscX, viscY);
  fieldOper.getGridInfoPassed(T, myGrid, viscTx, viscTy);
  fieldOper.getGridInfoPassed(F, myGrid, viscFx, viscFy);
  fieldOper.getGridInfoPassed(Z, myGrid, viscZx, viscZy);
}

void Variable::initializeLeftToRight(double m, double yMin, double yMax, double q, double xHS, double r0hs, double z0hs, double xMin, double xMax)
{
  PROFILE_FUNCTION();
  fieldOper.laminarFlow(U, m, yMin, yMax);
  fieldOper.initializeInternalField(V, 0);
  fieldOper.InitializeT(T, q, xHS, (yMin + yMax) / 2.0, xMin, xMax);
  fieldOper.InitializeF(F, xHS, xMin, xMax);
  fieldOper.InitializeZ(Z, z0hs, r0hs, xHS, (yMin + yMax) / 2.0, xMax);
}

void Variable::initializeRightToLeft(double m, double yMin, double yMax, double q, double xHS, double r0hs, double z0hs, double xMin, double xMax)
{
  PROFILE_FUNCTION();
  fieldOper.laminarFlow(U, -m, yMin, yMax);
  fieldOper.initializeInternalField(V, 0);
  fieldOper.InitializeT(T, q, xHS, (yMin + yMax) / 2.0, xMax, xMin);
  fieldOper.InitializeF(F, xHS, xMax, xMin);
  fieldOper.InitializeZ(Z, z0hs, r0hs, xHS, (yMin + yMax) / 2.0, xMax);
}

void Variable::initializeWall()
{
  PROFILE_FUNCTION();
  fieldOper.initializeInternalField(T, 1.0);
  fieldOper.linearExtrapolateCondition(T, Field::south);
  fieldOper.linearExtrapolateCondition(T, Field::north);
  fieldOper.linearExtrapolateCondition(T, Field::west);
  fieldOper.linearExtrapolateCondition(T, Field::east);
}

void Variable::setMassFluxes(Grid &myGrid)
{
  PROFILE_FUNCTION();

  fieldOper.getGridInfoPassed(massFluxE, myGrid, U.viscX[0], U.viscY[0]);
  fieldOper.getGridInfoPassed(massFluxN, myGrid, V.viscX[0], V.viscY[0]);

  fieldOper.computeEastMassFluxes(massFluxE, U);
  fieldOper.computeNorthMassFluxes(massFluxN, V);
}

void Variable::setInletBoundaryConditionLeftToRight()
{
  PROFILE_FUNCTION();
  fieldOper.inletBoundaryCondition(T, Field::west, 0.0);
  fieldOper.inletBoundaryCondition(F, Field::west, 1.0);
  fieldOper.inletBoundaryCondition(Z, Field::west, 0.0);
}

void Variable::setInletBoundaryConditionRightToLeft()
{
  PROFILE_FUNCTION();
  fieldOper.inletBoundaryCondition(T, Field::east, 0.0);
  fieldOper.inletBoundaryCondition(F, Field::east, 1.0);
  fieldOper.inletBoundaryCondition(Z, Field::east, 0.0);
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
  double error;

  FiniteMatrix::finiteMat ST(Teqn->sourceRelaxed);
  error = Teqn->solve(T, ST, alpha, niter, itersol, changeIter);

  FiniteMatrix::finiteMat SF(Feqn->sourceRelaxed);
  error = std::max(error, Feqn->solve(F, SF, alpha, niter, itersol, changeIter));

  FiniteMatrix::finiteMat SZ(Zeqn->sourceRelaxed);
  error = std::max(error, Zeqn->solve(Z, SZ, alpha, niter, itersol, changeIter));

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
  Field fieldOper;

  double laminarm = m;
  if (paralel.isRightToLeft())
    laminarm = -m;

  fieldOper.laminarFlow(U, laminarm, yMin, yMax);
  fieldOper.initializeInternalField(V, 0);

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

void Variable::exchangeTemperature(Paralel &paralel, double &exCte)
{
  PROFILE_FUNCTION();

  paralel.GatherWallTemperature(TWall, TNextToWall, T, exCte);

  if (paralel.myProc == 0)
    paralel.ExchangeWallTemperature(TWall, TNextToWall, exCte);
  MPI_Barrier(MPI_COMM_WORLD);

  paralel.ShareWallTemperatureInfo(TWall, T);
}

void Variable::writeTInWall(Paralel &paralel, int iter)
{
  PROFILE_FUNCTION();

  string file = "TW32";
  std::ostringstream filename;
  filename << file << " " << paralel.worldMyProc << " " << iter;
  char cstr[filename.str().size() + 2];
  strcpy(cstr, filename.str().c_str());

  std::ofstream outfile;
  outfile.open(cstr, std::ios::out);

  for (int i = 0; i < NI; i++)
  {
    outfile << i << "-" << TWall.value[i] << std::endl;
  }
  outfile.close();
}
