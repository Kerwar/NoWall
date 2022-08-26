#include "variable.hpp"

Variable::Variable()
    : solT(NTOTAL, MINPUT),
      solF(NTOTAL, MINPUT),
      solZ(NTOTAL, MINPUT),
      solU(NTOTAL, MINPUT),
      solV(NTOTAL, MINPUT),
      U(NI, NJ),
      V(NI, NJ),
      T(NI, NJ),
      F(NI, NJ),
      Z(NI, NJ),
      TWall(NTOTAL + 2, 1),
      TNextToWall(NTOTAL + 2, 1),
      DTinWall(NI, 0),
      massFluxE(NI, NJ),
      massFluxN(NI, NJ) {}

Variable::~Variable() {}

void Variable::passInfoSolutionGrid(const Grid &solGrid, double viscX, double viscY){
  double viscTx = viscX;
  double viscTy = viscY;

  double viscFx = viscX / LeF;
  double viscFy = viscY / LeF;

  double viscZx = viscX / LeZ;
  double viscZy = viscY / LeZ;

  DY = solGrid.DY;

  solU.getGridInfoPassed(solGrid, viscTx, viscTy);
  solV.getGridInfoPassed(solGrid, viscTx, viscTy);
  solT.getGridInfoPassed(solGrid, viscTx, viscTy);
  solF.getGridInfoPassed(solGrid, viscFx, viscFy);
  solZ.getGridInfoPassed(solGrid, viscZx, viscZy);
  TNextToWall.XC = solU.XC;
}
void Variable::passInfoMyGrid(const Grid &myGrid, double viscX,
                              double viscY) {
  PROFILE_FUNCTION();

  double viscTx = viscX;
  double viscTy = viscY;

  double viscFx = viscX / LeF;
  double viscFy = viscY / LeF;

  double viscZx = viscX / LeZ;
  double viscZy = viscY / LeZ;

  U.getGridInfoPassed(myGrid, viscX, viscY);
  V.getGridInfoPassed(myGrid, viscX, viscY);

  T.getGridInfoPassed(myGrid, viscTx, viscTy);
  F.getGridInfoPassed(myGrid, viscFx, viscFy);
  Z.getGridInfoPassed(myGrid, viscZx, viscZy);
}

void Variable::setChannelEquations(Equation *&Teqn, Equation *&Feqn,
                                   Equation *&Zeqn,
                                   const Paralel<NTOTAL, MINPUT, NPROCS> &paralel,
                                   const double &m, const double &q, double &DT,
                                   const double &exCte, const int &iter) {
  PROFILE_FUNCTION();
  if (iter == 1) {
    if (paralel.isLeftToRight())
      Teqn = new Equation(fvm::diffusiveTermS(T, DTinWall, exCte) +
                          fvm::convectiveTerm(T, massFluxE, massFluxN, m) +
                          fvm::heatProduction(T, Z, q));
    else
      Teqn = new Equation(fvm::diffusiveTermN(T, DTinWall, exCte) +
                          fvm::convectiveTerm(T, massFluxE, massFluxN, m) +
                          fvm::heatProduction(T, Z, q));

    Feqn =
        new Equation(fvm::diffusiveTerm(F) +
                     fvm::convectiveTerm(F, massFluxE, massFluxN, m) -
                     fvm::intermidiateReaction(F, Z, T, beta, gamma_reaction));

    Zeqn =
        new Equation(fvm::diffusiveTerm(Z) +
                     fvm::convectiveTerm(Z, massFluxE, massFluxN, m) +
                     fvm::intermidiateReaction(F, Z, T, beta, gamma_reaction) -
                     fvm::zComsumptium(Z));

    Teqn->DT = DT;
    Feqn->DT = DT;
    Zeqn->DT = DT;
  } else {
    if (paralel.isLeftToRight())
      Teqn->updateEquation(fvm::diffusiveTermS(T, DTinWall, exCte) +
                           fvm::convectiveTerm(T, massFluxE, massFluxN, m) +
                           fvm::heatProduction(T, Z, q));
    else
      Teqn->updateEquation(fvm::diffusiveTermN(T, DTinWall, exCte) +
                           fvm::convectiveTerm(T, massFluxE, massFluxN, m) +
                           fvm::heatProduction(T, Z, q));

    Feqn->updateEquation(
        fvm::diffusiveTerm(F) +
        fvm::convectiveTerm(F, massFluxE, massFluxN, m) -
        fvm::intermidiateReaction(F, Z, T, beta, gamma_reaction));

    Zeqn->updateEquation(
        fvm::diffusiveTerm(Z) +
        fvm::convectiveTerm(Z, massFluxE, massFluxN, m) +
        fvm::intermidiateReaction(F, Z, T, beta, gamma_reaction) -
        fvm::zComsumptium(Z));
  }
}

void Variable::setWallShear(Equation *&Teqn, Equation *&Feqn, Equation *&Zeqn,
                            Direction side) {
  PROFILE_FUNCTION();
  if (side == west || side == east) {
    Teqn->SetWallShearY(T, side);
    Feqn->SetWallShearY(F, side);
    Zeqn->SetWallShearY(Z, side);
  } else {
    Teqn->SetWallShearX(T, side);
    Feqn->SetWallShearX(F, side);
    Zeqn->SetWallShearX(Z, side);
  }
}

void Variable::setExchangeWallShear(Equation *&Teqn, Equation *&Feqn,
                                    Equation *&Zeqn, Direction side, int iStr,
                                    int iEnd, int mainexI1, int mainexI2,
                                    int exI1, int exI2) {
  PROFILE_FUNCTION();

  Teqn->SetWallShearTX(T, iStr, iEnd, mainexI1, mainexI2, exI1, exI2, side);
  Feqn->SetWallShearX(F, side);
  Zeqn->SetWallShearX(Z, side);
}

double Variable::solveEquations(Equation *&Teqn, Equation *&Feqn,
                                Equation *&Zeqn, double alpha, int &niter,
                                int &itersol, int changeIter) {
  PROFILE_FUNCTION();
  double error = 10;
  double newerror = 10;

  error = Teqn->solve(T, alpha, niter, itersol, changeIter);

  newerror = Feqn->solve(F, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  newerror = Zeqn->solve(Z, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  forAllInterior(NI, NJ) {
    F.value[i + j * NI] = std::max(F.value[i + j * NI], 0.0);
  }

  return error;
}

double Variable::calculateNewM(const double &alpha, const double &m,
                               const double &q) {
  PROFILE_FUNCTION();
  int index = ifix + jfix * NI;
  T.value[index] = 0.7;
  double diffusiveTermX =
      T.viscX[index] * T.Se[jfix] * (T.value[index + 1] - T.value[index]) /
          T.DXPtoE[ifix] -
      T.viscX[index - 1] * T.Se[jfix - 1] *
          (T.value[index] - T.value[index - 1]) / T.DXPtoE[ifix - 1];

  double diffusiveTermY =
      T.viscY[index] * T.Sn[ifix] * (T.value[index + NI] - T.value[index]) /
          T.DYPtoN[jfix] -
      T.viscY[index - NI] * T.Sn[ifix - 1] *
          (T.value[index] - T.value[index - NI]) / T.DYPtoN[jfix - 1];

  double convectiveTermX =
      massFluxE.value[index] *
          (T.value[index + 1] * T.FXE[ifix] + T.value[index] * T.FXP[ifix]) -
      massFluxE.value[index - 1] * (T.value[index] * T.FXE[ifix - 1] +
                                    T.value[index - 1] * T.FXP[ifix - 1]);

  double convectiveTermY =
      massFluxN.value[index] *
          (T.value[index + NI] * T.FYN[jfix] + T.value[index] * T.FYP[jfix]) -
      massFluxN.value[index - NI] * (T.value[index] * T.FYN[jfix - 1] +
                                     T.value[index - NI] * T.FYP[jfix - 1]);

  double reactionTerm = q * Z.value[index] * T.volume[index];

  double newM = alpha * (reactionTerm + diffusiveTermX + diffusiveTermY) /
                    (convectiveTermX + convectiveTermY) +
                (1 - alpha) * m;

  newM = std::min(m * upperBoundFactorm, newM);
  newM = std::max(m * lowerBoundFactorm, newM);
  newM = std::max(newM, 0.4);
  return newM;
}

double Variable::calculateNewQ(const double &alpha, const double &m,
                               const double &q) {
  PROFILE_FUNCTION();
  int index = ifix + jfix * NI;
  T.value[index] = 0.7;
  double diffusiveTermX =
      T.viscX[index] * T.Se[jfix] * (T.value[index + 1] - T.value[index]) /
          T.DXPtoE[ifix] -
      T.viscX[index - 1] * T.Se[jfix - 1] *
          (T.value[index] - T.value[index - 1]) / T.DXPtoE[ifix - 1];

  double diffusiveTermY =
      T.viscY[index] * T.Sn[ifix] * (T.value[index + NI] - T.value[index]) /
          T.DYPtoN[jfix] -
      T.viscY[index - NI] * T.Sn[ifix - 1] *
          (T.value[index] - T.value[index - NI]) / T.DYPtoN[jfix - 1];

  double convectiveTermX =
      massFluxE.value[index] *
          (T.value[index + 1] * T.FXE[ifix] + T.value[index] * T.FXP[ifix]) -
      massFluxE.value[index - 1] * (T.value[index] * T.FXE[ifix - 1] +
                                    T.value[index - 1] * T.FXP[ifix - 1]);

  double convectiveTermY =
      massFluxN.value[index] *
          (T.value[index + NI] * T.FYN[jfix] + T.value[index] * T.FYP[jfix]) -
      massFluxN.value[index - NI] * (T.value[index] * T.FYN[jfix - 1] +
                                     T.value[index - NI] * T.FYP[jfix - 1]);

  double reactionTerm = Z.value[index] * T.volume[index];
  if (Z.value[index] < 10e-12) reactionTerm = 10e-12 * T.volume[index];

  double newQ = alpha *
                    (m * (convectiveTermX + convectiveTermY) - diffusiveTermX -
                     diffusiveTermY) /
                    (reactionTerm) +
                (1 - alpha) * q;

  newQ = std::min(q * upperBoundFactorm, newQ);
  newQ = std::max(q * lowerBoundFactorm, newQ);
  // std::cout << "->" << T.NI << "\n";
  return newQ;
}

void Variable::setFixIndex(double xfix, double yfix) {
  PROFILE_FUNCTION();

  ifix = 0;
  jfix = 0;

  double error = 1000;
  double newerror = abs(T.XC[ifix] - xfix);

  while (newerror < error) {
    error = abs(T.XC[ifix] - xfix);
    ifix++;
    newerror = abs(T.XC[ifix] - xfix);
  }
  ifix--;

  error = 1000;
  newerror = abs(T.YC[jfix] - yfix);

  while (newerror < error) {
    error = abs(T.YC[jfix] - yfix);
    jfix++;
    newerror = abs(T.YC[jfix] - yfix);
  }
  jfix--;
}

void Variable::readFile(Paralel<NTOTAL, MINPUT, NPROCS> &paralel, int block,
                        double &m, double yMin, double yMax) {
  PROFILE_FUNCTION();

  FileReader fileread;

  double laminarm = m;
  if (paralel.isRightToLeft()) laminarm = -m;

  U.laminarFlow(laminarm, yMin, yMax);
  V.initializeInternalField(0);

  if (paralel.myProc == 0) {
    fileread.readField(initialsol, block, 1, solT, paralel.locIStr,
                       paralel.locIEnd);
    fileread.readField(initialsol, block, 2, solF, paralel.locIStr,
                       paralel.locIEnd);
    fileread.readField(initialsol, block, 3, solZ, paralel.locIStr,
                       paralel.locIEnd);
  }

  paralel.myComm.wait();

  paralel.distributeToProcs(solT, T);
  paralel.distributeToProcs(solF, F);
  paralel.distributeToProcs(solZ, Z);
}

bool Variable::isTOutEqualToQ(double q, const Paralel<NTOTAL, MINPUT, NPROCS> &paralel) {
  PROFILE_FUNCTION();

  bool result = false;
  double expectedTOut = q * M / 2.0;
  if (paralel.myProc == 0 && paralel.isRightToLeft()) {
    double TOut = 0;
    for (int j = 0; j < M; j++) TOut += solT.value[1 + j * NTOTAL];

    if (abs(TOut - expectedTOut) < 10e-3) result = true;
  } else if (paralel.myProc == 0 && paralel.isLeftToRight()) {
    double TOut = 0;
    for (int j = 0; j < M; j++) TOut += solT.value[NTOTAL - 2 + j * NTOTAL];

    if (abs(TOut - expectedTOut) < 10e-3) result = true;
  }

  bool buffer = result;

  MPI_Bcast(&buffer, 1, MPI_CXX_BOOL, paralel.UCMainProc, MPI_COMM_WORLD);

  result = buffer || result;

  MPI_Bcast(&buffer, 1, MPI_CXX_BOOL, paralel.DCMainProc, MPI_COMM_WORLD);

  result = buffer || result;

  return result;
}

void Variable::writeTInWall(Paralel<NTOTAL, MINPUT, NPROCS> &paralel,
                            const Grid &mainGrid, const Grid &myGrid,
                            int iter) {
  PROFILE_FUNCTION();

  string file = "TW";
  std::ostringstream filename;
  filename << file << "--" << paralel.worldMyProc << "--" << iter << ".txt";

  std::string trystr = filename.str();
  std::ofstream outfile;
  outfile.open(trystr, std::ios::out);
  outfile << myGrid.exI1 << " " << myGrid.exI2 << std::endl;
  for (int i = paralel.iStr; i < paralel.iEnd; i++) {
    if (paralel.isLeftToRight())
      outfile << i << " " << i - paralel.iStr + 1 << " " << mainGrid.XC[i + 1]
              << " " << myGrid.XC[i - paralel.iStr + 1] << " "
              << TWall.value[i + 1] << " " << DTinWall[i + 1] << " "
              << std::endl;
    else
      outfile << i << " " << i - paralel.iStr + 1 << " " << mainGrid.XC[i + 1]
              << " " << myGrid.XC[i - paralel.iStr + 1] << " "
              << TWall.value[i + 1] << " "
              << T.value[i - paralel.iStr + 1 + (NJ - 1) * NI] << " "
              << std::endl;
  }

  outfile.close();
}
