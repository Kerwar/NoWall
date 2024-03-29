#include "variable.hpp"

VariableManager::VariableManager()
    : sol(NTOTAL, MTOTAL),
      var(NI, NJ),
      TWall(NTOTAL + 2, 1),
      TNextToWall(NTOTAL + 2, 1),
      DTinWall(NI, 0),
      massFluxE(NI, NJ),
      massFluxN(NI, NJ) {}

VariableManager::~VariableManager() {}

void VariableManager::passInfoSolutionGrid(const Grid &solGrid, double viscX,
                                           double viscY) {
  double viscTx = viscX;
  double viscTy = viscY;

  double viscFx = viscX / LeF;
  double viscFy = viscY / LeF;

  double viscZx = viscX / LeZ;
  double viscZy = viscY / LeZ;

  DY = solGrid.DY;

  sol.U.getGridInfoPassed(solGrid, viscTx, viscTy);
  sol.V.getGridInfoPassed(solGrid, viscTx, viscTy);
  sol.T.getGridInfoPassed(solGrid, viscTx, viscTy);
  sol.F.getGridInfoPassed(solGrid, viscFx, viscFy);
  sol.Z.getGridInfoPassed(solGrid, viscZx, viscZy);
  TNextToWall.XC = sol.U.XC;
}

void VariableManager::passInfoMyGrid(const Grid &myGrid, double viscX,
                                     double viscY) {
  PROFILE_FUNCTION();

  double viscTx = viscX;
  double viscTy = viscY;

  double viscFx = viscX / LeF;
  double viscFy = viscY / LeF;

  double viscZx = viscX / LeZ;
  double viscZy = viscY / LeZ;

  var.U.getGridInfoPassed(myGrid, viscX, viscY);
  var.V.getGridInfoPassed(myGrid, viscX, viscY);

  var.T.getGridInfoPassed(myGrid, viscTx, viscTy);
  var.F.getGridInfoPassed(myGrid, viscFx, viscFy);
  var.Z.getGridInfoPassed(myGrid, viscZx, viscZy);
}

void VariableManager::setChannelEquations(SystemOfEquations &sys,
                                          const Paralel &paralel,
                                          const double &m, const double &q,
                                          double &DT, const int &iter) {
  PROFILE_FUNCTION();
  if (iter == 1) {
    if (paralel.is_left2right())
      sys.T = new Equation(fvm::diffusiveTermS(var.T, DTinWall, EXCTE) +
                           fvm::convectiveTerm(var.T, massFluxE, massFluxN, m) +
                           fvm::heatProduction(var.T, var.Z, q));
    else
      sys.T = new Equation(fvm::diffusiveTermN(var.T, DTinWall, EXCTE) +
                           fvm::convectiveTerm(var.T, massFluxE, massFluxN, m) +
                           fvm::heatProduction(var.T, var.Z, q));

    sys.F =
        new Equation(fvm::diffusiveTerm(var.F) +
                     fvm::convectiveTerm(var.F, massFluxE, massFluxN, m) -
                     fvm::intermidiateReaction(var.F, var.Z, var.T,
                                               beta_reaction, gamma_reaction));

    sys.Z =
        new Equation(fvm::diffusiveTerm(var.Z) +
                     fvm::convectiveTerm(var.Z, massFluxE, massFluxN, m) +
                     fvm::intermidiateReaction(var.F, var.Z, var.T,
                                               beta_reaction, gamma_reaction) -
                     fvm::zComsumptium(var.Z));

    sys.T->DT = DT;
    sys.F->DT = DT;
    sys.Z->DT = DT;
  } else {
    if (paralel.is_left2right())
      sys.T->updateEquation(
          fvm::diffusiveTermS(var.T, DTinWall, EXCTE) +
          fvm::convectiveTerm(var.T, massFluxE, massFluxN, m) +
          fvm::heatProduction(var.T, var.Z, q));
    else
      sys.T->updateEquation(
          fvm::diffusiveTermN(var.T, DTinWall, EXCTE) +
          fvm::convectiveTerm(var.T, massFluxE, massFluxN, m) +
          fvm::heatProduction(var.T, var.Z, q));

    sys.F->updateEquation(fvm::diffusiveTerm(var.F) +
                          fvm::convectiveTerm(var.F, massFluxE, massFluxN, m) -
                          fvm::intermidiateReaction(var.F, var.Z, var.T,
                                                    beta_reaction,
                                                    gamma_reaction));

    sys.Z->updateEquation(fvm::diffusiveTerm(var.Z) +
                          fvm::convectiveTerm(var.Z, massFluxE, massFluxN, m) +
                          fvm::intermidiateReaction(var.F, var.Z, var.T,
                                                    beta_reaction,
                                                    gamma_reaction) -
                          fvm::zComsumptium(var.Z));
  }
}

void VariableManager::setWallShear(SystemOfEquations &sys, Direction side) {
  PROFILE_FUNCTION();
  if (side == west || side == east) {
    sys.T->SetWallShearY(var.T, side);
    sys.F->SetWallShearY(var.F, side);
    sys.Z->SetWallShearY(var.Z, side);
  } else {
    sys.T->SetWallShearX(var.T, side);
    sys.F->SetWallShearX(var.F, side);
    sys.Z->SetWallShearX(var.Z, side);
  }
}

void VariableManager::setExchangeWallShear(SystemOfEquations &sys,
                                           Direction side)  {
  PROFILE_FUNCTION();

  sys.T->SetWallShearTX(var.T, end_ex_zone_1, begin_ex_zone_2, side);
  sys.F->SetWallShearX(var.F, side);
  sys.Z->SetWallShearX(var.Z, side);
}

double VariableManager::solveEquations(SystemOfEquations &sys, double alpha,
                                       int &niter, int &itersol,
                                       int changeIter) {
  PROFILE_FUNCTION();
  double error = 10;
  double newerror = 10;

  error = sys.T->solve(var.T, alpha, niter, itersol, changeIter);

  newerror = sys.F->solve(var.F, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  newerror = sys.Z->solve(var.Z, alpha, niter, itersol, changeIter);
  error = std::max(error, newerror);

  forAllInterior(NI, NJ) {
    var.F.value[i + j * NI] = std::max(var.F.value[i + j * NI], 0.0);
  }

  return error;
}

double VariableManager::calculateNewM(const double &alpha, const double &m,
                                      const double &q) {
  PROFILE_FUNCTION();
  int index = ifix + jfix * NI;

  Field *T = &var.T;

  T->value[index] = 0.7;
  double diffusiveTermX =
      T->viscX[index] * T->Se[jfix] * (T->value[index + 1] - T->value[index]) /
          T->DXPtoE[ifix] -
      T->viscX[index - 1] * T->Se[jfix - 1] *
          (T->value[index] - T->value[index - 1]) / T->DXPtoE[ifix - 1];

  double diffusiveTermY =
      T->viscY[index] * T->Sn[ifix] * (T->value[index + NI] - T->value[index]) /
          T->DYPtoN[jfix] -
      T->viscY[index - NI] * T->Sn[ifix - 1] *
          (T->value[index] - T->value[index - NI]) / T->DYPtoN[jfix - 1];

  double convectiveTermX =
      massFluxE.value[index] * (T->value[index + 1] * T->FXE[ifix] +
                                T->value[index] * T->FXP[ifix]) -
      massFluxE.value[index - 1] * (T->value[index] * T->FXE[ifix - 1] +
                                    T->value[index - 1] * T->FXP[ifix - 1]);

  double convectiveTermY =
      massFluxN.value[index] * (T->value[index + NI] * T->FYN[jfix] +
                                T->value[index] * T->FYP[jfix]) -
      massFluxN.value[index - NI] * (T->value[index] * T->FYN[jfix - 1] +
                                     T->value[index - NI] * T->FYP[jfix - 1]);

  double reactionTerm = q * var.Z.value[index] * T->volume[index];

  double newM = alpha * (reactionTerm + diffusiveTermX + diffusiveTermY) /
                    (convectiveTermX + convectiveTermY) +
                (1 - alpha) * m;

  newM = std::min(m * upperBoundFactorm, newM);
  newM = std::max(m * lowerBoundFactorm, newM);
  newM = std::max(newM, 0.4);
  return newM;
}

double VariableManager::calculateNewQ(const double &alpha, const double &m,
                                      const double &q) {
  PROFILE_FUNCTION();
  int index = ifix + jfix * NI;

  Field *T = &var.T;

  T->value[index] = 0.7;
  double diffusiveTermX =
      T->viscX[index] * T->Se[jfix] * (T->value[index + 1] - T->value[index]) /
          T->DXPtoE[ifix] -
      T->viscX[index - 1] * T->Se[jfix] *
          (T->value[index] - T->value[index - 1]) / T->DXPtoE[ifix - 1];

  double diffusiveTermY =
      T->viscY[index] * T->Sn[ifix] * (T->value[index + NI] - T->value[index]) /
          T->DYPtoN[jfix] -
      T->viscY[index - NI] * T->Sn[ifix] *
          (T->value[index] - T->value[index - NI]) / T->DYPtoN[jfix - 1];

  double convectiveTermX =
      massFluxE.value[index] * (T->value[index + 1] * T->FXE[ifix] +
                                T->value[index] * T->FXP[ifix]) -
      massFluxE.value[index - 1] * (T->value[index] * T->FXE[ifix - 1] +
                                    T->value[index - 1] * T->FXP[ifix - 1]);

  double convectiveTermY =
      massFluxN.value[index] * (T->value[index + NI] * T->FYN[jfix] +
                                T->value[index] * T->FYP[jfix]) -
      massFluxN.value[index - NI] * (T->value[index] * T->FYN[jfix - 1] +
                                     T->value[index - NI] * T->FYP[jfix - 1]);

  double reactionTerm = var.Z.value[index] * T->volume[index];
  if (var.Z.value[index] < 10e-10) reactionTerm = 10e-10 * T->volume[index];

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

void VariableManager::setFixIndex(double xfix, double yfix) {
  PROFILE_FUNCTION();

  ifix = 0;
  jfix = 0;

  double error = 1000;
  double newerror = abs(var.T.XC[ifix] - xfix);

  while (newerror < error && ifix < NI) {
    error = abs(var.T.XC[ifix] - xfix);
    ifix++;
    newerror = abs(var.T.XC[ifix] - xfix);
  }
  ifix--;

  error = 1000;
  newerror = abs(var.T.YC[jfix] - yfix);

  while (newerror < error) {
    error = abs(var.T.YC[jfix] - yfix);
    jfix++;
    newerror = abs(var.T.YC[jfix] - yfix);
  }
  jfix--;
}

void VariableManager::readFile(Paralel &paralel, int block, double &m,
                               double yMin, double yMax) {
  PROFILE_FUNCTION();

  FileReader fileread;

  double laminarm = m;
  if (paralel.is_right2left()) laminarm = -m;

  var.U.laminarFlow(laminarm, yMin, yMax);
  var.V.initializeInternalField(0);

  if (paralel.channel.is_main_proc()) {
    fileread.read_field(initialsol, block, 1, sol.T, paralel.locIStr,
                        paralel.locIEnd);
    fileread.read_field(initialsol, block, 2, sol.F, paralel.locIStr,
                        paralel.locIEnd);
    fileread.read_field(initialsol, block, 3, sol.Z, paralel.locIStr,
                        paralel.locIEnd);
  }

  paralel.channel.wait();

  paralel.distributeToProcs(sol.T, var.T);
  paralel.distributeToProcs(sol.F, var.F);
  paralel.distributeToProcs(sol.Z, var.Z);
}

bool VariableManager::isTOutEqualToQ(double q, const Paralel &paralel) {
  PROFILE_FUNCTION();

  bool result = false;
  double expectedTOut = q * M / 2.0;

  if (paralel.is_main_proc_down()) {
    double TOut = 0;
    for (int j = 0; j < M; j++) TOut += sol.T.value[1 + j * NTOTAL];

    if (abs(TOut - expectedTOut) < 10e-3) result = true;
  } else if (paralel.is_main_proc_up()) {
    double TOut = 0;
    for (int j = 0; j < M; j++) TOut += sol.T.value[NTOTAL - 2 + j * NTOTAL];

    if (abs(TOut - expectedTOut) < 10e-3) result = true;
  }

  bool buffer = result;

  MPI_Bcast(&buffer, 1, MPI_CXX_BOOL, paralel.UCMainProc, MPI_COMM_WORLD);

  result = buffer || result;

  MPI_Bcast(&buffer, 1, MPI_CXX_BOOL, paralel.DCMainProc, MPI_COMM_WORLD);

  result = buffer || result;

  return result;
}

void VariableManager::exchangeTemperature(const Paralel &paralel,
                                          const Grid &mainGrid) {
  PROFILE_FUNCTION();

  paralel.GatherWallTemperature(TWall, var.T);

  for (int i = 0; i < mainGrid.exI1; i++) TWall.value[i] = 0;
  for (int i = mainGrid.exI2; i < NTOTAL + 2; i++) TWall.value[i] = 0;

  if (paralel.channel.is_main_proc())
    TWall.value = paralel.ExchangeWallTemperature(TWall);
  MPI_Barrier(MPI_COMM_WORLD);

  paralel.ShareWallTemperatureInfo(TWall, var.T);

  // paralel.scatter(TWall.value[1], NI - 2, DTinWall[1], 0);
  MPI_Scatter(&TWall.value[1], NI - 2, MPI_DOUBLE, &DTinWall[1], NI - 2,
              MPI_DOUBLE, 0, paralel.channel.comm());
}

void VariableManager::set_mass_fluxes(const Grid &myGrid) {
  PROFILE_FUNCTION();

  massFluxE.getGridInfoPassed(myGrid, var.U.viscX[0], var.U.viscY[0]);
  massFluxN.getGridInfoPassed(myGrid, var.V.viscX[0], var.V.viscY[0]);

  massFluxE.computeEastMassFluxes(var.U);
  massFluxN.computeNorthMassFluxes(var.V);
}

void VariableManager::set_inlet_up_channel() {
  PROFILE_FUNCTION();
  var.T.inletBoundaryCondition(west, 0.0);
  var.F.inletBoundaryCondition(west, 1.0);
  var.Z.inletBoundaryCondition(west, 0.0);
}

void VariableManager::set_inlet_bot_channel() {
  PROFILE_FUNCTION();
  var.T.inletBoundaryCondition(east, 0.0);
  var.F.inletBoundaryCondition(east, 1.0);
  var.Z.inletBoundaryCondition(east, 0.0);
}

void VariableManager::assembleEquations(SystemOfEquations &sys) {
  PROFILE_FUNCTION();

  sys.T->assembleEquation();
  sys.F->assembleEquation();
  sys.Z->assembleEquation();

  sys.T->relax(var.T);
  sys.F->relax(var.F);
  sys.Z->relax(var.Z);
}

void VariableManager::sendInfoToCommMainProc(Paralel &paralel) {
  PROFILE_FUNCTION();
  paralel.SendInfoToCommMainProc(var.T, sol.T);
  paralel.SendInfoToCommMainProc(var.F, sol.F);
  paralel.SendInfoToCommMainProc(var.Z, sol.Z);
  paralel.SendInfoToCommMainProc(var.U, sol.U);
  paralel.SendInfoToCommMainProc(var.V, sol.V);
}

void VariableManager::sendInfoToNeighbours(Paralel &paralel) {
  PROFILE_FUNCTION();
  paralel.SendInfoToNeighbours(var.T);
  paralel.SendInfoToNeighbours(var.F);
  paralel.SendInfoToNeighbours(var.Z);
}

void VariableManager::setDirichlet(SystemOfEquations &sys, Direction side) {
  PROFILE_FUNCTION();
  sys.T->SetDirichlet(var.T, side);
  sys.F->SetDirichlet(var.F, side);
  sys.Z->SetDirichlet(var.Z, side);
}

void VariableManager::initialize_top_channel(const double &m, const double &q,
                                             const double &xHS) {
  PROFILE_FUNCTION();
  var.U.laminarFlow(m, y_top_min, 2.0 * y_top_max - y_top_min);
  var.V.initializeInternalField(0);
  var.T.InitializeT(q, xHS, channel_xmin, channel_xmax);
  var.F.InitializeF(xHS, channel_xmin, channel_xmax);
  var.Z.InitializeZ(z0hs, r0hs, xHS, y_top_max);
}

void VariableManager::initialize_bot_channel(const double &m, const double &q,
                                             const double &xHS) {
  PROFILE_FUNCTION();
  var.U.laminarFlow(-m, 2.0 * y_bot_min - y_bot_max, y_bot_max);
  var.V.initializeInternalField(0);
  var.T.InitializeT(q, xHS, channel_xmax, channel_xmin);
  var.F.InitializeF(xHS, channel_xmax, channel_xmin);
  var.Z.InitializeZ(z0hs, r0hs, xHS, y_bot_min);
}

void VariableManager::setExchangeIndexes(const Grid &mainGrid, const Grid &myGrid, const Paralel &paralel){
  if(paralel.iEnd <= mainGrid.exI1 || paralel.iStr >= mainGrid.exI2){
    end_ex_zone_1 = 1;
    begin_ex_zone_2 = NI;
  } 
  else if (paralel.iStr < mainGrid.exI1 && paralel.iEnd <= mainGrid.exI2){
    end_ex_zone_1 = 1;
    begin_ex_zone_2 = myGrid.exI1;
  } 
  else if (paralel.iStr >= mainGrid.exI1 && paralel.iEnd > mainGrid.exI2){
    end_ex_zone_1 = myGrid.exI2;
    begin_ex_zone_2 = NI;
  } 
  else if (paralel.iStr < mainGrid.exI1 && paralel.iEnd >= mainGrid.exI2) {
    end_ex_zone_1 = 1;
    begin_ex_zone_2 = NI;
  } 
}
