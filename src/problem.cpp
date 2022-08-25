#include "problem.hpp"

Problem::Problem(const double &prevm) : maxit(10E4), iShow(10E3), m(prevm) {
  PROFILE_FUNCTION();
  readFromFile = false;
  alpha = 0.95;

  q = 0.8;
  yWall = 0.0;
  yChannel = 0.5;
  yMin = 0;
  yMax = 1;
  xfix = 43.95;
  yfix = 0.75;
  xHS_U = 44;
  xHS_D = 76;
  z0hs = 0.5;
  alphaWall = 1;
  DT = 10e-3;
  itersol = 4;
  bK = 0.15;
  viscX = 1;
  viscY = a2;
}

Problem::~Problem() {}

void Problem::setUpProblem() {
  PROFILE_FUNCTION();
  mainGrid = Grid(NTOTAL, MINPUT, channel_xmin, channel_xmax, yMin, yMax);
  mainGrid.SetIEx(wall_xmin, wall_xmax);

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

  myGrid = Grid(N, M, myXMin, myXMax, myYMin, myYMax);
  myGrid.SetIEx(wall_xmin, wall_xmax);
}

void Problem::initializeVariables() {
  PROFILE_FUNCTION();
  variables.passInfoGridToAll(mainGrid, myGrid, viscX, viscY);
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
                                      r0hs, z0hs, channel_xmin, channel_xmax);
    else if (paralel.isRightToLeft())
      variables.initializeRightToLeft(m, myYMin - yChannel, myYMax, q, xHS_D,
                                      r0hs, z0hs, channel_xmin, channel_xmax);

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

void Problem::writeSolution(string &prefix, int i) {
  PROFILE_FUNCTION();

  FileWriter fileWriter;

  variables.sendInfoToCommMainProc(paralel);

  if (i != -1) {
    sufix = "";
    writefilename(sufix);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.isRightToLeft())
      fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.solU,
                            variables.solV, variables.solT, variables.solF,
                            variables.solZ, paralel.locIStr, paralel.locIEnd,
                            paralel.locJStr, paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.isLeftToRight())
      fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.solU,
                            variables.solV, variables.solT, variables.solF,
                            variables.solZ, paralel.locIStr, paralel.locIEnd,
                            paralel.locJStr, paralel.loc);

    // if (myProc == 0 && paralel.isRightToLeft())
    //   fileWriter.WriteTec(prefix, sufix, i, mainGrid, myGrid, variables.solU,
    //                       variables.solV, variables.solT, variables.solF,
    //                       variables.solZ, paralel.locIStr, paralel.locIEnd,
    //                       paralel.locJStr, paralel.loc);

    // MPI_Barrier(MPI_COMM_WORLD);

    // if (myProc == 0 && paralel.isLeftToRight())
    //   fileWriter.WriteTec(prefix, sufix, i, mainGrid, myGrid, variables.solU,
    //                       variables.solV, variables.solT, variables.solF,
    //                       variables.solZ, paralel.locIStr, paralel.locIEnd,
    //                       paralel.locJStr, paralel.loc);
  } else {
    sufix = "";
    writefilename(sufix);
    variables.writeTInWall(paralel, mainGrid, myGrid, 0);
    variables.sendInfoToCommMainProc(paralel);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.isRightToLeft())
      fileWriter.WriteInter(prefix, sufix, -1, mainGrid, myGrid, variables.solU,
                            variables.solV, variables.solT, variables.solF,
                            variables.solZ, paralel.locIStr, paralel.locIEnd,
                            paralel.locJStr, paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.isLeftToRight())
      fileWriter.WriteInter(prefix, sufix, -1, mainGrid, myGrid, variables.solU,
                            variables.solV, variables.solT, variables.solF,
                            variables.solZ, paralel.locIStr, paralel.locIEnd,
                            paralel.locJStr, paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);
    std::string lastfileprefix = "";
    std::string lastfilesufix = "";

    if (myProc == 0 && paralel.isRightToLeft())
      fileWriter.WriteInter(lastfileprefix, lastfilesufix, -1, mainGrid, myGrid,
                            variables.solU, variables.solV, variables.solT,
                            variables.solF, variables.solZ, paralel.locIStr,
                            paralel.locIEnd, paralel.locJStr, paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.isLeftToRight())
      fileWriter.WriteInter(lastfileprefix, lastfilesufix, -1, mainGrid, myGrid,
                            variables.solU, variables.solV, variables.solT,
                            variables.solF, variables.solZ, paralel.locIStr,
                            paralel.locIEnd, paralel.locJStr, paralel.loc);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

double Problem::mainIter(int i) {
  PROFILE_FUNCTION();

  double error = 1;

  variables.setChannelEquations(Teqn, Feqn, Zeqn, paralel, m, q, DT, exCte, i);
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

  // if (fixPointInThisProc()) m = variables.calculateNewM(0.75 * alpha, m, q);
  if (fixPointInThisProc()) q = variables.calculateNewQ(0.75 * alpha, m, q);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&q, 1, MPI_DOUBLE, paralel.fixPointProc, MPI_COMM_WORLD);

  MPI_Allreduce(&error, &error_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  variables.sendInfoToNeighbours(paralel);

  variables.exchangeTemperature(paralel, mainGrid);

  MPI_Barrier(MPI_COMM_WORLD);
  return error_result;
}

void Problem::freeComm() { paralel.freeComm(); }

bool Problem::fixPointInThisProc() {
  return isFixPointXCoordinateInThisProc() && isFixPointYCoordinateInThisProc();
}

bool Problem::isFixPointXCoordinateInThisProc() {
  return myXMin <= xfix && xfix <= myXMax;
}

bool Problem::isFixPointYCoordinateInThisProc() {
  return myYMin <= yfix && yfix <= myYMax;
}

void Problem::writefilename(string &filename) {
  PROFILE_FUNCTION();
  std::ostringstream temp;

  filename.append("_NxM-");
  temp.str("");
  temp.clear();
  temp << NTOTAL;
  filename.append(temp.str());

  filename.append("x");
  temp.str("");
  temp.clear();
  temp << MINPUT;
  filename.append(temp.str());

  filename.append("_Lxa-");
  temp.str("");
  temp.clear();
  temp << channel_xmax - channel_xmin;
  filename.append(temp.str());

  filename.append("x");
  temp.str("");
  temp.clear();
  temp << a;
  filename.append(temp.str());

  filename.append("_Ex-");
  temp.str("");
  temp.clear();
  temp << wall_xmax - wall_xmin;
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

void Problem::retrieveNandM(int &nxOut, int &nyOut, double &mOut) {
  PROFILE_FUNCTION();

  nxOut = N * 5;
  nyOut = M * 2;

  mOut = m;
}

void Problem::setExchangeConstant() {
  exCte = bK * a * a / 2.0;
  exCte /= (1.0 + mainGrid.DY * exCte);
}