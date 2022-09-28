#include "problem.hpp"

Problem::Problem(const double &prevm) : maxit(10E4), iShow(10E3) {
  PROFILE_FUNCTION();
  readFromFile = true;
  alpha = 0.95;

  m = 2.0;
  q = prevm;
  xfix = 43.95;
  yfix = 1.25;
  xHS_U = 44;
  xHS_D = 76;
  alphaWall = 1;
  DT = 10e-3;
  itersol = 4;
  viscX = 1;
  viscY = a2;
}

Problem::~Problem() {}

void Problem::setUpProblem() {
  PROFILE_FUNCTION();
  mainGrid =
      Grid(NTOTAL, MTOTAL, channel_xmin, channel_xmax, y_bot_min, y_top_max);
  mainGrid.SetIEx(wall_xmin, wall_xmax);

  myProc = paralel.channel.rank();
  nProcs = paralel.channel.size();

  myXMin = mainGrid.X[paralel.iStr];
  myXMax = mainGrid.X[paralel.iEnd];

  if (paralel.is_left2right()) {
    myYMin = y_top_min;
    myYMax = y_top_max;
  } else {
    myYMin = y_bot_min;
    myYMax = y_bot_max;
  }

  myGrid = Grid(N, M, myXMin, myXMax, myYMin, myYMax);
  myGrid.SetIEx(wall_xmin, wall_xmax);
}

void Problem::initializeVariables() {
  PROFILE_FUNCTION();

  variables.passInfoSolutionGrid(mainGrid, viscX, viscY);
  variables.passInfoMyGrid(myGrid, viscX, viscY);

  paralel.setProcWithFixPoint(fixPointInThisProc());
  if (fixPointInThisProc()) variables.setFixIndex(xfix, yfix);

  if (!readFromFile) {
    if (paralel.is_left2right())
      variables.initialize_top_channel(m, q, xHS_U);
    else if (paralel.is_right2left())
      variables.initialize_bot_channel(m, q, xHS_D);

    readFromFile = true;
  } else {
    if (paralel.is_left2right())
      variables.readFile(paralel, 2, m, y_top_min, 2.0 * y_top_max - y_top_min);

    MPI_Barrier(MPI_COMM_WORLD);
    if (paralel.is_right2left())
      variables.readFile(paralel, 1, m, 2.0 * y_bot_min - y_bot_max, y_bot_max);

    variables.sendInfoToNeighbours(paralel);
  }
  //  If density non constant this interpolate the Velocities and the mass
  //  fluxes have to be calculated with them
  // Field::vectorField UE(fieldOper.interpolatedFieldEast(U, myGrid));
  // Field::vectorField UN(fieldOper.interpolatedFieldNorth(V, myGrid));

  // fieldOper.getGridInfoPassed(UE, myGrid, viscX, viscY);
  // fieldOper.getGridInfoPassed(UN, myGrid, viscX, viscY);

  variables.set_mass_fluxes(myGrid);
  if (paralel.is_left2right() && paralel.isFirstProcOfRow()) {
    variables.set_inlet_up_channel();
  } else if (paralel.is_right2left() && paralel.isLastProcOfRow()) {
    variables.set_inlet_bot_channel();
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

    if (myProc == 0 && paralel.is_right2left())
      fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.sol,
                            paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.is_left2right())
      fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.sol,
                            paralel.loc);

  } else {
    sufix = "";
    writefilename(sufix);
    variables.sendInfoToCommMainProc(paralel);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.is_right2left())
      fileWriter.WriteInter(prefix, sufix, -1, mainGrid, myGrid, variables.sol,
                            paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.is_left2right())
      fileWriter.WriteInter(prefix, sufix, -1, mainGrid, myGrid, variables.sol,
                            paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);
    std::string lastfileprefix = "";
    std::string lastfilesufix = "";

    if (myProc == 0 && paralel.is_right2left())
      fileWriter.WriteInter(lastfileprefix, lastfilesufix, -1, mainGrid, myGrid,
                            variables.sol, paralel.loc);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myProc == 0 && paralel.is_left2right())
      fileWriter.WriteInter(lastfileprefix, lastfilesufix, -1, mainGrid, myGrid,
                            variables.sol, paralel.loc);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

double Problem::mainIter(int i) {
  PROFILE_FUNCTION();

  double error = 1;

  variables.setChannelEquations(sys, paralel, m, q, DT, i);
  Direction side = west;

  if (paralel.isFirstProcOfRow() && paralel.is_right2left())
    variables.setWallShear(sys, side);
  else
    variables.setDirichlet(sys, side);

  side = east;

  if (paralel.isLastProcOfRow() && paralel.is_left2right())
    variables.setWallShear(sys, side);
  else
    variables.setDirichlet(sys, side);

  side = south;

  if (paralel.is_right2left())
    variables.setWallShear(sys, side);
  else if (paralel.is_left2right())
    variables.setExchangeWallShear(sys, side, paralel.iStr, paralel.iEnd,
                                   mainGrid.exI1, mainGrid.exI2, myGrid.exI1,
                                   myGrid.exI2);

  side = north;

  if (paralel.is_right2left())
    variables.setExchangeWallShear(sys, side, paralel.iStr, paralel.iEnd,
                                   mainGrid.exI1, mainGrid.exI2, myGrid.exI1,
                                   myGrid.exI2);
  else if (paralel.is_left2right())
    variables.setWallShear(sys, side);

  variables.assembleEquations(sys);

  sys.T->EqnName = "T-Eqn";
  sys.F->EqnName = "F-Eqn";
  sys.Z->EqnName = "Z-Eqn";

  error = variables.solveEquations(sys, alpha, i, itersol, 0 * iShow);

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

bool Problem::fixPointInThisProc() {
  return isFixPointXCoordinateInThisProc() && isFixPointYCoordinateInThisProc();
}

bool Problem::isFixPointXCoordinateInThisProc() {
  return myXMin <= xfix && xfix < myXMax;
}

bool Problem::isFixPointYCoordinateInThisProc() {
  return myYMin <= yfix && yfix < myYMax;
}

void Problem::writefilename(string &filename) {
  PROFILE_FUNCTION();
  std::ostringstream temp;

  temp << filename;
  
  temp << "_NxM-" << NTOTAL << "x" << MTOTAL;
  temp << "_Lxa-" << channel_xmax - channel_xmin << "x" << a;
  temp << "_Ex-" << wall_xmax-wall_xmin;
  temp << "_q-" << q << "_m-" << m << "_beta-" << beta_reaction;
  temp << "_LeFxZ-" << LeF << "x" << LeZ << "_";

  filename = temp.str();
}

void Problem::retrieveNandM(int &nxOut, int &nyOut, double &mOut) {
  PROFILE_FUNCTION();

  nxOut = N * 5;
  nyOut = M * 2;

  mOut = m;
}