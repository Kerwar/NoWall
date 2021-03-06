#include "Problem.hpp"

template<int N, int M>
Problem<N,M>::Problem() : maxit(10), iShow(1), m(2), variables(1,1,1,1)
{
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
  xfix = 43;
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
  viscY = a2;
}

template<int N, int M>
Problem<N,M>::Problem(const int &ni, const int &nj, const double &prevm) :
 maxit(10E6), iShow(3E5), m(prevm), variables(ni, nj, N, M)
{
  PROFILE_FUNCTION();
  // if (nx > 20000)
  // {
	// tolerance = 10e-6;
	// maxit = 1000000;}
}

template<int N, int M>
Problem<N,M>::~Problem()
{
  // delete variables;
}

template<int N, int M>
void Problem<N,M>::setUpProblem()
{
  PROFILE_FUNCTION();
  paralel.setUpComm(N, M);
  mainGrid = Grid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);
  paralel.setUpMesh(N, M, mainGrid.exI1, mainGrid.exI2);

  myProc = paralel.myProc;
  nProcs = paralel.nProcs;

  myXMin = mainGrid.X[paralel.iStr];
  myXMax = mainGrid.X[paralel.iEnd];
  myYMin = 0.0;
  myYMax = yChannel;

  if (paralel.isLeftToRight())
  {
    myYMin += yChannel + yWall;
    myYMax += yChannel + yWall;
  }

  // std::cout << paralel.worldMyProc << " " << paralel.myNx << " " <<  paralel.myNy<< " " <<
  //  myXMin<< " " << myXMax<< " " << myYMin<< " " << myYMax << std::endl;
  myGrid = Grid(paralel.myNx, paralel.myNy, myXMin, myXMax, myYMin, myYMax);
  myGrid.SetIEx(xExMin, xExMax);

  solN = N;
  solM = M;
}

template<int N, int M>
void Problem<N,M>::initializeVariables()
{
  PROFILE_FUNCTION();
  variables.passInfoGridToAll(mainGrid, myGrid, viscX, viscY, LeF, LeZ);
  setExchangeConstant(myGrid.DY);
  int fixPointProcBuffer = 0;

  if (fixPointInThisProc())
  {
    fixPointProcBuffer += paralel.worldMyProc;
    variables.setFixIndex(xfix, yfix);
  }

  MPI_Allreduce(&fixPointProcBuffer, &paralel.fixPointProc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (!readFromFile)
  {
    if (paralel.isLeftToRight())
      variables.initializeLeftToRight(m, myYMin, myYMax + yChannel, q, xHS_U, r0hs, z0hs, xMin, xMax);
    else if (paralel.isRightToLeft())
      variables.initializeRightToLeft(m, myYMin - yChannel, myYMax, q, xHS_D, r0hs, z0hs, xMin, xMax);

    readFromFile = true;
  }
  else
  {
    if (paralel.isLeftToRight())
      variables.readFile(paralel, 2, m, myYMin, myYMax + yChannel);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (paralel.isRightToLeft())
      variables.readFile(paralel, 1, m, myYMin - yChannel, myYMax);

    variables.sendInfoToNeighbours(paralel);
  }
  //  If density non constant this interpolate the Velocities and the mass fluxes have to be calculated with them
  // Field::vectorField UE(fieldOper.interpolatedFieldEast(U, myGrid));
  // Field::vectorField UN(fieldOper.interpolatedFieldNorth(V, myGrid));

  // fieldOper.getGridInfoPassed(UE, myGrid, viscX, viscY);
  // fieldOper.getGridInfoPassed(UN, myGrid, viscX, viscY);

  variables.setMassFluxes(myGrid);

  if (paralel.isLeftToRight() && paralel.isProcNull(paralel.myLeft))
  {
    variables.setInletBoundaryConditionLeftToRight();
  }
  else if (paralel.isRightToLeft() && paralel.isProcNull(paralel.myRight))
  {
    variables.setInletBoundaryConditionRightToLeft();
  }

  variables.sendInfoToCommMainProc(paralel);

  variables.exchangeTemperature(paralel, exCte, mainGrid.exI1, mainGrid.exI2);
  
}

template<int N, int M>
void Problem<N,M>::writeSolution(string &prefix, int i)
{
  PROFILE_FUNCTION();

  FileWriter fileWriter;

  if (i == 0)
  {
    sufix = "";
    writefilename(sufix);
  }
  else if (i == -1)
  {
    sufix = "";
    variables.writeTInWall(paralel, mainGrid, myGrid, 0);
  }

  variables.sendInfoToCommMainProc(paralel);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myProc == 0 && paralel.isRightToLeft())
    fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.solU, variables.solV, variables.solT, variables.solF, variables.solZ, paralel.locIStr, paralel.locIEnd, paralel.locJStr, paralel.loc);
    // fileWriter.WriteBin(prefix, sufix, i, mainGrid, variables.solU, variables.solV, variables.solT, variables.solF, variables.solZ, paralel.locIStr, paralel.locIEnd, paralel.locJStr, paralel.locJEnd, paralel.loc);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myProc == 0 && paralel.isLeftToRight())
  fileWriter.WriteInter(prefix, sufix, i, mainGrid, myGrid, variables.solU, variables.solV, variables.solT, variables.solF, variables.solZ, paralel.locIStr, paralel.locIEnd, paralel.locJStr, paralel.loc);
    // fileWriter.WriteBin(prefix, sufix, i, mainGrid, variables.solU, variables.solV, variables.solT, variables.solF, variables.solZ, paralel.locIStr, paralel.locIEnd, paralel.locJStr, paralel.locJEnd, paralel.loc);
  MPI_Barrier(MPI_COMM_WORLD);

}

template<int N, int M>
double Problem<N,M>::mainIter(int i)
{
  PROFILE_FUNCTION();

  double error = 1;

  Equation *Teqn;
  Equation *Feqn;
  Equation *Zeqn;

  variables.setChannelEquations(Teqn, Feqn, Zeqn, m, q, beta, gamma, DT);

  Field::Direction side = Field::west;

  if (paralel.isFirstProcOfRow() && paralel.isRightToLeft())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);
  else
    variables.setDirichlet(Teqn, Feqn, Zeqn, side);

  side = Field::east;

  if (paralel.isLastProcOfRow() && paralel.isLeftToRight())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);
  else
    variables.setDirichlet(Teqn, Feqn, Zeqn, side);

  side = Field::south;

  if (paralel.isFirstProcOfCol() && paralel.isRightToLeft())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);
  else if (paralel.isFirstProcOfCol() && paralel.isLeftToRight())
    variables.setExchangeWallShear(Teqn, Feqn, Zeqn, side, paralel.iStr, paralel.iEnd, mainGrid.exI1, mainGrid.exI2,
                                    myGrid.exI1, myGrid.exI2);
  else
    variables.setDirichlet(Teqn, Feqn, Zeqn, side);

  side = Field::north;

  if (paralel.isLastProcOfCol() && paralel.isRightToLeft())
    variables.setExchangeWallShear(Teqn, Feqn, Zeqn, side, paralel.iStr, paralel.iEnd, mainGrid.exI1, mainGrid.exI2,
                                    myGrid.exI1, myGrid.exI2);
  else if (paralel.isLastProcOfCol() && paralel.isLeftToRight())
    variables.setWallShear(Teqn, Feqn, Zeqn, side);
  else
    variables.setDirichlet(Teqn, Feqn, Zeqn, side);

  variables.assembleEquations(Teqn, Feqn, Zeqn);

  Teqn->EqnName = "T-Eqn";
  Feqn->EqnName = "F-Eqn";
  Zeqn->EqnName = "Z-Eqn";

  error = variables.solveEquations(Teqn, Feqn, Zeqn, alpha, i, itersol, 10*iShow);

  double error_result = 0;

  if (fixPointInThisProc())
    m = variables.calculateNewM(0.75 * alpha, m, q);

  {
  PROFILE_SCOPE("Cast+Reduce");
  MPI_Bcast(&m, 1, MPI_DOUBLE, paralel.fixPointProc, MPI_COMM_WORLD);

  MPI_Allreduce(&error, &error_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  variables.sendInfoToNeighbours(paralel);
  variables.exchangeTemperature(paralel, exCte, mainGrid.exI1, mainGrid.exI2);

  delete Teqn;
  delete Feqn;
  delete Zeqn;

  MPI_Barrier(MPI_COMM_WORLD);
  return error_result;
}

template<int N, int M>
void Problem<N,M>::freeComm()
{
  paralel.freeComm();
}

template<int N, int M>
bool Problem<N,M>::fixPointInThisProc()
{
  PROFILE_FUNCTION();
  bool result = false;

  if (myXMin <= xfix && xfix <= myXMax)
    if (myYMin <= yfix && yfix <= myYMax)
      result = true;

  return result;
}

template<int N, int M>
void Problem<N,M>::writefilename(string &filename)
{
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

template<int N, int M>
void Problem<N,M>::retrieveNandM(int &nxOut, int &nyOut, double &mOut)
{
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

template<int N, int M>
void inline Problem<N,M>::setExchangeConstant(double &deltaY)
{
  exCte = bK * a * a * 1.5 * deltaY ;
}
