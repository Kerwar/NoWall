#include "Paralel.hpp"
#include <iostream>

using std::cout;
using std::endl;
using std::to_string;

Paralel::Paralel(int worldsize)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &worldMyProc);
  MPI_Comm_size(MPI_COMM_WORLD, &worldNProcs);
}

void Paralel::setUpComm(const int &NX, const int &NY)
{
  PROFILE_FUNCTION();

  cUSize = ceil((double)worldNProcs / 2);
  cDSize = worldNProcs - cUSize;

  // Associate the processor with the correspondant location
  if (worldMyProc < cUSize)
  {
    loc = up;
  }
  else
  {
    loc = down;
  }

  MPI_Comm helpComm;
  MPI_Comm_split(MPI_COMM_WORLD, loc, worldMyProc, &helpComm);
  MPI_Comm_size(helpComm, &nProcs);

  double nmRelation = sqrt((float)NX / (float)NY * nProcs);

  int dims[2] = {0, 1};

  nProcsInRow = SetXProcs(nmRelation);
  nProcsInCol = nProcs / nProcsInRow;

  nProcsInRow = nProcs;
  nProcsInCol = 1;

  dims[0] = nProcsInRow;
  dims[1] = nProcsInCol;

  int period[2] = {false, false};
  int reorder = true;

  MPI_Cart_create(helpComm, 2, dims, period, reorder, &myComm);
  MPI_Comm_free(&helpComm);

  MPI_Comm_rank(myComm, &myProc);

  int myCoords[2];
  MPI_Cart_coords(myComm, myProc, 2, myCoords);

  myRowId = myCoords[0];
  myColId = myCoords[1];

  MPI_Cart_shift(myComm, 0, 1, &myLeft, &myRight);
  MPI_Cart_shift(myComm, 1, 1, &myBot, &myTop);

  MainsProcs(up, UCMainProc);
  MainsProcs(down, DCMainProc);
}

void Paralel::setUpMesh(const int &NX, const int &NY, double &channelWidth, int &exi1, int &exi2)
{
  PROFILE_FUNCTION();
  NYChannel = floor(NY / 2.0);

  exI1 = exi1;
  exI2 = exi2;

  NXWall = exI2 - exI1;

  switch (loc)
  {
  case up:
    setProcMesh(NX, NYChannel);

    locIStr = 0;
    locIEnd = NX;
    locJStr = NYChannel;
    locJEnd = NY - 1;

    for (int i = 0; i < cUSize; i++)
    {
      jStrAllProc[i] += locJStr;
      jEndAllProc[i] += locJStr;
    }
    break;

  case down:
    setProcMesh(NX, NYChannel);

    locIStr = 0;
    locIEnd = NX;
    locJStr = 0;
    locJEnd = NYChannel - 1;
    break;
  }

  iStr = iStrAllProc[myProc];
  iEnd = iEndAllProc[myProc];
  jStr = jStrAllProc[myProc];
  jEnd = jEndAllProc[myProc];

  myNx = iEnd - iStr;
  myNy = jEnd - jStr;

  MPI_Type_vector(myNy + 2, 1, myNx + 2, MPI_DOUBLE_PRECISION, &column_type);
  MPI_Type_commit(&column_type);
}

Paralel::~Paralel()
{
  MPI_Type_free(&column_type);
}

void Paralel::freeComm()
{
  if (myComm != MPI_COMM_NULL)
    MPI_Comm_free(&myComm);
}

int Paralel::id(int x, int y)
{
  int myCoords[2] = {x, y};
  int result;

  MPI_Cart_rank(myComm, myCoords, &result);

  return result;
}

int Paralel::worldid(int row, int col)
{
  return col + row * nProcsInCol;
}

void Paralel::SendInfoToNeighbours(Field &vec)
{
  int NI = vec.NI;
  int NJ = vec.NJ;

  int leftTag = 10, rightTag = 11, topTag = 12, botTag = 13;

  leftTag *= (1 + loc);
  rightTag *= (1 + loc);
  topTag *= (1 + loc);
  botTag *= (1 + loc);

  MPI_Barrier(myComm);
  // SENDING INFO TO THE LEFT

  MPI_Sendrecv(&vec.value[1], 1, column_type, myLeft, leftTag,
               &vec.value[NI - 1], 1, column_type, myRight, leftTag, myComm, MPI_STATUS_IGNORE);

  if (!isProcNull(myRight))
    for (int j = 0; j < NJ; j++)
      vec.value[NI - 1 + j * NI] = (vec.value[NI - 1 + j * NI] + vec.value[NI - 2 + j * NI]) * 0.5;

  MPI_Barrier(myComm);
  // SENDING INFO TO THE RIGHT

  MPI_Sendrecv(&vec.value[NI - 1], 1, column_type, myRight, rightTag,
               &vec.value[0], 1, column_type, myLeft, rightTag, myComm, MPI_STATUS_IGNORE);

  MPI_Barrier(myComm);

  // SENDING INFO TO THE BOT
  MPI_Sendrecv(&vec.value[NI], NI, MPI_DOUBLE_PRECISION, myBot, botTag,
               &vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop, botTag, myComm, MPI_STATUS_IGNORE);

  if (!isProcNull(myTop))
    for (int i = 0; i < NI; i++)
      vec.value[i + (NJ - 1) * NI] = (vec.value[i + (NJ - 1) * NI] + vec.value[i + (NJ - 2) * NI]) * 0.5;
  MPI_Barrier(myComm);

  // SENDING INFO TO THE TOP
  MPI_Sendrecv(&vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop, topTag,
               &vec.value[0], NI, MPI_DOUBLE_PRECISION, myBot, topTag, myComm, MPI_STATUS_IGNORE);
}

void Paralel::MainsProcs(Loc location, int &proc)
{
  PROFILE_FUNCTION();
  int buffer;

  if (location == loc)
  {
    if (myProc == 0)
    {
      buffer = worldMyProc;
    }
    PMPI_Bcast(&buffer, 1, MPI_INT, 0, myComm);
  }

  int sendProc;
  switch (location)
  {
  case up:
    sendProc = 0;
    break;

  case down:
    sendProc = cUSize;
    break;
  }
  PMPI_Bcast(&buffer, 1, MPI_INT, sendProc, MPI_COMM_WORLD);

  proc = buffer;
}

void Paralel::ExchangeWallTemperature(Field &TWall, Field &TNextToWall, double &exCte, int &solExI1, int &solExI2)
{
  PROFILE_FUNCTION();

  int NI = TWall.NI;
  double wallRecv[NI] = {};
  double nextToWallRecv[NI] = {};

  int wallID = 22;
  int nextToWallID = 23;

  if (isLeftToRight())
  {
    PMPI_Sendrecv(&TWall.value[0], NI, MPI_DOUBLE_PRECISION, DCMainProc, wallID, &wallRecv, NI, MPI_DOUBLE_PRECISION, DCMainProc, wallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    PMPI_Sendrecv(&TNextToWall.value[0], NI, MPI_DOUBLE_PRECISION, DCMainProc, nextToWallID, &nextToWallRecv, NI, MPI_DOUBLE_PRECISION, DCMainProc, nextToWallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = solExI1; i < solExI2; i++)
      TWall.value[i] = ((8.0 + exCte) * (9.0 * TWall.value[i] - TNextToWall.value[i]) + 9.0 * exCte * wallRecv[i] - exCte * nextToWallRecv[i]) / (64.0 + 16 * exCte);
  }
  else if (isRightToLeft())
  {
    PMPI_Sendrecv(&TWall.value[0], NI, MPI_DOUBLE_PRECISION, UCMainProc, wallID, &wallRecv, NI, MPI_DOUBLE_PRECISION, UCMainProc, wallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    PMPI_Sendrecv(&TNextToWall.value[0], NI, MPI_DOUBLE_PRECISION, UCMainProc, nextToWallID, &nextToWallRecv, NI, MPI_DOUBLE_PRECISION, UCMainProc, nextToWallID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = solExI1; i < solExI2; i++)
      TWall.value[i] = ((8.0 + exCte) * (9.0 * TWall.value[i] - TNextToWall.value[i]) + 9.0 * exCte * wallRecv[i] - exCte * nextToWallRecv[i]) / (64.0 + 16 * exCte);
  }
}

void Paralel::ShareWallTemperatureInfo(Field &TWall, Field &T)
{
  PROFILE_FUNCTION();

  int NI = T.NI;
  int NJ = T.NJ;
  int n = T.NI - 1;

  switch (loc)
  {
  case up:
    SendWallTInTheChannel(TWall);
    for (int i = 1; i < n; i++)
      T.value[i] = TWall.value[i + iStr];
    break;

  case down:
    SendWallTInTheChannel(TWall);
    for (int i = 1; i < n; i++)
      T.value[i + (NJ - 1) * NI] = TWall.value[i + iStr];
    break;
  }
}

void Paralel::SendWallTInTheChannel(Field &TWall)
{
  PROFILE_FUNCTION();
  int NI = TWall.NI;

  PMPI_Bcast(&TWall.value[0], NI, MPI_DOUBLE, 0, myComm);
}

double Paralel::TWDinBorder(double (&T)[4], double &exCte)
{
  double newT = 9.0 * (exCte * T[2] + T[1]) - (T[0] + exCte * T[3]);
  newT = newT / (8.0 * (1.0 + exCte));

  return newT;
}

double Paralel::TWUinBorder(double (&T)[4], double &exCte)
{
  double newT = 9.0 * (T[2] + exCte * T[1]) - (exCte * T[0] + T[3]);
  newT = newT / (8.0 * (1.0 + exCte));

  return newT;
}

void Paralel::setProcMesh(const int &NX, const int &NY)
{
  PROFILE_FUNCTION();

  nProcsExtraXPoint = NX % nProcsInRow;
  nProcsExtraYPoint = NY % nProcsInCol;

  iStrAllProc = vector<int>(nProcs, 0);
  iEndAllProc = vector<int>(nProcs, 0);
  jStrAllProc = vector<int>(nProcs, 0);
  jEndAllProc = vector<int>(nProcs, 0);

  int xPointsPerProc = floor((float)NX / (float)nProcsInRow);
  int yPointsPerProc = floor((float)NY / (float)nProcsInCol);

  for (int RowId = 0; RowId < nProcsInRow; RowId++)
  {
    for (int ColId = 0; ColId < nProcsInCol; ColId++)
    {
      int idN = id(RowId, ColId);
      iStrAllProc[idN] = RowId * xPointsPerProc + std::min(nProcsExtraXPoint, RowId);
      iEndAllProc[idN] = (RowId + 1) * xPointsPerProc + std::min(nProcsExtraXPoint, RowId + 1);

      jStrAllProc[idN] = ColId * yPointsPerProc + std::min(nProcsExtraYPoint, ColId);
      jEndAllProc[idN] = (ColId + 1) * yPointsPerProc + std::min(nProcsExtraYPoint, ColId + 1);
      // std::cout << idN << " " << iStrAllProc[idN] << " " << iEndAllProc[idN] << std::endl;
      // std::cout << idN << " " << jStrAllProc[idN] << " " << jEndAllProc[idN] << std::endl;
    }
  }
}

void Paralel::SendInfoToCommMainProc(Field &vec, Field &sol)
{
  PROFILE_FUNCTION();

  int NI = vec.NI;
  int NJ = vec.NJ;

  int iIS = NI - 2;
  int jIS = NJ - 2;
  int sendCount = iIS * jIS;

  double sendBuffer[sendCount] = {};

  forAllInterior(NI, NJ)
  {
    sendBuffer[(i - 1) * jIS + (j - 1)] = vec.value[i + j * NI];
  }

  int recvCount = 0;
  int counts[nProcs] = {};
  int displacements[nProcs] = {};

  for (int k = 0; k < nProcs; k++)
  {
    displacements[k] = recvCount;
    counts[k] = (iEndAllProc[k] - iStrAllProc[k]) *
                (jEndAllProc[k] - jStrAllProc[k]);
    recvCount += counts[k];
  }
  double recvBuffer[recvCount];

  if (myProc == 0)
  {
    MPI_Gatherv(&sendBuffer[0], sendCount, MPI_DOUBLE, recvBuffer, counts, displacements,
                MPI_DOUBLE, 0, myComm);
  }
  else
  {
    MPI_Gatherv(&sendBuffer, sendCount, MPI_DOUBLE, NULL, NULL, NULL,
                MPI_DOUBLE, 0, myComm);
  }

  if (myProc == 0)
  {
    for (int k = 0; k < nProcs; k++)
    {
      int jSize = jEndAllProc[k] - jStrAllProc[k];
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++)
      {
        for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++)
        {
          int iRelative = i - iStrAllProc[k];
          int jRelative = j - jStrAllProc[k];

          int ij = displacements[k] + iRelative * jSize + jRelative;
          sol.value[i + j * sol.NI] = recvBuffer[ij];
        }
      }
    }
  }

  MPI_Barrier(myComm);
}

int Paralel::SetXProcs(double &nmRelation)
{
  PROFILE_FUNCTION();

  vector<int> divisors(nProcs, 1);

  for (int i = 1; i < nProcs; i++)
  {
    if (nProcs % (i + 1) == 0)
    {
      divisors[i] = i + 1;
    }
  }

  int result = divisors[0];
  for (long unsigned int i = 1; i < divisors.size(); i++)
  {
    if (std::abs(nmRelation - divisors[i]) <= std::abs(nmRelation - result))
    {
      result = divisors[i];
    }
  }
  return result;
}

string Paralel::PrintMyRank()
{
  // string result = to_string();
  return to_string(myProc) + "/" + to_string(nProcs - 1) +
         "-(" + to_string(myRowId) + "/" + to_string(nProcsInRow - 1) +
         "," + to_string(myColId) + "/" + to_string(nProcsInCol - 1) + ")" +
         "-" + to_string(worldMyProc) + "/" + to_string(worldNProcs - 1);
}

void Paralel::PrintProgress(string message, bool showProgress)
{
  if (showProgress == true)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    cout << message << PrintMyRank() << ".\n";
    if (worldMyProc == 0)
      cout << "------------------------------------------------------------\n";
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void Paralel::distributeToProcs(Field &sol, Field &vec)
{
  PROFILE_FUNCTION();

  int NI = vec.NI;

  int counts[nProcs];
  int sendCount = 0;
  int displacements[nProcs];

  double buffer[sol.NI * sol.NJ];

  for (int k = 0; k < nProcs; k++)
  {
    displacements[k] = sendCount;

    for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++)
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++)
      {
        int relative_i = i - iStrAllProc[k];
        int relative_j = j - jStrAllProc[k];
        int index = displacements[k] + relative_j * (iEndAllProc[k] - iStrAllProc[k]) + relative_i;
        buffer[index] = sol.value[i + j * sol.NI];
        // if(isWall()&&myProc == 0)
        // std::cout << index << " " << sol[i][j].value << std::endl;
        sendCount++;
      }
    counts[k] = sendCount - displacements[std::max(k, 0)];
    // std::cout << k << "·" << counts[k] << std::endl;
  }
  double recvValue[counts[myProc]];
  if (myProc == 0)
    MPI_Scatterv(buffer, counts, displacements, MPI_DOUBLE, &recvValue, counts[myProc], MPI_DOUBLE, 0, myComm);
  else
    MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvValue, counts[myProc], MPI_DOUBLE, 0, myComm);

  for (int i = 0; i < counts[myProc]; i++)
  {
    int vec_i = i % (NI - 2);
    int vec_j = floor(i / (NI - 2));
    // std::cout << i << " " << vec_i << " " << vec_j << " " << counts[myProc] << " " << recvValue[i]<< endl;
    vec.value[vec_i + 1 + (vec_j + 1) * NI] = recvValue[i];
  }
}

void Paralel::GatherWallTemperature(Field &TWall, Field &TNextToWall, Field &T, double &exCte)
{
  PROFILE_FUNCTION();
  // Field::vec1dfield TWall2Row(TWall.size());
  if (isLeftToRight())
  {
    GatherTemperature(T, TNextToWall, 2);
    GatherTemperature(T, TWall, 1);
    // TWall = TIntoDerivative(TWall, TWall2Row, 1, exCte);
  }
  else if (isRightToLeft())
  {
    GatherTemperature(T, TNextToWall, T.NJ - 3);
    GatherTemperature(T, TWall, T.NJ - 2);
    // TWall = TIntoDerivative(TWall, TWall2Row, 1, exCte);
  }
}

void Paralel::GatherTemperature(Field &T, Field &TWall, int J)
{
  int count[nProcs];
  int displacement[nProcs];

  for (int i = 0; i < nProcs; i++)
  {
    count[i] = iEndAllProc[i] - iStrAllProc[i];
    displacement[i] = iStrAllProc[i];
  }

  if (myProc == 0)
    MPI_Gatherv(&T.value[1 + J * T.NI], T.NI - 2, MPI_DOUBLE, &TWall.value[1], count, displacement, MPI_DOUBLE, 0, myComm);
  else
    MPI_Gatherv(&T.value[1 + J * T.NI], T.NI - 2, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, myComm);

  for (int i = 0; i < exI1; i++)
    TWall.value[i] = 0;

  for (int i = exI2; i < TWall.NI; i++)
    TWall.value[i] = 0;
}

// Field::vec1dfield Paralel::TIntoDerivative(const Field::vec1dfield &TWall, const Field::vec1dfield &TW2R, double exCte, double ex2)
// {
//   Field::vec1dfield TW(TWall.size());
//   for (unsigned long int i = 0; i < TWall.size(); i++)
//   {
//     TW[i].value = exCte * (9.0 * TWall[i].value - TW2R[i].value) / (8.0 * (ex2 + exCte));
//   }

//   return TW;
// }