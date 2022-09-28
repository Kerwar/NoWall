#include "paralel.hpp"

Paralel::Paralel() {
  PROFILE_FUNCTION();

  MPI_Comm_rank(MPI_COMM_WORLD, &worldMyProc);

  if (worldMyProc < NPROCS / 2)
    loc = up;
  else
    loc = down;

  MPI_Comm helpComm;
  MPI_Comm_split(MPI_COMM_WORLD, loc, worldMyProc, &helpComm);
  int dims[2] = {NPROCSINROW, NPROCSINCOL};

  int period[2] = {false, false};
  int reorder = true;

  MPI_Cart_create(helpComm, 2, dims, period, reorder, channel.comm_ref());
  MPI_Comm_free(&helpComm);
  channel.update();

  int myCoords[2];
  MPI_Cart_coords(channel.comm(), channel.rank(), 2, myCoords);

  myRowId = myCoords[0];
  myColId = myCoords[1];

  MPI_Cart_shift(channel.comm(), 0, 1, &myLeft, &myRight);
  MPI_Cart_shift(channel.comm(), 1, 1, &myBot, &myTop);

  MainsProcs(up, UCMainProc);
  MainsProcs(down, DCMainProc);

  switch (loc) {
    case up:
      setProcMesh();

      locIStr = 0;
      locIEnd = NTOTAL;
      locJStr = NJ - 2;
      locJEnd = MTOTAL - 1;

      for (int i = 0; i < NPROCS / 2; i++) {
        jStrAllProc[i] += locJStr;
        jEndAllProc[i] += locJStr;
      }
      break;

    case down:
      setProcMesh();

      locIStr = 0;
      locIEnd = NTOTAL;
      locJStr = 0;
      locJEnd = NJ - 2;
      break;
  }

  iStr = iStrAllProc[channel.rank()];
  iEnd = iEndAllProc[channel.rank()];
  jStr = jStrAllProc[channel.rank()];
  jEnd = jEndAllProc[channel.rank()];

  MPI_Type_vector(NJ, 1, NI, MPI_DOUBLE_PRECISION, &column_type);
  MPI_Type_commit(&column_type);
}

Paralel::~Paralel() {}

int Paralel::id(int x, int y) {
  int myCoords[2] = {x, y};
  int result;

  MPI_Cart_rank(channel.comm(), myCoords, &result);

  return result;
}

void Paralel::SendInfoToNeighbours(Field &vec) {
  int leftTag = 10, rightTag = 11, topTag = 12, botTag = 13;

  leftTag *= (1 + loc);
  rightTag *= (1 + loc);
  topTag *= (1 + loc);
  botTag *= (1 + loc);

  channel.wait();
  // SENDING INFO TO THE LEFT
  MPI_Sendrecv(&vec.value[1], 1, column_type, myLeft, leftTag,
               &vec.value[NI - 1], 1, column_type, myRight, leftTag,
               channel.comm(), MPI_STATUS_IGNORE);

  if (!isProcNull(myRight))
    for (int j = 0; j < NJ; j++)
      vec.value[NI - 1 + j * NI] =
          (vec.value[NI - 1 + j * NI] + vec.value[NI - 2 + j * NI]) * 0.5;

  channel.wait();
  // SENDING INFO TO THE RIGHT

  MPI_Sendrecv(&vec.value[NI - 1], 1, column_type, myRight, rightTag,
               &vec.value[0], 1, column_type, myLeft, rightTag, channel.comm(),
               MPI_STATUS_IGNORE);

  channel.wait();

  // SENDING INFO TO THE BOT
  MPI_Sendrecv(&vec.value[NI], NI, MPI_DOUBLE_PRECISION, myBot, botTag,
               &vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop,
               botTag, channel.comm(), MPI_STATUS_IGNORE);

  if (!isProcNull(myTop))
    for (int i = 0; i < NI; i++)
      vec.value[i + (NJ - 1) * NI] =
          (vec.value[i + (NJ - 1) * NI] + vec.value[i + (NJ - 2) * NI]) * 0.5;
  channel.wait();

  // SENDING INFO TO THE TOP
  MPI_Sendrecv(&vec.value[(NJ - 1) * NI], NI, MPI_DOUBLE_PRECISION, myTop,
               topTag, &vec.value[0], NI, MPI_DOUBLE_PRECISION, myBot, topTag,
               channel.comm(), MPI_STATUS_IGNORE);
}

void Paralel::MainsProcs(Loc location, int &proc) {
  PROFILE_FUNCTION();
  int buffer;

  if (location == loc) {
    if (channel.rank() == 0) buffer = worldMyProc;
    PMPI_Bcast(&buffer, 1, MPI_INT, 0, channel.comm());
  }

  int sendProc;
  switch (location) {
    case up:
      sendProc = 0;
      break;

    case down:
      sendProc = NPROCS / 2;
      break;
  }
  PMPI_Bcast(&buffer, 1, MPI_INT, sendProc, MPI_COMM_WORLD);

  proc = buffer;
}

vector<double> Paralel::ExchangeWallTemperature(const Field &TWall) const {
  PROFILE_FUNCTION();

  int NW = NTOTAL + 2;

  vector<double> recv(NW, 0);
  int wallID = 22;

  if (is_left2right())
    PMPI_Sendrecv(&TWall.value[0], NW, MPI_DOUBLE_PRECISION, DCMainProc, wallID,
                  &recv[0], NW, MPI_DOUBLE_PRECISION, DCMainProc, wallID,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  else if (is_right2left())
    PMPI_Sendrecv(&TWall.value[0], NW, MPI_DOUBLE_PRECISION, UCMainProc, wallID,
                  &recv[0], NW, MPI_DOUBLE_PRECISION, UCMainProc, wallID,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return recv;
}

void Paralel::ShareWallTemperatureInfo(Field &TWall, Field &T) const {
  PROFILE_FUNCTION();

  int n = NI - 1;

  switch (loc) {
    case up:
      SendWallTInTheChannel(TWall);
      for (int i = 1; i < n; i++) T.value[i] = TWall.value[i + iStr];
      break;

    case down:
      SendWallTInTheChannel(TWall);
      for (int i = 1; i < n; i++)
        T.value[i + (NJ - 1) * NI] = TWall.value[i + iStr];
      break;
  }
}

void Paralel::SendWallTInTheChannel(Field &TWall) const {
  PROFILE_FUNCTION();

  PMPI_Bcast(&TWall.value[0], NTOTAL + 2, MPI_DOUBLE, 0, channel.comm());
}

void Paralel::setProcMesh() {
  PROFILE_FUNCTION();

  iStrAllProc = vector<int>(channel.size(), 0);
  iEndAllProc = vector<int>(channel.size(), 0);
  jStrAllProc = vector<int>(channel.size(), 0);
  jEndAllProc = vector<int>(channel.size(), 0);

  for (int RowId = 0; RowId < NPROCSINROW; RowId++)
    for (int ColId = 0; ColId < NPROCSINCOL; ColId++) {
      int idN = id(RowId, ColId);
      iStrAllProc[idN] = RowId * (NI - 2);
      iEndAllProc[idN] = (RowId + 1) * (NI - 2);

      jStrAllProc[idN] = ColId * (NJ - 2);
      jEndAllProc[idN] = (ColId + 1) * (NJ - 2);
    }
}

void Paralel::SendInfoToCommMainProc(Field &vec, Field &sol) {
  PROFILE_FUNCTION();

  const int iIS = NI - 2;
  const int jIS = NJ - 2;
  const int sendCount = iIS * jIS;

  double sendBuffer[sendCount];

  forAllInterior(NI, NJ) {
    sendBuffer[(i - 1) + (j - 1) * iIS] = vec.value[i + j * NI];
  }

  int recvCount = 0;
  int counts[NPROCS / 2] = {};
  int displacements[NPROCS / 2] = {};

  for (int k = 0; k < channel.size(); k++) {
    displacements[k] = recvCount;
    counts[k] = iIS * jIS;
    recvCount += counts[k];
  }

  std::vector<double> recvBuffer(recvCount);
  if (channel.is_main_proc())
    MPI_Gatherv(&sendBuffer[0], sendCount, MPI_DOUBLE, &recvBuffer[0], counts,
                displacements, MPI_DOUBLE, 0, channel.comm());
  else
    MPI_Gatherv(&sendBuffer, sendCount, MPI_DOUBLE, NULL, NULL, NULL,
                MPI_DOUBLE, 0, channel.comm());

  if (channel.is_main_proc())
    for (int k = 0; k < channel.size(); k++) {
      int iSize = iEndAllProc[k] - iStrAllProc[k];
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++)
        for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++) {
          int iRelative = i - iStrAllProc[k];
          int jRelative = j - jStrAllProc[k];

          int ij = displacements[k] + iRelative + jRelative * iSize;
          sol.value[i + j * NTOTAL] = recvBuffer[ij];
        }
    }

  channel.wait();
}

string Paralel::PrintMyRank() {
  // string result = to_string();
  return to_string(channel.rank()) + "/" + to_string(channel.size()) + "-(" +
         to_string(myRowId) + "/" + to_string(NPROCSINROW - 1) + "," +
         to_string(myColId) + "/" + to_string(NPROCSINCOL - 1) + ")" + "-" +
         to_string(worldMyProc) + "/" + to_string(NPROCS);
}

void Paralel::distributeToProcs(Field &sol, Field &vec) {
  PROFILE_FUNCTION();

  int counts[NPROCS / 2];
  int displacements[NPROCS / 2];
  int sendCount = 0;

  double buffer[NTOTAL * MINPUT];

  for (int k = 0; k < channel.size(); k++) {
    displacements[k] = sendCount;

    for (int j = jStrAllProc[k]; j < jEndAllProc[k]; j++)
      for (int i = iStrAllProc[k]; i < iEndAllProc[k]; i++) {
        int relative_i = i - iStrAllProc[k];
        int relative_j = j - jStrAllProc[k];
        int index = displacements[k] +
                    relative_j * (iEndAllProc[k] - iStrAllProc[k]) + relative_i;
        buffer[index] = sol.value[i + j * NTOTAL];
        // if(isWall()&&myProc == 0)
        // std::cout << index << " " << sol[i][j].value << std::endl;
        sendCount++;
      }
    counts[k] = sendCount - displacements[std::max(k, 0)];
    // std::cout << k << "·" << counts[k] << std::endl;
  }
  vector<double> recvValue(counts[channel.rank()]);
  if (channel.is_main_proc())
    MPI_Scatterv(buffer, counts, displacements, MPI_DOUBLE, &recvValue[0],
                 counts[channel.rank()], MPI_DOUBLE, 0, channel.comm());
  else
    MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, &recvValue[0],
                 counts[channel.rank()], MPI_DOUBLE, 0, channel.comm());

  for (int i = 0; i < counts[channel.rank()]; i++) {
    int vec_i = i % (NI - 2);
    int vec_j = floor(i / (NI - 2));
    // std::cout << i << " " << vec_i << " " << vec_j << " " << counts[myProc]
    // << " " << recvValue[i]<< endl;
    vec.value[vec_i + 1 + (vec_j + 1) * NI] = recvValue[i];
  }
}

void Paralel::GatherWallTemperature(Field &TWall, Field &T) const {
  PROFILE_FUNCTION();
  if (is_left2right())
    GatherTemperature(T, TWall, 1);
  else if (is_right2left())
    GatherTemperature(T, TWall, T.NJ - 2);
}

void Paralel::GatherTemperature(Field &T, Field &TWall, int J) const {
  MPI_Gather(&T.value[1 + J * NI], NI - 2, MPI_DOUBLE, &TWall.value[1], NI - 2,
             MPI_DOUBLE, 0, channel.comm());
  // for (int i = 0; i < exI1; i++) TWall.value[i] = 0;

  // for (int i = exI2; i < N + 2; i++) TWall.value[i] = 0;
}

void Paralel::setProcWithFixPoint(const bool &fixPointInThisProc) {
  if (fixPointInThisProc) fixPointProc += worldMyProc;

  MPI_Allreduce(&fixPointProc, &fixPointProc, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
}

bool Paralel::is_left2right() const { return loc == up ? true : false; }

bool Paralel::is_right2left() const { return loc == down ? true : false; }

bool Paralel::is_main_proc_up() const {
  return is_left2right() && channel.is_main_proc();
}

bool Paralel::is_main_proc_down() const {
  return is_right2left() && channel.is_main_proc();
}

bool Paralel::isFirstProcOfCol() const { return isProcNull(myBot); }

bool Paralel::isLastProcOfCol() const { return isProcNull(myTop); }

bool Paralel::isFirstProcOfRow() const { return isProcNull(myLeft); }

bool Paralel::isLastProcOfRow() const { return isProcNull(myRight); }

bool Paralel::isProcNull(const int &proc) const {
  return proc == MPI_PROC_NULL ? true : false;
}