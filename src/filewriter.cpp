#include "filewriter.hpp"

FileWriter::FileWriter() {
  // ctor
}

FileWriter::~FileWriter() {
  // dtor
}

void FileWriter::WriteTec(const string &prefix, string sufix, int time,
                          const Grid &mainGrid, const Grid &myGrid,
                          const Variables &sol, int iStr, int iEnd, int jStr,
                          Loc loc) {
  int mainNI = mainGrid.NI;
  int NJ = myGrid.NJ;

  string timename;

  std::ostringstream temp;
  temp << time;

  timename = temp.str();

  string newname = "Sol";

  string name = prefix + newname;
  name.append(sufix);

  string newfilename = name.append(timename);

  std::ofstream outfile;

  string myType = ".dat";

  string name2 = newfilename.append(myType);

  int NXtemp = mainGrid.NI;
  int NYtemp = myGrid.NJ;
  if (loc == Loc::down) {
    outfile.open(name2);

    outfile << "VARIABLES=\"X\", \"Y\", \"U\", \"V\", \"T\", \"F\", \"Z\" \n";
    outfile << "ZONE T=\"Down\""
            << " ,I=" << NXtemp << ", J=" << NYtemp << ", DATAPACKING=POINT\n";
  } else {
    outfile.open(name2, std::ios::app);
    outfile << "ZONE T=\"Up\""
            << " ,I=" << NXtemp << ", J=" << NYtemp << ", DATAPACKING=POINT\n";
  }

  int jMax = myGrid.NJ + jStr;
  for (int i = iStr - 1; i < iEnd + 1; i++)
    for (int j = jStr; j < jMax; j++) {
      outfile << mainGrid.XC[i] << " ";
      outfile << mainGrid.YC[j] << " ";
      double value = sol.U.value[id(i, j, mainNI, NJ + jStr)];
      outfile << value << " ";
      value = sol.V.value[id(i, j, mainNI, NJ + jStr)];
      outfile << value << " ";
      value = sol.T.value[id(i, j, mainNI, NJ + jStr)];
      outfile << value << " ";
      value = sol.F.value[id(i, j, mainNI, NJ + jStr)];
      outfile << value << " ";
      value = sol.Z.value[id(i, j, mainNI, NJ + jStr)];
      outfile << value << "\n";
    }

  outfile.close();
}

void FileWriter::WriteInter(const string &prefix, string sufix, int time,
                            const Grid &mainGrid, const Grid &myGrid,
                            const Variables &sol, int iStr, int iEnd, int jStr,
                            Loc loc) {
  // Needed for q file
  double mach = 1;
  double alpha = 1;
  double reyn = 1;
  double qtime = 0;

  int mainNI = NTOTAL + 2;

  string timename;

  if (time == 0) {
    if (loc == Loc::down) {
      string name2 = prefix + "Grid" + sufix + ".xyz";

      std::ofstream outfile;

      outfile.open(name2, std::ios::binary);

      write_grid_header(outfile);

      for (int j = jStr; j < NJ + jStr; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char *)(&mainGrid.XC[i]), sizeof(mainGrid.XC[i]));

      for (int j = jStr; j < NJ + jStr - 1; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char *)(&mainGrid.YC[j]), sizeof(mainGrid.YC[j]));
      for (int i = 0; i < mainNI; i++)
        outfile.write((char *)(&myGrid.YC[NJ - 1]), sizeof(myGrid.YC[NJ - 1]));

      outfile.close();
    } else if (loc == Loc::up) {
      string myGridType = ".xyz";

      string newname = "Grid";

      string name2 = prefix + newname;
      name2.append(sufix);
      name2.append(myGridType);

      std::ofstream outfile;
      outfile.open(name2, std::ios::app | std::ios::binary);

      int jMax = myGrid.NJ + jStr;

      for (int j = jStr; j < jMax; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char *)(&mainGrid.XC[i]), sizeof(mainGrid.XC[i]));

      for (int i = 0; i < mainNI; i++)
        outfile.write((char *)(&myGrid.YC[0]), sizeof(myGrid.YC[0]));
      for (int j = jStr + 1; j < jMax; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char *)(&mainGrid.YC[j]), sizeof(mainGrid.YC[j]));

      outfile.close();
    }
  }

  std::ostringstream temp;
  temp << time;

  timename = temp.str();

  string name = prefix + "Sol" + sufix + timename;

  std::ofstream outfile;

  string myType = ".q";

  string name2 = name + myType;

  int Blocks = 2;
  int NXtemp = mainGrid.NI;
  int NYtemp = myGrid.NJ;
  if (loc == Loc::down) {
    outfile.open(name2, std::ios::binary);

    outfile.write((char *)&Blocks, sizeof(Blocks));

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
  } else {
    outfile.open(name2, std::ios::app | std::ios::binary);
  }
  // In q file there are 4 variables:
  outfile.write((char *)&mach, sizeof(mach));
  outfile.write((char *)&alpha, sizeof(alpha));
  outfile.write((char *)&reyn, sizeof(reyn));
  outfile.write((char *)&qtime, sizeof(qtime));

  // 1) Density
  // 2) U
  // 3) V
  // 4) Energy
  double rho = 1;
  int jMax = myGrid.NJ + jStr;

  write_variable(outfile, vector<double>(NTOTAL * MINPUT, rho), iStr, iEnd,
                 jStr, jMax);
  write_variable(outfile, sol.U.value, iStr, iEnd, jStr, jMax);
  write_variable(outfile, sol.V.value, iStr, iEnd, jStr, jMax);
  write_variable(outfile, sol.T.value, iStr, iEnd, jStr, jMax);

  outfile.close();

  string myFileType = ".f";

  string name3 = name + myFileType;

  if (loc == Loc::down) {
    outfile.open(name3, std::ios::out | std::ios::binary);

    int NVar = 3;

    outfile.write((char *)&Blocks, sizeof(Blocks));

    NXtemp = mainGrid.NI;
    NYtemp = myGrid.NJ;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
    outfile.write((char *)&NVar, sizeof(NVar));

    NXtemp = mainGrid.NI;
    NYtemp = myGrid.NJ;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
    outfile.write((char *)&NVar, sizeof(NVar));
  } else
    outfile.open(name3, std::ios::app | std::ios::binary);

  write_variable(outfile, sol.T.value, iStr, iEnd, jStr, jMax);
  write_variable(outfile, sol.F.value, iStr, iEnd, jStr, jMax);
  write_variable(outfile, sol.Z.value, iStr, iEnd, jStr, jMax);

  outfile.close();
}

void FileWriter::write_grid_header(std::ofstream &file) {
  int mainNI = NTOTAL + 2;
  int Blocks = 2;

  file.write((char *)&Blocks, sizeof(Blocks));

  file.write((char *)&mainNI, sizeof(mainNI));
  file.write((char *)&NJ, sizeof(NJ));

  file.write((char *)&mainNI, sizeof(mainNI));
  file.write((char *)&NJ, sizeof(NJ));
}

void FileWriter::write_variable(std::ofstream &file, const vector<double> &vec,
                                int iStr, int iEnd, int jStr, int jEnd) {
  int mainNI = NTOTAL + 2;
  for (int j = jStr; j < jEnd; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      double value = vec[id(i, j, mainNI, NJ + jStr)];
      file.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
}