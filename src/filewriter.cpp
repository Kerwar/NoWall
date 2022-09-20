#include "filewriter.hpp"

FileWriter::FileWriter() {
  // ctor
}

FileWriter::~FileWriter() {
  // dtor
}

void FileWriter::WriteInter(const string &prefix, string sufix, int time,
                            const Grid &mainGrid, const Grid &myGrid,
                            const Variables &sol, Loc loc) {
  // Needed for q file
  double mach = 1;
  double alpha = 1;
  double reyn = 1;
  double qtime = 0;

  int mainNI = NTOTAL + 2;

  string timename;

  if (time == 0) {
    std::string name2 =
        prefix + "Grid_" + to_string(NTOTAL) + "x" + to_string(MINPUT) + ".xyz";
    std::ofstream outfile;
    if (loc == Loc::down) {
      outfile.open(name2, std::ios::binary);

      write_grid_header(outfile);

      write_grid_x(outfile, mainGrid.XC);
      write_grid_y(outfile, myGrid.YC);

    } else if (loc == Loc::up) {
      outfile.open(name2, std::ios::app | std::ios::binary);

      write_grid_x(outfile, mainGrid.XC);
      write_grid_y(outfile, myGrid.YC);
    }
    outfile.close();
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
  int jStr = loc == Loc::down ? 0 : myGrid.NJ - 2;
  int jMax = myGrid.M + jStr;

  write_variable(outfile, vector<double>(NTOTAL * MINPUT, rho), jStr, jMax);
  write_variable(outfile, sol.U.value, jStr, jMax);
  write_variable(outfile, sol.V.value, jStr, jMax);
  write_variable(outfile, sol.T.value, jStr, jMax);

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

  write_variable(outfile, sol.T.value, jStr, jMax);
  write_variable(outfile, sol.F.value, jStr, jMax);
  write_variable(outfile, sol.Z.value, jStr, jMax);

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
                                int jStr, int jEnd) {
  int mainNI = NTOTAL + 2;
  for (int j = jStr; j < jStr + NJ; j++)
    for (int i = 0; i < mainNI; i++) {
      double value = vec[id(i, j, mainNI, jStr, jEnd)];
      file.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
}

void FileWriter::write_grid_x(std::ofstream &file, const vector<double> &grid) {
  int mainNI = NTOTAL + 2;
  for (int j = 0; j < NJ; j++)
    for (int i = 0; i < mainNI; i++)
      file.write((char *)(&grid[i]), sizeof(grid[i]));
}

void FileWriter::write_grid_y(std::ofstream &file, const vector<double> &grid) {
  int mainNI = NTOTAL + 2;
  for (int j = 0; j < NJ; j++)
    for (int i = 0; i < mainNI; i++)
      file.write((char *)(&grid[j]), sizeof(grid[j]));
}

int FileWriter::id(const int &I, const int &J, const int &NX, const int &jStr,
                   const int &jEnd) {
  return std::min(std::max(I - 1, 0), NTOTAL - 1) +
         std::min(std::max(J - 1, jStr), jEnd - 1) * NTOTAL;
};