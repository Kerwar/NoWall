#include "FileWriter.hpp"

FileWriter::FileWriter() {
  // ctor
}

FileWriter::~FileWriter() {
  // dtor
}

void FileWriter::WriteInter(const string &prefix, string sufix, int time,
                            const Grid &mainGrid, const Grid &myGrid,
                            const Field &Utemp, const Field &Vtemp,
                            const Field &Ttemp, const Field &Ftemp,
                            const Field &Ztemp, int iStr, int iEnd, int jStr,
                            Loc loc) {
  // Needed for q file
  double mach = 1;
  double alpha = 1;
  double reyn = 1;
  double qtime = 0;
  int ttime = time;

  int mainNI = mainGrid.NI;
  int NJ = myGrid.NJ;

  string timename;

  if (time == 0) {
    if (loc == Loc::down) {
      string myGridType = ".xyz";

      string newname = "Grid";

      string name2 = prefix + newname;

      name2.append(sufix);
      name2.append(myGridType);

      char cstr[name2.size() + 2];
      strcpy(cstr, name2.c_str());

      std::ofstream outfile;

      outfile.open(cstr, std::ios::binary);

      int Blocks = 2;
      outfile.write((char *)&Blocks, sizeof(Blocks));

      int NXtemp = mainNI;
      int NYtemp = NJ;

      outfile.write((char *)&NXtemp, sizeof(NXtemp));
      outfile.write((char *)&NYtemp, sizeof(NYtemp));

      NXtemp = mainNI;
      NYtemp = NJ;

      outfile.write((char *)&NXtemp, sizeof(NXtemp));
      outfile.write((char *)&NYtemp, sizeof(NYtemp));

      for (int j = jStr; j < NJ + jStr; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char*)(&mainGrid.XC[i]),
                        sizeof(mainGrid.XC[i]));

      for (int j = jStr; j < NJ + jStr - 1; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char*)(&mainGrid.YC[j]),
                        sizeof(mainGrid.YC[j]));
      for (int i = 0; i < mainNI; i++)
        outfile.write((char*)(&myGrid.YC[NJ - 1]),
                      sizeof(myGrid.YC[NJ - 1]));

      outfile.close();
    } else if (loc == Loc::up) {
      string myGridType = ".xyz";

      string newname = "Grid";

      string name2 = prefix + newname;
      name2.append(sufix);
      name2.append(myGridType);

      char cstr[name2.size() + 2];
      strcpy(cstr, name2.c_str());

      std::ofstream outfile;
      outfile.open(cstr, std::ios::app | std::ios::binary);

      int jMax = myGrid.NJ + jStr;

      for (int j = jStr; j < jMax; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char*)(&mainGrid.XC[i]),
                        sizeof(mainGrid.XC[i]));

      for (int i = 0; i < mainNI; i++)
        outfile.write((char*)(&myGrid.YC[0]),
                      sizeof(myGrid.YC[0]));
      for (int j = jStr + 1; j < jMax; j++)
        for (int i = 0; i < mainNI; i++)
          outfile.write((char*)(&mainGrid.YC[j]),
                        sizeof(mainGrid.YC[j]));

      outfile.close();
    }
  }

  std::ostringstream temp;
  temp << ttime;

  timename = temp.str();

  string newname = "Sol";

  string name = prefix + newname;
  name.append(sufix);

  string newfilename = name.append(timename);

  std::ofstream outfile;

  string myType = ".q";

  string name2 = newfilename.append(myType);
  char qstr[name2.size() + 2];
  strcpy(qstr, name2.c_str());

  int Blocks = 2;
  int NXtemp = mainGrid.NI;
  int NYtemp = myGrid.NJ;
  if (loc == Loc::down) {
    outfile.open(qstr, std::ios::binary);

    outfile.write((char *)&Blocks, sizeof(Blocks));

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
  } else {
    outfile.open(qstr, std::ios::app | std::ios::binary);
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
  for (int j = jStr; j < jMax; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      outfile.write(reinterpret_cast<char *>(&rho), sizeof(rho));
    }
  for (int j = jStr; j < jMax; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      double value = Utemp.value[id(i, j, mainNI, NJ + jStr)];
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  for (int j = jStr; j < jMax; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      double value = Vtemp.value[id(i, j, mainNI, NJ + jStr)];
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }

  for (int j = jStr; j < jMax; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      double value = Ttemp.value[id(i, j, mainNI, NJ + jStr)];
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  outfile.close();

  newname = "Sol";

  name = prefix + newname;
  name.append(sufix);

  newfilename = name.append(timename);

  string myFileType = ".f";

  string name3 = newfilename.append(myFileType);
  char cstr[name3.size() + 2];
  strcpy(cstr, newfilename.c_str());
  if (loc == Loc::down) {
    outfile.open(cstr, std::ios::out | std::ios::binary);

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
    outfile.open(cstr, std::ios::app | std::ios::binary);

  for (int j = jStr; j < jMax; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      double value = Ttemp.value[id(i, j, mainNI, NJ + jStr)];
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  for (int j = jStr; j < jMax; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      double value = Ftemp.value[id(i, j, mainNI, NJ + jStr)];
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  for (int j = jStr; j < jMax; j++)
    for (int i = iStr - 1; i < iEnd + 1; i++) {
      double value = Ztemp.value[id(i, j, mainNI, NJ + jStr)];
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  outfile.close();
}
