#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <string>
#include <cstring>
#include <sstream>

#include "Grid.hpp"
#include "Field.hpp"
#include "Paralel.hpp"

using std::string;

class FileWriter
{
public:
  FileWriter();
  virtual ~FileWriter();

  // void writeTFZ(string &name, int time, Grid &myGrid, Field::vectorField &Ttemp, Field::vectorField &Ftemp, Field::vectorField &Ztemp, int iStr, int iEnd, int jStr, int jEnd, Paralel::Loc loc);
  void WriteBin(string prefix, string sufix, int time, Grid &myGrid, Field &Utemp, Field &Vtemp, Field &Ttemp, Field &Ftemp, Field &Ztemp, int iStr, int iEnd, int jStr, int jEnd, Paralel::Loc loc);
  void WriteInter(string prefix, string sufix, int time,
                  const Grid &mainGrid, const Grid &myGrid,
                  const Field &Utemp, const Field &Vtemp,
                  const Field &Ttemp, const Field &Ftemp, const Field &Ztemp,
                  int iStr, int iEnd, int jStr, Paralel::Loc loc);

  const inline int id(const int &I, const int &J, const int &NI, const int &NJ) { return std::min(std::max(I, 0), NI-3) + std::min(std::max(J, 0), NJ-3) * (NI-2); };

private:
  string prd(const double x, const int decDigits, const int width);
};

#endif // FILEWRITER_H
