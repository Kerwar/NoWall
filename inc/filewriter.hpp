#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

#include "field.hpp"
#include "grid.hpp"
#include "paralel.hpp"
#include "systemofequations.hpp"

using std::string;

class FileWriter {
 public:
  FileWriter();
  virtual ~FileWriter();

  void WriteInter(const string &prefix, string sufix, int time,
                  const Grid &mainGrid, const Grid &myGrid,
                  const Variables &sol, int iStr, int iEnd, int jStr, Loc loc);

  void WriteTec(const string &prefix, string sufix, int time,
                const Grid &mainGrid, const Grid &myGrid, const Variables &sol,
                int iStr, int iEnd, int jStr, Loc loc);

  inline int id(const int &I, const int &J, const int &NI, const int &NJ) {
    return std::min(std::max(I, 0), NI - 3) +
           std::min(std::max(J, 0), NJ - 3) * (NI - 2);
  };

 private:
  string prd(const double x, const int decDigits, const int width);
  void write_grid_header(std::ofstream &file);
  void write_variable(std::ofstream &file, const vector<double> &vec, int iStr,
                      int iEnd, int jStr, int jEnd);
};

#endif  // FILEWRITER_H
