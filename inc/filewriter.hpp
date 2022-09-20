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
                  const Variables &sol, Loc loc);

  int id(const int &I, const int &J, const int &NX, const int &jStr, const int &jEnd);

 private:
  string prd(const double x, const int decDigits, const int width);
  void write_grid_header(std::ofstream &file);
  void write_variable(std::ofstream &file, const vector<double> &vec, int jStr, int jEnd);
  void write_grid_x(std::ofstream &file, const vector<double> &grid);
  void write_grid_y(std::ofstream &file, const vector<double> &grid);
};

#endif  // FILEWRITER_H
