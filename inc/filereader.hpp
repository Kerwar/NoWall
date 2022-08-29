#ifndef FILEREADER_H
#define FILEREADER_H

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "field.hpp"
#include "grid.hpp"

using std::cout, std::endl;
using std::string;
using std::vector;

class FileReader {
 public:
  FileReader();
  ~FileReader();
  void read_field(const string &name, int blockWanted, int variableWanted, Field &vec,
                 int locIStr, int locIEnd);
  void read_grid(const string &name, vector<double> &XC, vector<double> YC);                 

  private:

  void read_grid_header(std::ifstream &file, vector<int> &NX, vector<int> &NY);
  vector<double> read_variable(std::ifstream &file, int NX, int NY);
};

#endif