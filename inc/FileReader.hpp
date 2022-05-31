#ifndef FILEREADER_H
#define FILEREADER_H

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Field.hpp"
#include "Grid.hpp"

using std::cout, std::endl;
using std::string;
using std::vector;

class FileReader {
 public:
  FileReader();
  ~FileReader();
  void readField(string &name, int blockWanted, int variableWanted, Field &vec,
                 int locIStr, int locIEnd);
};

#endif