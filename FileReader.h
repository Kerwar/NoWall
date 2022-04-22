#ifndef FILEREADER_H
#define FILEREADER_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <cmath>

#include "Grid.h"
#include "Field.h"

using std::cout, std::endl;
using std::vector;
using std::string;

class FileReader
{
  public:

  FileReader();
  ~FileReader();
  void readField(string &name, int blockWanted, int variableWanted, Field &vec, int locIStr, int locIEnd);
};

#endif