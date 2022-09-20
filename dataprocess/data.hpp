#ifndef _DATA_HPP_
#define _DATA_HPP_

#include<vector>
using std::vector;
using std::string;

struct Data {
  Data(){};
  double a, q, m;
  double T_before1, T_after1, T;
  double Z_Max_prev, Z_Max, Z_Max_post, Z;
  int indexes;
  int N, M;

  string dims;
};

#endif