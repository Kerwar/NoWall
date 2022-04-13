#ifndef GRID_H
#define GRID_H

#include <vector>
#include "ForAllOperators.h"

using std::vector;

class Grid
{
public:
  Grid();
  Grid(int &n, int &m, double &xmin, double &xmax, double &ymin, double &ymax);
  virtual ~Grid();

  //void setGrid(int &n, int &m, double &xmin, double &xmax, double &ymin, double &ymax);
  void setX(vector<vector<double>> &);
  void setY(vector<vector<double>> &);

  void setXC(vector<vector<double>> &, vector<vector<double>> &);
  void setYC(vector<vector<double>> &, vector<vector<double>> &);

  void setXF(vector<vector<double>> &, vector<vector<double>> &, vector<vector<double>> &);
  void setYF(vector<vector<double>> &, vector<vector<double>> &, vector<vector<double>> &);

  // Number of points introduced by the User
  int N, M;

  // Number of points that are needed for FVM
  int NI, NJ;
  int exI1, exI2;

  // Size of the chanels
  // (THIS IS FALSE RIGHT NOW) DISCLAIMER=> THIS GRID GOES FROM -XLENGTH TO XLENGTH AND FROM 0 TO YLENGTH
  double xMin, xMax;
  double yMin, yMax;
  double Xlength, Ylength;

  //  ___________ ___________ ___________ ___________
  // |           |           |           |           |
  // |     o     |     o     |     o     |     o     |
  // |___________|___________|___________|___________|
  // |           |           |           |           |
  // |     o     |     o     |     o     |     o     |
  // |___________|___________|___________|___________|
  // |           |           |           |           |
  // |     o     |     o     |     o     |     o     |
  // |___________|___________|___________|___________|
  // |           |           |           |           |
  // |     o     |     o     |     o     |     o     |
  // |___________|___________|___________|___________|
  // |           |           |           |           |
  // |     o     |     o     |     o     |     o     |
  // |___________|___________|___________|___________|
  // x1          x2          x3          x4        x5&x6
  // xc1  xc2         xc3         xc4         xc5   xc6

  // Values of the X, Y coordinates.
  vector<vector<double>> X, Y;

  // Values of the Center of the Volumes
  vector<vector<double>> XC, YC;

  // (Distance of the face to the node) / (Distance Between nodes)
  vector<vector<double>> XF, YF;

  // First Aprox of the distances
  double DX, DY;

  void SetIEx(double &ExMin, double &ExMax);
};

#endif