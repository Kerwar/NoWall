#ifndef GRID_H
#define GRID_H

#include <memory>

#include "ForAllOperators.hpp"

using std::vector;

class Grid {
 public:
  Grid();
  Grid(const int &n, const int &m, const double &xmin, const double &xmax,
       const double &ymin, const double &ymax);
  virtual ~Grid();

  void setX(vector<double> &X);
  void setY(vector<double> &Y);

  void setXC(vector<double> &X, vector<double> &XC);
  void setYC(vector<double> &Y, vector<double> &YC);

  void setXF(vector<double> &X, vector<double> &XC,
             vector<double> &XF);
  void setYF(vector<double> &Y, vector<double> &YC,
             vector<double> &YF);

  // Number of points introduced by the User
  int N, M;

  // Number of points that are needed for FVM
  int NI, NJ;
  int exI1, exI2;

  // Size of the chanels
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
  vector<double> X, Y;

  // Values of the Center of the Volumes
  vector<double> XC, YC;

  // (Distance of the face to the node) / (Distance Between nodes)
  vector<double> XF, YF;

  // First Aprox of the distances
  double DX, DY;

  void SetIEx(const double &ExMin, const double &ExMax);

 public:
  friend void swap(Grid &first, Grid &second) {
    using std::swap;

    swap(first.N, second.N);
    swap(first.M, second.M);

    swap(first.NI, second.NI);
    swap(first.NJ, second.NJ);

    swap(first.exI1, second.exI1);
    swap(first.exI2, second.exI2);

    swap(first.xMin, second.xMin);
    swap(first.xMax, second.xMax);

    swap(first.yMin, second.yMin);
    swap(first.yMax, second.yMax);

    swap(first.Xlength, second.Xlength);
    swap(first.Ylength, second.Ylength);

    swap(first.X, second.X);
    swap(first.Y, second.Y);
    swap(first.XC, second.XC);
    swap(first.YC, second.YC);
    swap(first.XF, second.XF);
    swap(first.YF, second.YF);

    swap(first.DX, second.DX);
    swap(first.DY, second.DY);
  };

  Grid &operator=(Grid other) {
    swap(*this, other);

    return *this;
  };
};

#endif