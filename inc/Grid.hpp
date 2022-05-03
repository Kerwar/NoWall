#ifndef GRID_H
#define GRID_H

#include "ForAllOperators.hpp"

using std::vector;

class Grid
{
public:
  Grid();
  Grid(const int &n,const int &m,const double &xmin,const double &xmax,const double &ymin,const double &ymax);
  virtual ~Grid();


  void setX(double* & X);
  void setY(double* & Y);

  void setXC(double* &X, double* &XC);
  void setYC(double* &Y, double* &YC);

  void setXF(double* &X, double* &XC, double* &XF);
  void setYF(double* &Y, double* &YC, double* &YF);

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
  double *X, *Y;

  // Values of the Center of the Volumes
  double *XC, *YC;

  // (Distance of the face to the node) / (Distance Between nodes)
  double *XF, *YF;

  // First Aprox of the distances
  double DX, DY;

  void SetIEx(double &ExMin, double &ExMax);
public:
  friend void swap(Grid &first, Grid &second)
  {
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
  
  Grid& operator=(Grid other)
  {
    swap(*this, other);

    return *this;
  };
};

#endif