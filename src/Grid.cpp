#include "Grid.hpp"

#include <iostream>

Grid::Grid() {}

Grid::Grid(const int &n, const int &m, const double &xmin, const double &xmax,
           const double &ymin, const double &ymax)
    : N(n),
      M(m),
      NI(n + 2),
      NJ(m + 2),
      xMin(xmin),
      xMax(xmax),
      yMin(ymin),
      yMax(ymax),
      Xlength(xMax - xMin),
      Ylength(yMax - yMin),
      X(vector<double>(NI,0)),
      Y(vector<double>(NJ,0)),
      XC(vector<double>(NI,0)),
      YC(vector<double>(NJ,0)),
      XF(vector<double>(NI,0)),
      YF(vector<double>(NJ,0)) {
  DX = Xlength / N;
  DY = Ylength / M;

  setX(X);
  setY(Y);
  setXC(X, XC);
  setYC(Y, YC);
  setXF(X, XC, XF);
  setYF(Y, YC, YF);
}

void Grid::setX(vector<double> &vecX) {
  vecX[0] = xMin;

  for (int i = 1; i < NI - 1; i++) vecX[i] = vecX[i - 1] + DX;

  vecX[NI - 1] = vecX[NI - 2];
}

void Grid::setY(vector<double> &vecY) {
  vecY[0] = yMin;

  for (int j = 1; j < NJ - 1; j++) vecY[j] = vecY[j - 1] + DY;

  vecY[NJ - 1] = vecY[NJ - 2];
}

void Grid::setXC(vector<double> &vecX,
                 vector<double> &vecXC) {
  vecXC[0] = xMin;

  for (int i = 1; i < NI; i++) vecXC[i] = (vecX[i] + vecX[i - 1]) * 0.5;
}

void Grid::setYC(vector<double> &vecY,
                 vector<double> &vecYC) {
  vecYC[0] = yMin;

  for (int j = 1; j < NJ; j++) vecYC[j] = (vecY[j] + vecY[j - 1]) * 0.5;
}

void Grid::setXF(vector<double> &vecX,
                 vector<double> &vecXC,
                 vector<double> &vecXF) {
  for (int i = 0; i < NI - 1; i++)
    vecXF[i] = (vecX[i] - vecXC[i]) / (vecXC[i + 1] - vecXC[i]);

  vecXF[0] = 0.0;
  vecXF[NI - 1] = 0.0;
}

void Grid::setYF(vector<double> &vecY,
                 vector<double> &vecYC,
                 vector<double> &vecYF) {
  for (int j = 0; j < NJ - 1; j++)
    vecYF[j] = (vecY[j] - vecYC[j]) / (vecYC[j + 1] - vecYC[j]);

  vecYF[0] = 0.0;
  vecYF[NJ - 1] = 0.0;
}

void Grid::SetIEx(const double &ExMin, const double &ExMax) {
  int NXWall = 0;

  while (NXWall < NI - 1 && XC[NXWall] < ExMin) NXWall++;

  exI1 = NXWall >= 0 && NXWall < NI - 1 ? NXWall : 0;

  while (NXWall < NI - 1 && XC[NXWall] < ExMax) NXWall++;

  if (NXWall == 0)
    exI2 = exI1;
  else
    exI2 = XC[NXWall - 1] < ExMax && XC[exI1] >= ExMin ? NXWall : exI1;
}

Grid::~Grid() {}
