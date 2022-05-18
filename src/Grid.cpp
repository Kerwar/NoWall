#include "Grid.hpp"
#include <iostream>

Grid::Grid() : X(NULL), Y(NULL), XC(NULL), YC(NULL), XF(NULL), YF(NULL)
{
}

Grid::Grid(const int &n, const int &m,
           const double &xmin, const double &xmax,
           const double &ymin, const double &ymax) : N(n), M(m), NI(n + 2), NJ(m + 2),
                                                     xMin(xmin), xMax(xmax), yMin(ymin), yMax(ymax), Xlength(xMax - xMin), Ylength(yMax - yMin),
                                                     X(new double[NI]), Y(new double[NJ]), XC(new double[NI]), YC(new double[NJ]),
                                                     XF(new double[NI]), YF(new double[NJ])
{
  DX = Xlength / N;
  DY = Ylength / M;

  setX(X);
  setY(Y);
  setXC(X, XC);
  setYC(Y, YC);
  setXF(X, XC, XF);
  setYF(Y, YC, YF);
}

void Grid::setX(std::shared_ptr<double[]> &vecX)
{
  vecX[0] = xMin;

  for (int i = 1; i < NI - 1; i++)
    vecX[i] = vecX[i - 1] + DX;

  vecX[NI - 1] = vecX[NI - 2];
}

void Grid::setY(std::shared_ptr<double[]> &vecY)
{
  vecY[0] = yMin;

  for (int j = 1; j < NJ - 1; j++)
    vecY[j] = vecY[j - 1] + DY;

  vecY[NJ - 1] = vecY[NJ - 2];
}

void Grid::setXC(std::shared_ptr<double[]> &vecX, std::shared_ptr<double[]> &vecXC)
{
  vecXC[0] = xMin;

  for (int i = 1; i < NI; i++)
    vecXC[i] = (vecX[i] + vecX[i - 1]) * 0.5;
}

void Grid::setYC(std::shared_ptr<double[]> &vecY, std::shared_ptr<double[]> &vecYC)
{
  vecYC[0] = yMin;

  for (int j = 1; j < NJ; j++)
    vecYC[j] = (vecY[j] + vecY[j - 1]) * 0.5;
}

void Grid::setXF(std::shared_ptr<double[]> &vecX, std::shared_ptr<double[]> &vecXC, std::shared_ptr<double[]> &vecXF)
{
  for (int i = 0; i < NI - 1; i++)
    vecXF[i] = (vecX[i] - vecXC[i]) / (vecXC[i + 1] - vecXC[i]);

  vecXF[0] = 0.0;
  vecXF[NI - 1] = 0.0;
}

void Grid::setYF(std::shared_ptr<double[]> &vecY, std::shared_ptr<double[]> &vecYC, std::shared_ptr<double[]> &vecYF)
{
  for (int j = 0; j < NJ - 1; j++)
    vecYF[j] = (vecY[j] - vecYC[j]) / (vecYC[j + 1] - vecYC[j]);

  vecYF[0] = 0.0;
  vecYF[NJ - 1] = 0.0;
}

void Grid::SetIEx(const double &ExMin, const double &ExMax)
{
  int NXWall = 0;

  while (NXWall < NI - 1 && XC[NXWall] < ExMin)
    NXWall++;

  exI1 = NXWall >= 0 && NXWall < NI - 1 ? NXWall : 0;

  while (NXWall < NI - 1 && XC[NXWall] < ExMax)
    NXWall++;

  if (NXWall == 0)
    exI2 = exI1;
  else
    exI2 = XC[NXWall - 1] < ExMax && XC[exI1] >= ExMin ? NXWall : exI1;
}

Grid::~Grid()
{
}
