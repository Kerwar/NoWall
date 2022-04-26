#include "Grid.hpp"
#include <iostream>

Grid::Grid() : X(NULL), Y(NULL), XC(NULL), YC(NULL), XF(NULL), YF(NULL)
{
}

Grid::Grid(int &n, int &m, double &xmin, double &xmax, double &ymin, double &ymax) : N(n), M(m), NI(n + 2), NJ(m + 2),
                                                                                     xMin(xmin), xMax(xmax), yMin(ymin), yMax(ymax), Xlength(xMax - xMin), Ylength(yMax - yMin),
                                                                                     X(new double[NI * NJ]), Y(new double[NI * NJ]), XC(new double[NI * NJ]), YC(new double[NI * NJ]),
                                                                                     XF(new double[NI * NJ]), YF(new double[NI * NJ])
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

void Grid::setX(double *&vecX)
{
  for (int j = 0; j < NJ; j++)
    vecX[j * NI] = xMin;

  for (int i = 1; i < NI - 1; i++)
    for (int j = 0; j < NJ; j++)
      vecX[i + j * NI] = vecX[i - 1 + j * NI] + DX;

  for (int j = 0; j < NJ; j++)
    vecX[NI - 1 + j * NI] = vecX[NI - 2 + j * NI];
}

void Grid::setY(double *&vecY)
{
  for (int i = 0; i < NI; i++)
    vecY[i] = yMin;

  for (int i = 0; i < NI; i++)
    for (int j = 1; j < NJ - 1; j++)
      vecY[i + j * NI] = vecY[i + (j - 1) * NI] + DY;

  for (int i = 0; i < NI; i++)
    vecY[i + (NJ - 1) * NI] = vecY[i + (NJ - 2) * NI];
}

void Grid::setXC(double *&vecX, double *&vecXC)
{
  for (int j = 0; j < NJ; j++)
    vecXC[j * NI] = xMin;

  for (int i = 1; i < NI; i++)
    for (int j = 0; j < NJ; j++)
      vecXC[i + j * NI] = (vecX[i + j * NI] + vecX[i - 1 + j * NI]) * 0.5;
}

void Grid::setYC(double *&vecY, double *&vecYC)
{
  for (int i = 0; i < NI; i++)
    vecYC[i] = yMin;

  for (int i = 0; i < NI; i++)
    for (int j = 1; j < NJ; j++)
      vecYC[i + j * NI] = (vecY[i + j * NI] + vecY[i + (j - 1) * NI]) * 0.5;
}

void Grid::setXF(double *&vecX, double *&vecXC, double *&vecXF)
{
  for (int i = 0; i < NI - 1; i++)
    for (int j = 0; j < NJ; j++)
      vecXF[i + j * NI] = (vecX[i + j * NI] - vecXC[i + j * NI]) / (vecXC[i + 1 + j * NI] - vecXC[i + j * NI]);

  for (int j = 0; j < NJ; j++)
  {
    vecXF[j * NI] = 0.0;
    vecXF[NI - 1 + j * NI] = 0.0;
  }
}

void Grid::setYF(double *&vecY, double *&vecYC, double *&vecYF)
{
  for (int i = 0; i < NI; i++)
    for (int j = 0; j < NJ - 1; j++)
      vecYF[i + j * NI] = (vecY[i + j * NI] - vecYC[i + j * NI]) / (vecYC[i + (j + 1) * NI] - vecYC[i + j * NI]);

  for (int i = 0; i < NI; i++)
  {
    vecYF[i] = 0.0;
    vecYF[i + (NJ - 1) * NI] = 0.0;
  }
}

void Grid::SetIEx(double &ExMin, double &ExMax)
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
  delete[] X;
  delete[] Y;
  delete[] XC;
  delete[] YC;
  delete[] XF;
  delete[] YF;
}
