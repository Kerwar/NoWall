#include "Grid.h"
#include <iostream>

Grid::Grid()
{}

Grid::Grid(int &n, int &m, double &xmin, double &xmax, double &ymin, double &ymax) : N(n), M(m), xMin(xmin), xMax(xmax), yMin(ymin), yMax(ymax), Xlength(xMax - xMin), Ylength(yMax - yMin),
                                                                                     X(n + 2, vector<double>(m + 2, xmin)), Y(n + 2, vector<double>(m + 2, ymin)),
                                                                                     XC(n + 2, vector<double>(m + 2, xmin)), YC(n + 2, vector<double>(m + 2, ymin)),
                                                                                     XF(n + 2, vector<double>(m + 2, 0.0)), YF(n + 2, vector<double>(m + 2, 0.0))
{
  NI = N + 2;
  NJ = M + 2;
  // NIM = NI - 1; // REMOVABLE?
  // NJM = NJ - 1; // REMOVABLE?

  DX = Xlength / N;
  DY = Ylength / M;

  setX(X);
  setY(Y);
  setXC(X, XC);
  setYC(Y, YC);
  setXF(X, XC, XF);
  setYF(Y, YC, YF);
}

/*void Grid::setGrid(int &n, int &m, double &xmin, double &xmax, double &ymin, double &ymax):                X(n + 2, vector<double>(m + 2, xmin)), Y(n + 2, vector<double>(m + 2, ymin))
{
  N = n;
  M = m;
  xMin =xmin;
  xMax = xmax;
  yMin = ymin;
  yMax = ymax;
  Xlength  = xMax - xMin;
  Ylength = yMax - yMin;                                                                   
  
  vector<vector<double>> x(N + 2, vector<double>(M + 2, xMin));
  vector<vector<double>> y(N + 2, vector<double>(M + 2, xMin));
  vector<vector<double>> xc(N + 2, vector<double>(M + 2, xMin));
  vector<vector<double>> yc(N + 2, vector<double>(M + 2, yMin));
  vector<vector<double>> xf(N + 2, vector<double>(M + 2, 0.0));
  vector<vector<double>> yf(N + 2, vector<double>(M + 2, 0.0));
  XC = xc;
  YC = yc;
  XF(n + 2, vector<double>(m + 2, 0.0));
  YF(n + 2, vector<double>(m + 2, 0.0));

  NI = N + 2;
  NJ = M + 2;
  // NIM = NI - 1; // REMOVABLE?
  // NJM = NJ - 1; // REMOVABLE?

  DX = Xlength / N;
  DY = Ylength / M;

  setX(X);
  setY(Y);
  setXC(X, XC);
  setYC(Y, YC);
  setXF(X, XC, XF);
  setYF(Y, YC, YF);
}
*/
void Grid::setX(vector<vector<double>> &vecX)
{
  for (unsigned int i = 1; i < vecX.size() - 1; i++)
  {
    for (unsigned int j = 0; j < vecX[0].size(); j++)
    {
      vecX[i][j] = vecX[i - 1][j] + DX;
    }
  }

  for (unsigned int j = 0; j < vecX[0].size(); j++)
  {
    vecX[vecX.size() - 1][j] = vecX[vecX.size() - 2][j];
  }
}

void Grid::setY(vector<vector<double>> &vecY)
{
  for (unsigned int i = 0; i < vecY.size(); i++)
  {
    for (unsigned int j = 1; j < vecY[0].size() - 1; j++)
    {
      vecY[i][j] = vecY[i][j - 1] + DY;
    }
  }

  for (unsigned int i = 0; i < vecY.size(); i++)
  {
    vecY[i][vecY[i].size() - 1] = vecY[i][vecY[i].size() - 2];
  }
}

void Grid::setXC(vector<vector<double>> &vecX, vector<vector<double>> &vecXC)
{
  for (unsigned int i = 1; i < vecX.size(); i++)
  {
    for (unsigned int j = 0; j < vecX[0].size(); j++)
    {
      vecXC[i][j] = (vecX[i][j] + vecX[i - 1][j]) * 0.5;
    }
  }
}

void Grid::setYC(vector<vector<double>> &vecY, vector<vector<double>> &vecYC)
{
  for (unsigned int i = 0; i < vecY.size(); i++)
  {
    for (unsigned int j = 1; j < vecY[0].size(); j++)
    {
      vecYC[i][j] = (vecY[i][j] + vecY[i][j - 1]) * 0.5;
    }
  }
}

void Grid::setXF(vector<vector<double>> &vecX, vector<vector<double>> &vecXC, vector<vector<double>> &vecXF)
{
  for (unsigned int i = 0; i < vecX.size() - 1; i++)
  {
    for (unsigned int j = 0; j < vecX[0].size(); j++)
    {
      vecXF[i][j] = (vecX[i][j] - vecXC[i][j]) / (vecXC[i + 1][j] - vecXC[i][j]);
    }
  }

  // for(unsigned int j = 0; j < vecX[0].size(); j++)
  //   {
  //     vecXF[vecX.size()-1][j] = 0.0;
  //   }
}

void Grid::setYF(vector<vector<double>> &vecY, vector<vector<double>> &vecYC, vector<vector<double>> &vecYF)
{
  for (unsigned int i = 0; i < vecY.size(); i++)
  {
    for (unsigned int j = 0; j < vecY[0].size() - 1; j++)
    {
      vecYF[i][j] = (vecY[i][j] - vecYC[i][j]) / (vecYC[i][j + 1] - vecYC[i][j]);
    }
  }
  //  for(unsigned int i = 0; i < vecY.size(); i++)
  //   {
  //       vecYF[i][vecY[i].size()-1] = 0.0;
  //   }
}

void Grid::SetIEx(double &ExMin, double &ExMax)
{
  int NXWall = 0;

  while (NXWall < NI - 1 && XC[NXWall][0] < ExMin)
  {
    // std::cout << NXWall << " " << XC[NXWall][0] << std::endl;
    NXWall++;
  }
  exI1 = NXWall >= 0 && NXWall < NI - 1 ? NXWall : 0;

  while (NXWall < NI - 1 && XC[NXWall][0] < ExMax)
  {
    NXWall++;
  }
  if (NXWall == 0)
  {
    exI2 = exI1;
  }
  else
  {
    exI2 = XC[NXWall - 1][0] < ExMax && XC[exI1][0] >= ExMin ? NXWall: exI1;
  }
};

Grid::~Grid()
{
}