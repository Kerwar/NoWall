#include "Field.h"
#include <iostream>
#include <cmath>

Field::Field() : value(0.0), density(1.0)
{
}

Field::Field(int &NI_, int &NJ_) : value(0.0), NI(NI_), NJ(NJ_), density(1.0) // IF EVER NON CONST DENSITY THIS HAS TO CHANGE
{
}

Field::~Field()
{
}

void Field::getGridInfoPassed(Field::vectorField &f, Grid &myGrid, double &viscX, double &viscY)
{
  NI = myGrid.NI;
  NJ = myGrid.NJ;

  forAll(f)
  {
    f[i][j].X = myGrid.X[i + j*NI];
    f[i][j].Y = myGrid.Y[i + j*NI];
    f[i][j].XC = myGrid.XC[i + j*NI];
    f[i][j].YC = myGrid.YC[i + j*NI];
    f[i][j].FXE = myGrid.XF[i + j*NI];
    f[i][j].FYN = myGrid.YF[i + j*NI];
    f[i][j].FXP = 1.0 - myGrid.XF[i + j*NI];
    f[i][j].FYP = 1.0 - myGrid.YF[i + j*NI];

    f[i][j].viscX = viscX;
    f[i][j].viscY = viscY;
  }

  forAllInternal(f)
  {
    f[i][j].DXPtoE = std::abs(myGrid.XC[i + 1 + j*NI] - myGrid.XC[i + j*NI]);
    f[i][j].DYPtoN = std::abs(myGrid.YC[i + (j+1)*NI] - myGrid.YC[i + j*NI]);

    f[i][j].Se = std::abs(myGrid.Y[i + j*NI] - myGrid.Y[i + (j-1)*NI]);
    f[i][j].Sn = std::abs(myGrid.X[i + j*NI] - myGrid.X[i - 1 + j*NI]);

    f[i][j].volume = f[i][j].Se * f[i][j].Sn;
  }
}

void Field::inletBoundaryCondition(Field::vectorField &vec, Direction side, double bvalue)
{
  if (side == west)
  {
    forWestBoundary(vec)
    {
      vec[i][j].value = bvalue;
    }
  }
  else if (side == east)
  {
    forEastBoundary(vec)
    {
      vec[i][j].value = bvalue;
    }
  }
  else if (side == south)
  {
    forSouthBoundary(vec)
    {
      vec[i][j].value = bvalue;
    }
  }
  else if (side == north)
  {
    forNorthBoundary(vec)
    {
      vec[i][j].value = bvalue;
    }
  }
}

void Field::laminarFlow(Field::vectorField &vec, double m, double yMin, double yMax)
{
  forAll(vec)
  {
    double distanceToBot = vec[i][j].YC - yMin;
    double distanceToTop = yMax - vec[i][j].YC;
    // std::cout << i << " " << j << " " << distanceToBot << " " << distanceToTop << std::endl;
    vec[i][j].value = 6.0 * m / std::abs(m) * distanceToBot * distanceToTop;
  }
}

void Field::InitializeT(Field::vectorField &vec, double &q,
                        double xHotSpot, double yHotSpot, double xMin, double xMax)
{
  forAll(vec)
  {
    double pos = vec[i][j].XC - xHotSpot;
    double simmetricpos = vec[i][j].XC - (xMin + xMax - xHotSpot);

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));

    vec[i][j].value = q * 1.5 * tanV;
    if (std::abs(pos) > std::abs(simmetricpos))
    {
      tanV = (atan(simmetricpos) - atan(-(xMin - xMax) / 2)) /
             (atan((xMin - xMax) / 2) - atan(-(xMin - xMax) / 2));
      vec[i][j].value = std::max(q, q * 1.5 * tanV);
    }
  }
}

void Field::InitializeF(Field::vectorField &vec, double xHotSpot, double xMin, double xMax)
{
  forAll(vec)
  {
    double pos = vec[i][j].XC - xHotSpot;

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));
    vec[i][j].value = std::max(1.0 - tanV, 0.0);
  }
}

void Field::InitializeZ(Field::vectorField &vec, double T0hs, double r0hs,
                        double xHotSpot, double yHotSpot, double xMax)
{
  forAll(vec)
  {
    double rr = sqrt(pow(vec[i][j].XC - xHotSpot, 2) +
                     pow(vec[i][j].YC - yHotSpot, 2));

    vec[i][j].value = T0hs * exp(-rr / r0hs);
  }
}

void Field::initializeInternalField(Field::vectorField &vec, double val)
{
  forAllInternal(vec)
  {
    vec[i][j].value = val;
  }
}

void Field::linearExtrapolateCondition(Field::vectorField &vec, Field::Direction wallname)
{
  if (wallname == west)
  {
    forWestBoundary(vec)
    {
      vec[i][j].value = vec[i + 1][j].value + (vec[i + 1][j].value - vec[i + 2][j].value) * vec[i + 1][j].FXE;
    }
  }
  else if (wallname == east)
  {
    forEastBoundary(vec)
    {
      vec[i][j].value = vec[i - 1][j].value + (vec[i - 1][j].value - vec[i - 2][j].value) * vec[i - 1][j].FXP;
    }
  }
  else if (wallname == south)
  {
    forSouthBoundary(vec)
    {
      vec[i][j].value = vec[i][j + 1].value + (vec[i][j + 1].value - vec[i][j + 2].value) * vec[i][j + 1].FYN;
    }
  }
  else if (wallname == north)
  {
    forNorthBoundary(vec)
    {
      vec[i][j].value = vec[i][j - 1].value + (vec[i][j - 1].value - vec[i][j - 2].value) * vec[i][j - 1].FYP;
    }
  }
}

Field::vectorField Field::interpolatedFieldEast(Field::vectorField &vec, Grid &myGrid)
{
  Field::vectorField temp(vec.size(), vector<Field>(vec[0].size()));

  forAllInternalUCVs(vec)
  {
    double FXE = myGrid.XF[i + j*NI];
    double FXP = 1.0 - FXE;

    temp[i][j].value = vec[i + 1][j].value * FXE + vec[i][j].value * FXP;
  }

  return temp;
}

Field::vectorField Field::interpolatedFieldNorth(Field::vectorField &vec, Grid &myGrid)
{
  Field::vectorField temp(vec.size(), vector<Field>(vec[0].size()));

  forAllInternalVCVs(vec)
  {
    double FYN = myGrid.YF[i + j *NI];
    double FYP = 1.0 - FYN;

    temp[i][j].value = vec[i][j + 1].value * FYN + vec[i][j].value * FYP;
  }

  return temp;
}

void Field::computeEastMassFluxes(Field::vectorField &vec, Field::vectorField &corrU)
{
  // For non constant density forAllInternalUCVs
  forAllInternal(vec)
  {
    double sArea = vec[i][j].Se;
    double density = vec[i][j].density;

    vec[i][j].value = sArea * density * corrU[i][j].value;
  }
}

void Field::computeNorthMassFluxes(Field::vectorField &vec, Field::vectorField &corrV)
{
  // For non constant density forAllInternalVCVs
  forAllInternal(vec)
  {
    double sArea = vec[i][j].Sn;
    double density = vec[i][j].density;
    vec[i][j].value = sArea * density * corrV[i][j].value;
  }
}