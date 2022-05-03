#include "Field.hpp"

// Field::Field(): value(NULL),
//   X(NULL), XC(NULL),
//   FXE(NULL), FXP(NULL),
//   Y(NULL), YC(NULL),
//   FYN(NULL), FYP(NULL),
//   DXPtoE(NULL), DYPtoN(NULL),
//   Se(NULL), Sn(NULL),
//   viscX(NULL), viscY(NULL),
//   density(NULL), volume(NULL)
// {
// }

Field::Field(const int &_NI, const int &_NJ) : NI(_NI), NJ(_NJ), value(new double[NI * NJ]),
                                               X(new double[NI * NJ]), XC(new double[NI * NJ]),
                                               FXE(new double[NI * NJ]), FXP(new double[NI * NJ]),
                                               Y(new double[NI * NJ]), YC(new double[NI * NJ]),
                                               FYN(new double[NI * NJ]), FYP(new double[NI * NJ]),
                                               DXPtoE(new double[NI * NJ]), DYPtoN(new double[NI * NJ]),
                                               Se(new double[NI * NJ]), Sn(new double[NI * NJ]),
                                               viscX(new double[NI * NJ]), viscY(new double[NI * NJ]),
                                               density(new double[NI * NJ]), volume(new double[NI * NJ])
// IF EVER NON CONST DENSITY THIS HAS TO CHANGE
{
}

Field::~Field()
{
  delete[] value;
  delete[] X;
  delete[] XC;
  delete[] FXE;
  delete[] FXP;
  delete[] Y;
  delete[] YC;
  delete[] FYN;
  delete[] FYP;
  delete[] DXPtoE;
  delete[] DYPtoN;
  delete[] Se;
  delete[] Sn;
  delete[] viscX;
  delete[] viscY;
  delete[] density;
  delete[] volume;
}

void Field::getGridInfoPassed(const Grid &myGrid, double &viscX_, double &viscY_)
{
  forAllN(NI, NJ)
  {
    int index = i + j * NI;
    X[index] = myGrid.X[index];
    Y[index] = myGrid.Y[index];
    XC[index] = myGrid.XC[index];
    YC[index] = myGrid.YC[index];
    FXE[index] = myGrid.XF[index];
    FYN[index] = myGrid.YF[index];
    FXP[index] = 1.0 - myGrid.XF[index];
    FYP[index] = 1.0 - myGrid.YF[index];
    density[index] = 1.0;
    viscX[index] = viscX_;
    viscY[index] = viscY_;
  }

  forAllInterior(NI, NJ)
  {
    int index = i + j * NI;
    DXPtoE[index] = std::abs(myGrid.XC[index + 1] - myGrid.XC[index]);
    DYPtoN[index] = std::abs(myGrid.YC[index + NI] - myGrid.YC[index]);

    Se[index] = std::abs(myGrid.Y[index] - myGrid.Y[index - NI]);
    Sn[index] = std::abs(myGrid.X[index] - myGrid.X[index - 1]);

    volume[index] = Se[index] * Sn[index];
  }
}

void Field::inletBoundaryCondition(Direction side, double bvalue)
{
  if (side == west)
  {
    forWBoundary(NI, NJ)
    {
      value[i + j * NI] = bvalue;
    }
  }
  else if (side == east)
  {
    forEBoundary(NI, NJ)
    {
      value[i + j * NI] = bvalue;
    }
  }
  else if (side == south)
  {
    forSBoundary(NI, NJ)
    {
      value[i + j * NI] = bvalue;
    }
  }
  else if (side == north)
  {
    forNBoundary(NI, NJ)
    {
      value[i + j * NI] = bvalue;
    }
  }
}

void Field::laminarFlow(double m, double yMin, double yMax)
{
  forAllN(NI, NJ)
  {
    double distanceToBot = YC[i + j * NI] - yMin;
    double distanceToTop = yMax - YC[i + j * NI];
    // std::cout << i << " " << j << " " << distanceToBot << " " << distanceToTop << std::endl;
    value[i + j * NI] = 6.0 * m / std::abs(m) * distanceToBot * distanceToTop;
  }
}

void Field::InitializeT(double &q, double xHotSpot, double yHotSpot,
                        double xMin, double xMax)
{
  forAllN(NI, NJ)
  {
    double pos = XC[i + j * NI] - xHotSpot;
    double simmetricpos = XC[i + j * NI] - (xMin + xMax - xHotSpot);

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));

    value[i + j * NI] = q * 1.5 * tanV;
    if (std::abs(pos) > std::abs(simmetricpos))
    {
      tanV = (atan(simmetricpos) - atan(-(xMin - xMax) / 2)) /
             (atan((xMin - xMax) / 2) - atan(-(xMin - xMax) / 2));
      value[i + j * NI] = std::max(q, q * 1.5 * tanV);
    }
  }
}

void Field::InitializeF(double xHotSpot, double xMin, double xMax)
{
  forAllN(NI, NJ)
  {
    double pos = XC[i + j * NI] - xHotSpot;

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));
    value[i + j * NI] = std::max(1.0 - tanV, 0.0);
  }
}

void Field::InitializeZ(double T0hs, double r0hs,
                        double xHotSpot, double yHotSpot, double xMax)
{
  forAllN(NI, NJ)
  {
    double rr = sqrt(pow(XC[i + j * NI] - xHotSpot, 2) +
                     pow(YC[i + j * NI] - yHotSpot, 2));

    value[i + j * NI] = T0hs * exp(-rr / r0hs);
  }
}

void Field::initializeInternalField(double val)
{
  forAllInterior(NI, NJ)
  {
    value[i + j * NI] = val;
  }
}

void Field::linearExtrapolateCondition(const Field::Direction &wallname)
{
  if (wallname == west)
  {
    forWBoundary(NI, NJ)
    {
      int index = i + j * NI;
      value[index] = value[index + 1] + (value[index + 1] - value[index + 2]) * FXE[index + 1];
    }
  }
  else if (wallname == east)
  {
    forEBoundary(NI, NJ)
    {
      int index = i + j * NI;
      value[index] = value[index - 1] + (value[index - 1] - value[index - 2]) * FXP[index - 1];
    }
  }
  else if (wallname == south)
  {
    forSBoundary(NI, NJ)
    {
      int index = i + j * NI;
      value[index] = value[index + NI] + (value[index + NI] - value[index + 2 * NI]) * FYN[index + NI];
    }
  }
  else if (wallname == north)
  {
    forNBoundary(NI, NJ)
    {
      int index = i + j * NI;
      value[index] = value[index - NI] + (value[index - NI] - value[index - 2 * NI]) * FYP[index - NI];
    }
  }
}

void Field::interpolatedFieldEast(const Field &vec, const Grid &myGrid)
{
  forAllInteriorUCVs(NI, NJ)
  {
    double FXE = myGrid.XF[i + j * NI];
    double FXP = 1.0 - FXE;

    value[i + j * NI] = vec.value[i + 1 + j * NI] * FXE + vec.value[i + j * NI] * FXP;
  }
}

void Field::interpolatedFieldNorth(const Field &vec, const Grid &myGrid)
{
  forAllInteriorVCVs(NI, NJ)
  {
    double FYN = myGrid.YF[i + j * NI];
    double FYP = 1.0 - FYN;

    value[i + j * NI] = vec.value[i + (j + 1) * NI] * FYN + vec.value[i + j * NI] * FYP;
  }
}

void Field::computeEastMassFluxes(const Field &U)
{
  // For non constant density forAllInteriorUCVs
  forAllInterior(NI, NJ)
  {
    value[i + j * NI] = Se[i + j * NI] * density[i + j * NI] * U.value[i + j * NI];
  }
}

void Field::computeNorthMassFluxes(const Field &V)
{
  // For non constant density forAllInteriorVCVs
  forAllInterior(NI, NJ)
  {
    value[i + j * NI] = Sn[i + j * NI] * density[i + j * NI] * V.value[i + j * NI];
  }
}