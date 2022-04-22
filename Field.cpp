#include "Field.h"

Field::Field()
{
}

Field::Field(int NI_, int NJ_) : NI(NI_), NJ(NJ_) // IF EVER NON CONST DENSITY THIS HAS TO CHANGE
{
  value = new double[NI * NJ];
  X = new double[NI * NJ];
  XC = new double[NI * NJ];
  FXE = new double[NI * NJ];
  FXP = new double[NI * NJ];
  Y = new double[NI * NJ];
  YC = new double[NI * NJ];
  FYN = new double[NI * NJ];
  FYP = new double[NI * NJ];
  DXPtoE = new double[NI * NJ];
  DYPtoN = new double[NI * NJ];
  Se = new double[NI * NJ];
  Sn = new double[NI * NJ];
  viscX = new double[NI * NJ];
  viscY = new double[NI * NJ];
  density = new double[NI * NJ];
  volume = new double[NI * NJ];
}

Field::~Field()
{
}

void Field::getGridInfoPassed(Field &f, Grid &myGrid, double &viscX, double &viscY)
{
  NI = myGrid.NI;
  NJ = myGrid.NJ;

  forAllN(NI, NJ)
  {
    int index = i + j * NI;
    f.X[index] = myGrid.X[index];
    f.Y[index] = myGrid.Y[index];
    f.XC[index] = myGrid.XC[index];
    f.YC[index] = myGrid.YC[index];
    f.FXE[index] = myGrid.XF[index];
    f.FYN[index] = myGrid.YF[index];
    f.FXP[index] = 1.0 - myGrid.XF[index];
    f.FYP[index] = 1.0 - myGrid.YF[index];

    f.viscX[index] = viscX;
    f.viscY[index] = viscY;
  }

  forAllInterior(NI, NJ)
  {
    int index = i + j * NI;
    f.DXPtoE[index] = std::abs(myGrid.XC[index + 1] - myGrid.XC[index]);
    f.DYPtoN[index] = std::abs(myGrid.YC[index + NI] - myGrid.YC[index]);

    f.Se[index] = std::abs(myGrid.Y[index] - myGrid.Y[index - NI]);
    f.Sn[index] = std::abs(myGrid.X[index] - myGrid.X[index - 1]);

    f.volume[index] = f.Se[index] * f.Sn[index];
  }
}

void Field::inletBoundaryCondition(Field &vec, Direction side, double bvalue)
{
  NI = vec.NI;
  NJ = vec.NJ;

  if (side == west)
  {
    forWBoundary(NI, NJ)
    {
      vec.value[i + j * NI] = bvalue;
    }
  }
  else if (side == east)
  {
    forEBoundary(NI, NJ)
    {
      vec.value[i + j * NI] = bvalue;
    }
  }
  else if (side == south)
  {
    forSBoundary(NI, NJ)
    {
      vec.value[i + j * NI] = bvalue;
    }
  }
  else if (side == north)
  {
    forNBoundary(NI, NJ)
    {
      vec.value[i + j * NI] = bvalue;
    }
  }
}

void Field::laminarFlow(Field &vec, double m, double yMin, double yMax)
{

  forAllN(NI, NJ)
  {
    double distanceToBot = vec.YC[i + j * NI] - yMin;
    double distanceToTop = yMax - vec.YC[i + j * NI];
    // std::cout << i << " " << j << " " << distanceToBot << " " << distanceToTop << std::endl;
    vec.value[i + j * NI] = 6.0 * m / std::abs(m) * distanceToBot * distanceToTop;
  }
}

void Field::InitializeT(Field &vec, double &q,
                        double xHotSpot, double yHotSpot, double xMin, double xMax)
{
  NI = vec.NI;
  NJ = vec.NJ;

  forAllN(NI, NJ)
  {
    double pos = vec.XC[i + j * NI] - xHotSpot;
    double simmetricpos = vec.XC[i + j * NI] - (xMin + xMax - xHotSpot);

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));

    vec.value[i + j * NI] = q * 1.5 * tanV;
    if (std::abs(pos) > std::abs(simmetricpos))
    {
      tanV = (atan(simmetricpos) - atan(-(xMin - xMax) / 2)) /
             (atan((xMin - xMax) / 2) - atan(-(xMin - xMax) / 2));
      vec.value[i + j * NI] = std::max(q, q * 1.5 * tanV);
    }
  }
}

void Field::InitializeF(Field &vec, double xHotSpot, double xMin, double xMax)
{
  NI = vec.NI;
  NJ = vec.NJ;

  forAllN(NI, NJ)
  {
    double pos = vec.XC[i + j * NI] - xHotSpot;

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));
    vec.value[i + j * NI] = std::max(1.0 - tanV, 0.0);
  }
}

void Field::InitializeZ(Field &vec, double T0hs, double r0hs,
                        double xHotSpot, double yHotSpot, double xMax)
{
  NI = vec.NI;
  NJ = vec.NJ;

  forAllN(NI, NJ)
  {
    double rr = sqrt(pow(vec.XC[i + j * NI] - xHotSpot, 2) +
                     pow(vec.YC[i + j * NI] - yHotSpot, 2));

    vec.value[i + j * NI] = T0hs * exp(-rr / r0hs);
  }
}

void Field::initializeInternalField(Field &vec, double val)
{
  NI = vec.NI;
  NJ = vec.NJ;

  forAllInterior(NI, NJ)
  {
    vec.value[i + j * NI] = val;
  }
}

void Field::linearExtrapolateCondition(Field &vec, Field::Direction wallname)
{
  NI = vec.NI;
  NJ = vec.NJ;

  if (wallname == west)
  {
    forWBoundary(NI, NJ)
    {
      int index = i + j * NI;
      vec.value[index] = vec.value[index + 1] + (vec.value[index + 1] - vec.value[index + 2]) * vec.FXE[index + 1];
    }
  }
  else if (wallname == east)
  {
    forEBoundary(NI, NJ)
    {
      int index = i + j * NI;
      vec.value[index] = vec.value[index - 1] + (vec.value[index - 1] - vec.value[index - 2]) * vec.FXP[index - 1];
    }
  }
  else if (wallname == south)
  {
    forSBoundary(NI, NJ)
    {
      int index =i + j * NI;
      vec.value[index] = vec.value[index + NI] + (vec.value[index + NI] - vec.value[index + 2 * NI]) * vec.FYN[index + NI];
    }
  }
  else if (wallname == north)
  {
    forNBoundary(NI, NJ)
    {
      int index = i + j * NI;
      vec.value[index] = vec.value[index - NI] + (vec.value[index - NI] - vec.value[index - 2 * NI]) * vec.FYP[index - NI];
    }
  }
}

void Field::interpolatedFieldEast(Field &interpolated, Field &vec, Grid &myGrid)
{
  NI = vec.NI;
  NJ = vec.NJ;

  forAllInteriorUCVs(NI, NJ)
  {
    double FXE = myGrid.XF[i + j * NI];
    double FXP = 1.0 - FXE;

    interpolated.value[i + j * NI] = vec.value[i + 1 + j * NI] * FXE + vec.value[i + j * NI] * FXP;
  }
}

void Field::interpolatedFieldNorth(Field &interpolated, Field &vec, Grid &myGrid)
{
  NI = vec.NI;
  NJ = vec.NJ;

  forAllInteriorVCVs(NI, NJ)
  {
    double FYN = myGrid.YF[i + j * NI];
    double FYP = 1.0 - FYN;

    interpolated.value[i + j * NI] = vec.value[i + (j + 1) * NI] * FYN + vec.value[i + j * NI] * FYP;
  }
}

void Field::computeEastMassFluxes(Field &vec, Field &corrU)
{
  NI = vec.NI;
  NJ = vec.NJ;

  // For non constant density forAllInteriorUCVs
  forAllInterior(NI, NJ)
  {
    vec.value[i + j * NI] = vec.Se[i + j * NI] * vec.density[i + j * NI] * corrU.value[i + j * NI];
  }
}

void Field::computeNorthMassFluxes(Field &vec, Field &corrV)
{
  NI = vec.NI;
  NJ = vec.NJ;

  // For non constant density forAllInteriorVCVs
  forAllInterior(NI, NJ)
  {
    vec.value[i + j * NI] = vec.Sn[i + j * NI] * vec.density[i + j * NI] * corrV.value[i + j * NI];
  }
}