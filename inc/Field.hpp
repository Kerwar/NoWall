#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>
#include <cmath>

#include "Grid.hpp"

using std::string;
using std::vector;

class Field
{
public:
  const int NI, NJ;

  double *value;

  std::shared_ptr<double[]> X, XC, FXE, FYN, Y, YC;

  double *FXP, *FYP, *DXPtoE, *DYPtoN;
  double *Se, *Sn, *viscX, *viscY, *density, *volume;

  // Field();

  // Overlad constructor
  Field(const int &_NI, const int &_NJ);
  virtual ~Field();

  enum Direction
  {
    west,
    east,
    south,
    north,
    point
  };

  void getGridInfoPassed(const Grid &myGrid, double &viscX, double &viscY);
  void inletBoundaryCondition(Direction side, double bvalue);
  void laminarFlow(const double &m, const double &yMin, const double &yMax);
  void InitializeT(const double &q, const double &xHotSpot, const double &xMin, const double &xMax);
  void InitializeZ(const double &T0hs, const double &r0hs, const double &xHotSpot, const double &yHotSpot);
  void InitializeF(const double &xHotSpot, const double &xMin, const double &xMax);
  void initializeInternalField(double);
  void linearExtrapolateCondition(const Direction &wallname);

  void interpolatedFieldEast(const Field &vec);
  void interpolatedFieldNorth(const Field &vec);

  void computeEastMassFluxes(const Field &U);
  void computeNorthMassFluxes(const Field &V);

  double operator[](const int &i) { return value[i]; };
};

#endif