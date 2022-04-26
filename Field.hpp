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
  Field();

  // Overlad constructor
  Field(int NI, int NJ);
  virtual ~Field();

  enum Direction
  {
    west,
    east,
    south,
    north,
    point
  };

  int NI, NJ;

  double *value;

  double *X, *XC, *FXE, *FXP, *Y, *YC, *FYN, *FYP, *DXPtoE, *DYPtoN;
  double *Se, *Sn, *viscX, *viscY, *density, *volume;

  void getGridInfoPassed(Field &f, Grid &myGrid, double &viscX, double &viscY);
  void inletBoundaryCondition(Field &vec, Direction side, double bvalue);
  void laminarFlow(Field &vec, double m, double yMin, double yMax);
  void InitializeT(Field &vec, double &q, double xHotSpot, double yHotSpot, double xMin, double xMax);
  void InitializeZ(Field &vec, double T0hs, double r0hs, double xHotSpot, double yHotSpot, double xMax);
  void InitializeF(Field &vec, double xHotSpot, double xMin, double xMax);
  void initializeInternalField(Field &vec, double);
  void linearExtrapolateCondition(Field &vec, Direction);

  void interpolatedFieldEast(Field &interpolated, Field &vec, Grid &myGrid);
  void interpolatedFieldNorth(Field &interpolated, Field &vec, Grid &myGrid);

  void computeEastMassFluxes(Field &vec, Field &corrU);
  void computeNorthMassFluxes(Field &vec, Field &corrV);

  friend void swap(Field &first, Field &second)
  {
    using std::swap;

    swap(first.NI, second.NI);
    swap(first.NJ, second.NJ);

    swap(first.X, second.X);
    swap(first.XC, second.XC);
    
    swap(first.Y, second.Y);
    swap(first.YC, second.YC);
    
    swap(first.FXE, second.FXE);
    swap(first.FXP, second.FXP);
    
    swap(first.FYN, second.FYN);
    swap(first.FYP, second.FYP);

    swap(first.DXPtoE, second.DXPtoE);
    swap(first.DYPtoN, second.DYPtoN);

    swap(first.Se, second.Se);
    swap(first.Sn, second.Sn);
    
    swap(first.viscX, second.viscX);
    swap(first.viscY, second.viscY);
    
    swap(first.density, second.density);
    swap(first.volume, second.volume);

    swap(first.value, second.value);
  };

  Field operator=(Field other)
  {
    swap(*this, other);

    return *this;
  };
};

#endif