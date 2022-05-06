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

  double *X, *XC, *FXE, *FXP, *Y, *YC, *FYN, *FYP, *DXPtoE, *DYPtoN;
  double *Se, *Sn, *viscX, *viscY, *density, *volume;

  // Field();

  // Overlad constructor
  Field(const int &_NI,const  int& _NJ);
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

  // friend void swap(Field &first, Field &second)
  // {
  //   using std::swap;

  //   swap(first.X, second.X);
  //   swap(first.XC, second.XC);
    
  //   swap(first.Y, second.Y);
  //   swap(first.YC, second.YC);
    
  //   swap(first.FXE, second.FXE);
  //   swap(first.FXP, second.FXP);
    
  //   swap(first.FYN, second.FYN);
  //   swap(first.FYP, second.FYP);

  //   swap(first.DXPtoE, second.DXPtoE);
  //   swap(first.DYPtoN, second.DYPtoN);

  //   swap(first.Se, second.Se);
  //   swap(first.Sn, second.Sn);
    
  //   swap(first.viscX, second.viscX);
  //   swap(first.viscY, second.viscY);
    
  //   swap(first.density, second.density);
  //   swap(first.volume, second.volume);

  //   swap(first.value, second.value);
  // };

  // constexpr Field operator=(Field other)
  // {
  //   swap(*this, other);

  //   return *this;
  // };
};

#endif