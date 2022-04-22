#ifndef FIELD_H
#define FIELD_H

#include<vector>
#include<string>
#include <cmath>

#include "Grid.h"

using std::vector;
using std::string;

class Field
{
public:
  Field();

  //Overlad constructor
	Field(int NI, int NJ);
  virtual ~Field();
  
  enum Direction{west, east, south, north, point};

  double *value;

  int N, M, NI, NJ;

	double *X, *XC, *Y, *YC, *FXE, *FXP, *FYN, *FYP, *DXPtoE, *DYPtoN;
	double *Se, *Sn, *viscX, *viscY, *density, *volume;

  void getGridInfoPassed(Field &f, Grid &myGrid, double &viscX, double &viscY);
  void inletBoundaryCondition(Field& vec, Direction side, double bvalue);
  void laminarFlow(Field &vec, double m, double yMin, double yMax);
  void InitializeT(Field &vec, double &q, double xHotSpot, double yHotSpot, double xMin, double xMax);
  void InitializeZ(Field &vec, double T0hs, double r0hs, double xHotSpot, double yHotSpot, double xMax);
  void InitializeF(Field &vec, double xHotSpot, double xMin, double xMax);
  void initializeInternalField(Field &vec, double);
  void linearExtrapolateCondition(Field &vec, Direction);

  void interpolatedFieldEast(Field &interpolated, Field& vec, Grid& myGrid);
  void interpolatedFieldNorth(Field &interpolated, Field& vec, Grid& myGrid);

  void computeEastMassFluxes(Field& vec, Field& corrU);
  void computeNorthMassFluxes(Field& vec, Field& corrV);
};

#endif