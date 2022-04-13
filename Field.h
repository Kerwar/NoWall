#ifndef FIELD_H
#define FIELD_H

#include<vector>
#include<string>

#include "Grid.h"
#include "ForAllOperators.h"

using std::vector;
using std::string;

class Field
{
public:
  Field();

  //Overlad constructor
	Field(int&, int&);
  virtual ~Field();
  
  enum Direction{west, east, south, north, point};
  typedef vector<Field> vec1dfield;
  typedef vector<vector<Field>> vectorField;

  double value;

  int N, M, NI, NJ; 
	double X, XC, Y, YC, FXE, FXP, FYN, FYP, DXPtoE, DYPtoN;
	double Se, Sn, viscX, viscY, density, volume;

  void getGridInfoPassed(Field::vectorField&, Grid&, double&, double&);
  void inletBoundaryCondition(Field::vectorField&, Direction, double);
  void laminarFlow(Field::vectorField&, double, double, double);
  void InitializeT(Field::vectorField &vec, double &q, double xHotSpot, double yHotSpot, double xMin, double xMax);
  void InitializeZ(Field::vectorField &vec, double T0hs, double r0hs, double xHotSpot, double yHotSpot, double xMax);
  void InitializeF(Field::vectorField &vec, double xHotSpot, double xMin, double xMax);
  void initializeInternalField(Field::vectorField&, double);
  void linearExtrapolateCondition(Field::vectorField&, Direction);

  Field::vectorField interpolatedFieldEast(Field::vectorField& vec, Grid& myGrid);
  Field::vectorField interpolatedFieldNorth(Field::vectorField& vec, Grid& myGrid);

  void computeEastMassFluxes(Field::vectorField& vec, Field::vectorField& corrU);
  void computeNorthMassFluxes(Field::vectorField& vec, Field::vectorField& corrV);

  friend Field operator+(const Field&, const Field&);
  friend Field operator-(const Field&, const Field&);
  friend Field operator*(const double, const Field&);
  friend Field operator*(const Field&, const double);
  friend Field::vectorField operator+(const Field::vectorField&, const Field::vectorField&);
  friend Field::vectorField operator-(const Field::vectorField&, const Field::vectorField&);
  friend Field::vectorField operator*(const double, const Field::vectorField&);
  friend Field::vectorField operator*(const Field::vectorField&, const double);
};

#endif