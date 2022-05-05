#ifndef FINITEMATRIX_H
#define FINITEMATRIX_H

#include <iostream>
#include <vector>
#include <iomanip>
#include "ForAllOperators.hpp"
#include "Field.hpp"

using std::vector;

class FiniteMatrix
{
public:
  FiniteMatrix();
  virtual ~FiniteMatrix();

  typedef vector<vector<FiniteMatrix>> finiteMat;

  double value, aw, ae, as, an, ap, svalue;

  void print2dmat(finiteMat&);

  /* FUNCTIONS IN CASE U V AND P ARE IMPLEMENTED
  FiniteMatrix::finiteMat FiniteMatrix::intepolatedFieldEast(FiniteMatrix::finiteMat&, Grid&);
  FiniteMatrix::finiteMat FiniteMatrix::intepolatedFieldNorth(FiniteMatrix::finiteMat&, Grid&);

  Field::vectorField correctFaceVelocityEast(Field::vectorField&, Field::vectorField&,
Field::vectorField&, FiniteMatrix::finiteMat&, Grid&);
  Field::vectorField correctFaceVelocityNorth(Field::vectorField&, Field::vectorField&,
Field::vectorField&, FiniteMatrix::finiteMat&, Grid&);

  void correctEastMassFluxes(Field::vectorField&, Field::vectorField&, FiniteMatrix::finiteMat&);
  void correctNorthMassFluxes(Field::vectorField&, Field::vectorField&, FiniteMatrix::finiteMat&);
  */
  friend FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
  friend FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
  friend FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
  friend FiniteMatrix::finiteMat operator*(const double, const FiniteMatrix::finiteMat&);
};

#endif