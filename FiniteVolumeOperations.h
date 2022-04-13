#ifndef FINITEVOLUMEOPERATIONS_H
#define FINITEVOLUMEOPERATIONS_H

#include <math.h>
#include "Grid.h"
#include "Field.h"
#include "FiniteMatrix.h"
#include "ForAllOperators.h"

namespace fvm
{ 

  inline double plusupwind(double &v)
  {
    return (v + std::abs(v))/2;
  }

  inline double minusupwind(double &v)
  {
    return (v - std::abs(v))/2;
  }

  inline FiniteMatrix::finiteMat diffusiveTerm(Field::vectorField &vec)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));

    // Towards east side
    forAllInternalUCVs(vec)
    {
      APtemp[i][j].aevalue = -(vec[i][j].viscX * vec[i][j].Se) / vec[i][j].DXPtoE;
      APtemp[i + 1][j].awvalue = APtemp[i][j].aevalue;
    }

    for (long unsigned int j = 1; j < vec[0].size()-1; j++)
    {
      int i = 1;
      double DXPtoE = std::abs(vec[1][j].XC - vec[0][j].XC);
      double Se = std::abs(vec[0][j].Y - vec[0][j-1].Y);

      APtemp[i][j].awvalue = -(vec[i][j].viscX * Se) / DXPtoE;
      i = vec.size() - 2;
      APtemp[i][j].aevalue = -(vec[i][j].viscX * vec[i][j].Se) / vec[i][j].DXPtoE;
    }

    // Towards north side
    forAllInternalVCVs(vec)
    {
      APtemp[i][j].anvalue = -(vec[i][j].viscY * vec[i][j].Sn) / vec[i][j].DYPtoN;
      APtemp[i][j + 1].asvalue = APtemp[i][j].anvalue;
    }

    for (long unsigned int i = 1; i < vec.size()-1; i++)
    {
      int j = 1;
      double DYPtoN = std::abs(vec[i][1].YC - vec[i][0].YC);
      double Sn = std::abs(vec[i][0].X - vec[i-1][0].X);
      APtemp[i][j].asvalue = -(vec[i][j].viscY * Sn) / DYPtoN;
      j = vec[0].size() - 2;
      APtemp[i][j].anvalue = -(vec[i][j].viscY * vec[i][j].Sn) / vec[i][j].DYPtoN;
    }

    return APtemp;
  }

  inline FiniteMatrix::finiteMat convectiveTerm(Field::vectorField &vec, Field::vectorField &massFluxEast,
                                         Field::vectorField &massFluxNorth, double m)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));
    // Towards east side
    forAllInternalUCVs(vec)
    {
      APtemp[i][j].aevalue = m  *massFluxEast[i][j].value * vec[i][j].FXE;//* minusupwind(massFluxEast[i][j].value);
      APtemp[i + 1][j].awvalue = -m * massFluxEast[i][j].value * vec[i][j].FXP;// * plusupwind(massFluxEast[i][j].value);
      // double resultvalue = 0.0;
      // APtemp[i][j].svalue = APtemp[i][j].svalue + resultvalue;
    }

    for (long unsigned int j = 1; j < vec[0].size(); j++)
    {
      int i = 1;
      APtemp[i][j].awvalue =  -m *massFluxEast[i][j].value * vec[i][j].FXP;//* plusupwind(massFluxEast[i][j].value);
      i = vec.size() - 2;
      APtemp[i][j].aevalue =  m *massFluxEast[i][j].value * vec[i][j].FXE;//* minusupwind(massFluxEast[i][j].value);
    }

    // Towards north side
    forAllInternalVCVs(vec)
    {
      APtemp[i][j].anvalue = m * massFluxNorth[i][j].value * vec[i][j].FYN;//*minusupwind(massFluxNorth[i][j].value);
      APtemp[i][j + 1].asvalue = - m * massFluxNorth[i][j].value * vec[i][j].FYP;//* plusupwind(massFluxNorth[i][j].value);

      // double resultvalue = 0.0;
      // APtemp[i][j].svalue = APtemp[i][j].svalue + resultvalue;
    }

    for (long unsigned int i = 1; i < vec.size(); i++)
    {
      int j = 1;
      APtemp[i][j].asvalue =  - m * massFluxNorth[i][j].value * vec[i][j].FYP;//* plusupwind(massFluxNorth[i][j].value);
      j = vec[0].size() - 2;
      APtemp[i][j].anvalue =  m * massFluxNorth[i][j].value * vec[i][j].FYN;//* minusupwind(massFluxNorth[i][j].value);
    }

    return APtemp;
  }

  inline FiniteMatrix::finiteMat heatProduction(Field::vectorField &vec, Field::vectorField &Z, double q)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));

    forAllInternal(vec)
    {
      APtemp[i][j].svalue = q * Z[i][j].value * vec[i][j].Se * vec[i][j].Sn;
    }
    return APtemp;
  }

  inline FiniteMatrix::finiteMat intermidiateReaction(Field::vectorField &vec, Field::vectorField &vec2,
                                               Field::vectorField &T, double beta, double gamma)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));

    forAllInternal(vec)
    {
      APtemp[i][j].svalue = beta * beta * exp(beta * (T[i][j].value - 1.0) / (1.0 + gamma * (T[i][j].value - 1.0))) * vec[i][j].value * vec2[i][j].value * vec[i][j].Se * vec[i][j].Sn;
    }
    return APtemp;
  }

  inline FiniteMatrix::finiteMat zComsumptium(Field::vectorField &vec)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));

    forAllInternal(vec)
    {
      APtemp[i][j].svalue = vec[i][j].Se * vec[i][j].Sn * vec[i][j].value;
    }
    return APtemp;
  }

}
#endif