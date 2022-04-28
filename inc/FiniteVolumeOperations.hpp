#ifndef FINITEVOLUMEOPERATIONS_H
#define FINITEVOLUMEOPERATIONS_H

#include <math.h>
#include "Field.hpp"
#include "FiniteMatrix.hpp"
#include "ForAllOperators.hpp"

namespace fvm
{

  inline double plusupwind(double &v)
  {
    return (v + std::abs(v)) / 2;
  }

  inline double minusupwind(double &v)
  {
    return (v - std::abs(v)) / 2;
  }

  inline FiniteMatrix::finiteMat diffusiveTerm(Field &vec)
  {
    int NI = vec.NI;
    int NJ = vec.NJ;

    FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

    // Towards east side
    forAllInteriorUCVs(vec.NI, vec.NJ)
    {
      int index = i + j * NI;
      APtemp[i][j].aevalue = -(vec.viscX[index] * vec.Se[index]) / vec.DXPtoE[index];
      APtemp[i + 1][j].awvalue = APtemp[i][j].aevalue;
    }

    for (int j = 1; j < NJ - 1; j++)
    {
      int i = 1;
      double DXPtoE = std::abs(vec.XC[i + j * NI] - vec.XC[i -1 +j * NI]);
      double Se = std::abs(vec.Y[j * NI] - vec.Y[(j - 1) * NI]);

      APtemp[i][j].awvalue = -(vec.viscX[i + j * NI] * Se) / DXPtoE;
      i = NI - 2;
      APtemp[i][j].aevalue = -(vec.viscX[i + j * NI] * vec.Se[i + j * NI]) / vec.DXPtoE[i + j * NI];
    }

    // Towards north side
    forAllInteriorVCVs(NI, NJ)
    { 
      int index = i + j * NI;
      APtemp[i][j].anvalue = -(vec.viscY[index] * vec.Sn[index]) / vec.DYPtoN[index];
      APtemp[i][j + 1].asvalue = APtemp[i][j].anvalue;
    }

    for (int i = 1; i < NI - 1; i++)
    {
      int j = 1;
      double DYPtoN = std::abs(vec.YC[i + NI] - vec.YC[i]);
      double Sn = std::abs(vec.X[i] - vec.X[i - 1]);
      APtemp[i][j].asvalue = -(vec.viscY[i + j * NI] * Sn) / DYPtoN;
      j = NJ - 2;
      APtemp[i][j].anvalue = -(vec.viscY[i + j * NI] * vec.Sn[i + j * NI]) / vec.DYPtoN[i + j * NI];
    }

    return APtemp;
  }

  inline FiniteMatrix::finiteMat convectiveTerm(Field &vec, Field &massFluxEast,
                                                Field &massFluxNorth, double m)
  {
    int NI = vec.NI;
    int NJ = vec.NJ;

    FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));
    // Towards east side
    forAllInteriorUCVs(NI, NJ)
    {
      int index = i + j * NI;
      APtemp[i][j].aevalue = m * massFluxEast.value[index] * vec.FXE[index];      //* minusupwind(massFluxEast.value[i + j *NI]);
      APtemp[i + 1][j].awvalue = -m * massFluxEast.value[index] * vec.FXP[index]; // * plusupwind(massFluxEast.value[i + j *NI]);
      // double resultvalue = 0.0;
      // APtemp[i][j].svalue = APtemp[i][j].svalue + resultvalue;
    }

    for (int j = 1; j < NJ; j++)
    {
      int i = 1;
      APtemp[i][j].awvalue = -m * massFluxEast.value[i + j * NI] * vec.FXP[i + j * NI]; //* plusupwind(massFluxEast.value[i + j *NI]);
      i = NI - 2;
      APtemp[i][j].aevalue = m * massFluxEast.value[i + j * NI] * vec.FXE[i + j * NI]; //* minusupwind(massFluxEast.value[i + j *NI]);
    }

    // Towards north side
    forAllInteriorVCVs(NI, NJ)
    {
      int index = i + j * NI;
      APtemp[i][j].anvalue = m * massFluxNorth.value[index] * vec.FYN[index];      //*minusupwind(massFluxNorth.value[i + j *NI]);
      APtemp[i][j + 1].asvalue = - m * massFluxNorth.value[index] * vec.FYP[index]; //* plusupwind(massFluxNorth.value[i + j *NI]);

      // double resultvalue = 0.0;
      // APtemp[i][j].svalue = APtemp[i][j].svalue + resultvalue;
    }

    for (int i = 1; i < NI; i++)
    {
      int j = 1;
      APtemp[i][j].asvalue = - m * massFluxNorth.value[i + j * NI] * vec.FYP[i + j * NI]; //* plusupwind(massFluxNorth.value[i + j *NI]);
      j = NJ - 2;
      APtemp[i][j].anvalue = m * massFluxNorth.value[i + j * NI] * vec.FYN[i + j * NI]; //* minusupwind(massFluxNorth.value[i + j *NI]);
    }

    return APtemp;
  }

  inline FiniteMatrix::finiteMat heatProduction(Field &vec, Field &Z, double q)
  {
    int NI = vec.NI;
    int NJ = vec.NJ;

    FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

    forAllInterior(NI, NJ)
    {
      int index = i + j * NI;
      APtemp[i][j].svalue = q * Z.value[index] * vec.Se[index] * vec.Sn[index];
    }
    return APtemp;
  }

  inline FiniteMatrix::finiteMat intermidiateReaction(Field &vec, Field &vec2,
                                                      Field &T, double beta, double gamma)
  {
    int NI = vec.NI;
    int NJ = vec.NJ;
    FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

    forAllInterior(NI, NJ)
    {
      int index = i + j * NI;
      double expPortion =  exp(beta * (T.value[index] - 1.0) / (1.0 + gamma * (T.value[index] - 1.0)));
      APtemp[i][j].svalue = beta * beta * expPortion * vec.value[index] * vec2.value[index] * vec.Se[index] * vec.Sn[index];
    }
    return APtemp;
  }

  inline FiniteMatrix::finiteMat zComsumptium(Field &vec)
  {
    int NI = vec.NI;
    int NJ = vec.NJ;

    FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

    forAllInterior(NI, NJ)
    {
      int index = i + j * NI;
      APtemp[i][j].svalue = vec.Se[index] * vec.Sn[index] * vec.value[index];
    }
    return APtemp;
  }

}
#endif