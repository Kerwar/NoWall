#ifndef FINITEVOLUMEOPERATIONS_H
#define FINITEVOLUMEOPERATIONS_H

#include <math.h>

#include "Field.hpp"
#include "FiniteMatrix.hpp"
#include "ForAllOperators.hpp"

namespace fvm {

inline double plusupwind(const double &v) { return (v + std::abs(v)) / 2; };

inline double minusupwind(const double &v) { return (v - std::abs(v)) / 2; };

inline FiniteMatrix::finiteMat diffusiveTerm(const Field &vec) {
  int NI = vec.NI;
  int NJ = vec.NJ;

  FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

  // Towards east side
  forAllInteriorUCVs(vec.NI, vec.NJ) {
    int index = i + j * NI;
    APtemp[i][j].ae = -(vec.viscX[index] * vec.Se[j]) / vec.DXPtoE[i];
    APtemp[i + 1][j].aw = APtemp[i][j].ae;
  }

  for (int j = 1; j < NJ - 1; j++) {
    int i = 1;
    double DXPtoE = std::abs(vec.XC[1] - vec.XC[0]);

    APtemp[i][j].aw = -(vec.viscX[i + j * NI] * vec.Se[j]) / DXPtoE;
    i = NI - 2;
    APtemp[i][j].ae = -(vec.viscX[i + j * NI] * vec.Se[j]) / vec.DXPtoE[i];
  }

  // Towards north side
  forAllInteriorVCVs(NI, NJ) {
    int index = i + j * NI;
    APtemp[i][j].an = -(vec.viscY[index] * vec.Sn[i]) / vec.DYPtoN[j];
    APtemp[i][j + 1].as = APtemp[i][j].an;
  }

  for (int i = 1; i < NI - 1; i++) {
    int j = 1;
    double DYPtoN = std::abs(vec.YC[1] - vec.YC[0]);
    APtemp[i][j].as = -(vec.viscY[i + j * NI] * vec.Sn[i]) / DYPtoN;
    j = NJ - 2;
    APtemp[i][j].an = -(vec.viscY[i + j * NI] * vec.Sn[i]) / vec.DYPtoN[j];
  }

  return APtemp;
}

inline FiniteMatrix::finiteMat diffusiveTermS(const Field &vec,
                                              const vector<double> &TWall,
                                              const double &exCte) {
  int NI = vec.NI;
  int NJ = vec.NJ;

  FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

  APtemp = diffusiveTerm(vec);

  for (int i = 1; i < NI - 1; i++) {
    int j = 1;
    APtemp[i][j].as = 0;
    if (std::abs(TWall[i]) >= 10e-8) {
      APtemp[i][j].ap = exCte * vec.Sn[i] * vec.viscY[i + j * NI];
      APtemp[i][j].svalue =
          exCte * vec.Sn[i] * vec.viscY[i + j * NI] * TWall[i];
    }
  }

  return APtemp;
}

inline FiniteMatrix::finiteMat diffusiveTermN(const Field &vec,
                                              const vector<double> &TWall,
                                              const double &exCte) {
  int NI = vec.NI;
  int NJ = vec.NJ;

  FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

  APtemp = diffusiveTerm(vec);

  for (int i = 1; i < NI - 1; i++) {
    int j = NJ - 2;
    APtemp[i][j].an = 0;
    if (std::abs(TWall[i]) >= 10e-8) {
      APtemp[i][j].ap = exCte * vec.Sn[i] * vec.viscY[i + j * NI];
      APtemp[i][j].svalue =
          exCte * vec.Sn[i] * vec.viscY[i + j * NI] * TWall[i];
    }
    // if (Tdiff[i] != 0)
    // std::cout << APtemp[i][j].svalue << endl;
  }

  return APtemp;
}
inline FiniteMatrix::finiteMat convectiveTerm(const Field &vec,
                                              const Field &massFluxEast,
                                              const Field &massFluxNorth,
                                              const double &m) {
  int NI = vec.NI;
  int NJ = vec.NJ;

  FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));
  // Towards east side
  forAllInteriorUCVs(NI, NJ) {
    int index = i + j * NI;
    // APtemp[i][j].ae = m * minusupwind(massFluxEast.value[index]);
    APtemp[i][j].ae =
        m * massFluxEast.value[index] *
        vec.FXE[i];  //* minusupwind(massFluxEast.value[i + j *NI]);
    // APtemp[i + 1][j].aw = m * plusupwind(massFluxEast.value[index]);
    APtemp[i + 1][j].aw =
        -m * massFluxEast.value[index] *
        vec.FXP[i];  // * plusupwind(massFluxEast.value[i + j *NI]);
    // double resultvalue = 0.0;
    // APtemp[i][j].ap += ;
  }

  for (int j = 1; j < NJ; j++) {
    int i = 1;
    // APtemp[i][j].ae = m * minusupwind(massFluxEast.value[i + j *NI]);
    APtemp[i][j].aw =
        -m * massFluxEast.value[i + j * NI] *
        vec.FXP[0];  //* plusupwind(massFluxEast.value[i + j *NI]);
    i = NI - 2;
    // APtemp[i + 1][j].aw = m * plusupwind(massFluxEast.value[i + j *NI]);
    APtemp[i][j].ae =
        m * massFluxEast.value[i + j * NI] *
        vec.FXE[i];  //* minusupwind(massFluxEast.value[i + j *NI]);
  }

  // Towards north side
  forAllInteriorVCVs(NI, NJ) {
    int index = i + j * NI;
    // APtemp[i][j].an = m * minusupwind(massFluxNorth.value[index]);
    APtemp[i][j].an =
        m * massFluxNorth.value[index] *
        vec.FYN[j];  //*minusupwind(massFluxNorth.value[i + j *NI]);
    // APtemp[i + 1][j].as = m * plusupwind(massFluxNorth.value[index]);
    APtemp[i][j + 1].as =
        -m * massFluxNorth.value[index] *
        vec.FYP[j];  //* plusupwind(massFluxNorth.value[i + j *NI]);

    // double resultvalue = 0.0;
    // APtemp[i][j].svalue = APtemp[i][j].svalue + resultvalue;
  }

  for (int i = 1; i < NI; i++) {
    int j = 1;
    // APtemp[i][j].an = m * minusupwind(massFluxNorth.value[i + j * NI]);
    APtemp[i][j].as =
        -m * massFluxNorth.value[i + j * NI] *
        vec.FYP[0];  //* plusupwind(massFluxNorth.value[i + j *NI]);
    j = NJ - 2;
    // APtemp[i + 1][j].as = m * plusupwind(massFluxNorth.value[i + j * NI]);
    APtemp[i][j].an =
        m * massFluxNorth.value[i + j * NI] *
        vec.FYN[j];  //* minusupwind(massFluxNorth.value[i + j *NI]);
  }

  return APtemp;
}

inline FiniteMatrix::finiteMat heatProduction(const Field &vec, const Field &Z,
                                              const double &q) {
  int NI = vec.NI;
  int NJ = vec.NJ;

  FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

  forAllInterior(NI, NJ) {
    int index = i + j * NI;
    APtemp[i][j].svalue = q * Z[index] * vec.volume[index];
  }
  return APtemp;
}

inline FiniteMatrix::finiteMat intermidiateReaction(const Field &vec,
                                                    const Field &vec2,
                                                    const Field &T,
                                                    const double &beta,
                                                    const double &gamma) {
  int NI = vec.NI;
  int NJ = vec.NJ;
  FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

  forAllInterior(NI, NJ) {
    int index = i + j * NI;
    double expPortion = exp(beta * (T.value[index] - 1.0) /
                            (1.0 + gamma * (T.value[index] - 1.0)));
    APtemp[i][j].svalue = beta * beta * expPortion * vec[index] *
                          vec2[index] * vec.volume[index];
  }
  return APtemp;
}

inline FiniteMatrix::finiteMat zComsumptium(const Field &vec) {
  int NI = vec.NI;
  int NJ = vec.NJ;

  FiniteMatrix::finiteMat APtemp(NI, vector<FiniteMatrix>(NJ));

  forAllInterior(NI, NJ) {
    APtemp[i][j].svalue = vec.volume[i + j * NI] * vec[i + j * NI];
  }
  return APtemp;
}

}  // namespace fvm
#endif
