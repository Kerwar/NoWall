#include "Field.hpp"

Field::Field(const int &_NI, const int &_NJ)
    : NI(_NI),
      NJ(_NJ),
      value(vector<double>(NI * NJ, 0)),
      FXP(vector<double>(NI, 0)),
      FYP(vector<double>(NJ, 0)),
      DXPtoE(vector<double>(NI, 0)),
      DYPtoN(vector<double>(NJ, 0)),
      Se(vector<double>(NJ, 0)),
      Sn(vector<double>(NI, 0)),
      viscX(vector<double>(NI * NJ, 0)),
      viscY(vector<double>(NI * NJ, 0)),
      density(vector<double>(NI * NJ, 0)),
      volume(vector<double>(NI * NJ, 0))
// IF EVER NON CONST DENSITY THIS HAS TO CHANGE
{}

Field::~Field() {}

void Field::getGridInfoPassed(const Grid &myGrid,const  double &viscX_,
                              const double &viscY_) {
  X = myGrid.X;
  Y = myGrid.Y;
  XC = myGrid.XC;
  YC = myGrid.YC;
  FXE = myGrid.XF;
  FYN = myGrid.YF;
  forAllN(NI, NJ) {
    int index = i + j * NI;
    FXP[i] = 1.0 - myGrid.XF[i];
    FYP[j] = 1.0 - myGrid.YF[j];
    density[index] = 1.0;
    viscX[index] = viscX_;
    viscY[index] = viscY_;
  }

  forAllInterior(NI, NJ) {
    int index = i + j * NI;
    DXPtoE[i] = std::abs(myGrid.XC[i + 1] - myGrid.XC[i]);
    DYPtoN[j] = std::abs(myGrid.YC[j + 1] - myGrid.YC[j]);

    Se[j] = std::abs(myGrid.Y[j] - myGrid.Y[j - 1]);
    Sn[i] = std::abs(myGrid.X[i] - myGrid.X[i - 1]);

    volume[index] = Se[j] * Sn[i];
  }
}

void Field::inletBoundaryCondition(Direction side, double bvalue) {
  if (side == west) {
    forWBoundary(NI, NJ) { value[i + j * NI] = bvalue; }
  } else if (side == east) {
    forEBoundary(NI, NJ) { value[i + j * NI] = bvalue; }
  } else if (side == south) {
    forSBoundary(NI, NJ) { value[i + j * NI] = bvalue; }
  } else if (side == north) {
    forNBoundary(NI, NJ) { value[i + j * NI] = bvalue; }
  }
}

void Field::laminarFlow(const double &m, const double &yMin,
                        const double &yMax) {
  forAllN(NI, NJ) {
    double distanceToBot = YC[j] - yMin;
    double distanceToTop = yMax - YC[j];
    // std::cout << i << " " << j << " " << distanceToBot << " " <<
    // distanceToTop << std::endl;
    value[i + j * NI] = 6.0 * m / std::abs(m) * distanceToBot * distanceToTop;
  }
}

void Field::InitializeT(const double &q, const double &xHotSpot,
                        const double &xMin, const double &xMax) {
  forAllN(NI, NJ) {
    double pos = XC[i] - xHotSpot;
    double simmetricpos = XC[i] - (xMin + xMax - xHotSpot);

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));

    value[i + j * NI] = q * 1.5 * tanV;
    if (std::abs(pos) > std::abs(simmetricpos)) {
      tanV = (atan(simmetricpos) - atan(-(xMin - xMax) / 2)) /
             (atan((xMin - xMax) / 2) - atan(-(xMin - xMax) / 2));
      value[i + j * NI] = std::max(q, q * 1.5 * tanV);
    }
  }
}

void Field::InitializeF(const double &xHotSpot, const double &xMin,
                        const double &xMax) {
  forAllN(NI, NJ) {
    double pos = XC[i] - xHotSpot;

    double tanV = (atan(pos) - atan(-(xMax - xMin) / 2)) /
                  (atan((xMax - xMin) / 2) - atan(-(xMax - xMin) / 2));
    value[i + j * NI] = std::max(1.0 - tanV, 0.0);
  }
}

void Field::InitializeZ(const double &T0hs, const double &r0hs,
                        const double &xHotSpot, const double &yHotSpot) {
  forAllN(NI, NJ) {
    double rr = sqrt(pow(XC[i] - xHotSpot, 2) + pow(YC[j] - yHotSpot, 2));

    value[i + j * NI] = T0hs * exp(-rr / r0hs);
  }
}

void Field::initializeInternalField(double val) {
  forAllInterior(NI, NJ) { value[i + j * NI] = val; }
}

void Field::linearExtrapolateCondition(const Direction &wallname) {
  if (wallname == west) {
    forWBoundary(NI, NJ) {
      int index = i + j * NI;
      value[index] =
          value[index + 1] + (value[index + 1] - value[index + 2]) * FXE[i + 1];
    }
  } else if (wallname == east) {
    forEBoundary(NI, NJ) {
      int index = i + j * NI;
      value[index] =
          value[index - 1] + (value[index - 1] - value[index - 2]) * FXP[i - 1];
    }
  } else if (wallname == south) {
    forSBoundary(NI, NJ) {
      int index = i + j * NI;
      value[index] = value[index + NI] +
                     (value[index + NI] - value[index + 2 * NI]) * FYN[j + 1];
    }
  } else if (wallname == north) {
    forNBoundary(NI, NJ) {
      int index = i + j * NI;
      value[index] = value[index - NI] +
                     (value[index - NI] - value[index - 2 * NI]) * FYP[j - 1];
    }
  }
}

void Field::interpolatedFieldEast(const Field &vec) {
  forAllInteriorUCVs(NI, NJ) {
    value[i + j * NI] =
        vec.value[i + 1 + j * NI] * FXE[i] + vec.value[i + j * NI] * FXP[i];
  }
}

void Field::interpolatedFieldNorth(const Field &vec) {
  forAllInteriorVCVs(NI, NJ) {
    value[i + j * NI] =
        vec.value[i + (j + 1) * NI] * FYN[j] + vec.value[i + j * NI] * FYP[j];
  }
}

void Field::computeEastMassFluxes(const Field &U) {
  // For non constant density forAllInteriorUCVs
  forAllInterior(NI, NJ) {
    value[i + j * NI] = Se[j] * density[i + j * NI] * U.value[i + j * NI];
  }
}

void Field::computeNorthMassFluxes(const Field &V) {
  // For non constant density forAllInteriorVCVs
  forAllInterior(NI, NJ) {
    value[i + j * NI] = Sn[i] * density[i + j * NI] * V.value[i + j * NI];
  }
}