#include "FiniteMatrix.h"

FiniteMatrix::FiniteMatrix() : value(0.0), awvalue(0.0), aevalue(0.0), asvalue(0.0), anvalue(0.0), apvalue(0.0), svalue(0.0)
{
}

FiniteMatrix::~FiniteMatrix()
{
}

FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat &lhs, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat result(lhs);

  forAllInternal(lhs)
  {
    result[i][j].value += rhs[i][j].value;
    result[i][j].awvalue += rhs[i][j].awvalue;
    result[i][j].aevalue += rhs[i][j].aevalue;
    result[i][j].asvalue += rhs[i][j].asvalue;
    result[i][j].anvalue += rhs[i][j].anvalue;
    result[i][j].svalue += rhs[i][j].svalue;
  }

  return result;
}

/* FUNCTIONS IN CASE U V AND P ARE IMPLEMENTED
FiniteMatrix::finiteMat FiniteMatrix::intepolatedFieldEast(FiniteMatrix::finiteMat& vec, Grid& myGrid)
{
  FiniteMatrix::finiteMat temp(vec.size(), vector<FiniteMatrix>(vec[0].size()));

  forAllInternalUCVs(vec)
  {
    double FXE = myGrid.XF[i][j];
    double FXP = 1.0 - FXE;

    temp[i][j].value = vec[i+1][j].value * FXE + vec[i][j].value * FXP;
  }
  return temp;
}

FiniteMatrix::finiteMat FiniteMatrix::intepolatedFieldNorth(FiniteMatrix::finiteMat& vec, Grid& myGrid)
{
  FiniteMatrix::finiteMat temp(vec.size(), vector<FiniteMatrix>(vec[0].size()));

  forAllInternalVCVs(vec)
  {
    double FYN = myGrid.YF[i][j];
    double FYP = 1.0 - FYN;

    temp[i][j].value = vec[i][j+1].value * FYN + vec[i][j].value * FYP;
  }
  return temp;
}

Field::vectorField FiniteMatrix::correctFaceVelocityEast(Field::vectorField& interpolatedCellFaceVelocity, Field::vectorField& cellFacePressureGrad,
Field::vectorField& DPXField, FiniteMatrix::finiteMat& APinterpolated, Grid& myGrid)
{
  Field::vectorField temp(interpolatedCellFaceVelocity.size(), vector<Field>(interpolatedCellFaceVelocity[0].size()));

  forAllInternalUCVs(temp)
  {
    double DXPE = myGrid.XC[i+1][j] - myGrid.XC[i][j];
    double sArea = myGrid.Y[i][j] - myGrid.Y[i][j-1];
    double volume = DXPE*sArea;

    temp[i][j].value = interpolatedCellFaceVelocity[i][j].value - (APinterpolated[i][j].value * volume
     * (cellFacePressureGrad[i][j].value - DPXField[i][j].value));
  }

  return temp;
}

Field::vectorField FiniteMatrix::correctFaceVelocityNorth(Field::vectorField& interpolatedCellFaceVelocity, Field::vectorField& cellFacePressureGrad,
Field::vectorField& DPYField, FiniteMatrix::finiteMat& APinterpolated, Grid& myGrid)
{
  Field::vectorField temp(interpolatedCellFaceVelocity.size(), vector<Field>(interpolatedCellFaceVelocity[0].size()));

  forAllInternalVCVs(temp)
  {
    double DYPN = myGrid.YC[i][j+1] - myGrid.YC[i][j];
    double sArea = myGrid.X[i][j] - myGrid.X[i-1][j];
    double volume = DYPN*sArea;

    temp[i][j].value = interpolatedCellFaceVelocity[i][j].value - (APinterpolated[i][j].value * volume
     * (cellFacePressureGrad[i][j].value - DPYField[i][j].value));
  }

  return temp;
}

void FiniteMatrix::correctEastMassFluxes(Field::vectorField& massFE, Field::vectorField& PressureCorr, FiniteMatrix::finiteMat& AEmat)
{
  forAllInternalUCVs(massFE)
  {
    massFE[i][j].value += (AEmat[i][j].value * (PressureCorr[i+1][j].value - PressureCorr[i][j].value));
    }
}

void FiniteMatrix::correctNorthMassFluxes(Field::vectorField& massFN, Field::vectorField& PressureCorr, FiniteMatrix::finiteMat& ANmat)
{
  forAllInternalVCVs(massFN)
  {
    massFN[i][j].value += (ANmat[i][j].value * (PressureCorr[i][j+1].value - PressureCorr[i][j].value));
  }
}
*/

void FiniteMatrix::print2dmat(finiteMat &vec)
{
  for (unsigned int j = vec[0].size() - 1; j >= 0; j--)
  {
    for (unsigned int i = 0; i < vec.size(); i++)
    {
      std::cout << std::setprecision(3) << vec[i][j].value << " ";
    }
    std::cout << std::endl;
  }
}

FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat &lhs, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat result(lhs);

  forAllInternal(lhs)
  {
    result[i][j].value -= rhs[i][j].value;
    result[i][j].awvalue -= rhs[i][j].awvalue;
    result[i][j].aevalue -= rhs[i][j].aevalue;
    result[i][j].asvalue -= rhs[i][j].asvalue;
    result[i][j].anvalue -= rhs[i][j].anvalue;
    result[i][j].svalue -= rhs[i][j].svalue;
  }

  return result;
}

FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat &lhs, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat result(lhs);

  forAllInternal(lhs)
  {
    result[i][j].value *= rhs[i][j].value;
    result[i][j].awvalue *= rhs[i][j].awvalue;
    result[i][j].aevalue *= rhs[i][j].aevalue;
    result[i][j].asvalue *= rhs[i][j].asvalue;
    result[i][j].anvalue *= rhs[i][j].anvalue;
    result[i][j].svalue *= rhs[i][j].svalue;
  }

  return result;
}

FiniteMatrix::finiteMat operator*(const double dblvalue, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat result(rhs);

  forAllInternal(rhs)
  {
    result[i][j].value *= dblvalue;
  }

  return result;
}