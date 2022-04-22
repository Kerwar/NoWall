#include "Equation.h"

Equation::Equation(const FiniteMatrix::finiteMat &Fmatrix) : APinitial(Fmatrix),
                                                             value(0.0), Residual(0.0), RSM(0.0), RESOR(4), URF(0.8), // 0.8),
                                                             EqnName("No Name"), SOR(0.2),
                                                             AP(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AW(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AE(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AS(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AN(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             rAP(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             APU(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             APReaction(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceInitial(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceB(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceFinal(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceRelaxed(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             DT(10e-6),
                                                             UE(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             UN(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             LW(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             LS(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             LPR(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             RES(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             NI(AP.size()),
                                                             NJ(AP[0].size()),
                                                             NIM(NI - 1),
                                                             NJM(NJ - 1),
                                                             Literations(5)
{
  PROFILE_FUNCTION();

  forAllInterior(NI, NJ)
  {
    AW[i][j].value = APinitial[i][j].awvalue;
    AE[i][j].value = APinitial[i][j].aevalue;
    AS[i][j].value = APinitial[i][j].asvalue;
    AN[i][j].value = APinitial[i][j].anvalue;
    APReaction[i][j].value = APinitial[i][j].apvalue;
    AP[i][j].value = 0.0;
    sourceInitial[i][j].value = APinitial[i][j].svalue;
  }
}

Equation::~Equation()
{
}

void Equation::assembleEquation()
{
  PROFILE_FUNCTION();

  forAllInternal(AP)
  {
    AP[i][j].value = -AW[i][j].value - AE[i][j].value - AS[i][j].value - AN[i][j].value + APU[i][j].value + APReaction[i][j].value; // APU is what comes from BD
    // Source Initial does not contain boundaries
    sourceFinal[i][j].value = sourceInitial[i][j].value + sourceB[i][j].value;
  }
}

void Equation::relax(Field &vec)
{
  PROFILE_FUNCTION();
  if (URF != 1.0)
    forAllInternal(AP)
    {
      AP[i][j].value = (1.0 / URF) * AP[i][j].value;
      sourceRelaxed[i][j].value = sourceFinal[i][j].value + (1.0 - URF) * (AP[i][j].value * vec.value[i + j * NI]);

      rAP[i][j].value = 1.0 / AP[i][j].value;
    }
}

void Equation::noWallShearXBoundaryConditions(Field &vec, int start, int end, Field::Direction side) // Only U Velocity
{
  PROFILE_FUNCTION();

  int j;
  switch (side)
  {
  case Field::south:
    j = 1;
    for (int i = start; i < end; i++)
    {
      APU[i][j].value -= AS[i][j].value;
      sourceB[i][j].value -= AS[i][j].value * vec.value[i + j * NI];
      AS[i][j].value = 0;
    }
    break;
  case Field::north:
    j = NJ - 2;
    for (int i = start; i < end; i++)
    {
      APU[i][j].value -= AN[i][j].value;
      sourceB[i][j].value -= AN[i][j].value * vec.value[i + j * NI];
      AN[i][j].value = 0;
    }
    break;
  default:
    std::cout << "NoWallShearXBoundary received a bad direction" << std::endl;
  }
}

void Equation::noWallShearYBoundaryConditions(
    Field &vec, int start, int end, Field::Direction side) // Only V velocity
{
  PROFILE_FUNCTION();
  int i;
  switch (side)
  {
  case Field::west:

    i = 1;

    for (int j = start; j < end; j++)
    {
      APU[i][j].value -= AW[i][j].value;
      sourceB[i][j].value -= AW[i][j].value * vec.value[i + j * NI];
      AW[i][j].value = 0;
    }
    break;

  case Field::east:

    i = NI - 2;

    for (int j = start; j < end; j++)
    {
      APU[i][j].value -= AE[i][j].value;
      sourceB[i][j].value -= AE[i][j].value * vec.value[i + j * NI];
      AE[i][j].value = 0;
    }
    break;
  default:
    std::cout << "NoWallShearYBoundary received a bad direction" << std::endl;
  }
}

double Equation::solve(Field &phi, FiniteMatrix::finiteMat &sourceinput, double &alpha, int &niter, int &iterations, int iterChange)
{
  // SOlver for the SIP(Strongly Implicit Procedure) method
  // Iterative methods: A = M - N -> (M * phi(n) = N * phi(n-1) + S) if phi(n)=phi(n-1) there is no residual
  // In this case we want M = LU with L U only non zero when A is no zero
  // [\\ \        ]             [\           ]                 [\\ \        ]
  // [\\\ \       ]             [\\          ]                 [ \\ \       ]
  // [ \\\ \      ]             [ \\         ]                 [  \\ \      ]
  // [\ \\\ \     ]             [\ \\        ]                 [   \\ \     ]
  // [ \ \\\ \    ]             [ \ \\       ]                 [    \\ \    ]
  // [  \ \\\ \   ]             [  \ \\      ]                 [     \\ \   ]
  // [   \ \\\ \  ]             [   \ \\     ]                 [      \\ \  ]
  // [    \ \\\ \ ]             [    \ \\    ]                 [       \\ \ ]
  // [     \ \\\ \]             [     \ \\   ]                 [        \\ \]
  //        A                          L                              U
  // We can call Lp middle diagonal of L
  // Up we set it to 1
  // The farther diagonal of L is Lw and the close one Ls
  // The farther diagonal of U is Ue and the close one Un
  // When LU is operated appear 2 more diagonal Mse and Mnw
  // Coefficients of Upper and Lower Triangulation Matrices
  PROFILE_FUNCTION();
  if (niter > iterChange)
  {
    iterations = 1;
    return solveGaussSeidel(phi, sourceinput, alpha, iterations);
  }
  iterations = 4;
  // return solveExplicit(phi, sourceinput, alpha, iterations);
  for (int i = 1; i < NI - 1; i++)
  {
    for (int j = 1; j < NJ - 1; j++)
    {
      LW[i][j].value = AW[i][j].value / (1.0 + alpha * UN[i - 1][j].value);
      LS[i][j].value = AS[i][j].value / (1.0 + alpha * UE[i][j - 1].value);

      double P1 = alpha * LW[i][j].value * UN[i - 1][j].value;
      double P2 = alpha * LS[i][j].value * UE[i][j - 1].value;

      LPR[i][j].value = 1.0 / (AP[i][j].value + P1 + P2 - LW[i][j].value * UE[i - 1][j].value - LS[i][j].value * UN[i][j - 1].value);

      UN[i][j].value = (AN[i][j].value - P1) * LPR[i][j].value;
      UE[i][j].value = (AE[i][j].value - P2) * LPR[i][j].value;
    }
  }
  // Calculate Residual and overwriting with an intermidiate vector

  for (int L = 0; L < iterations; L++)
  {
    Residual = 0.0;
    for (int i = 1; i < NI - 1; i++)
    {
      for (int j = 1; j < NJ - 1; j++)
      {
        RES[i][j].value = sourceinput[i][j].value - (AN[i][j].value * phi.value[i + (j + 1) * NI] + AS[i][j].value * phi.value[i + (j - 1) * NI] + AE[i][j].value * phi.value[i + 1 + j * NI] + AW[i][j].value * phi.value[i - 1 + j * NI] + AP[i][j].value * phi.value[i + j * NI]);
        Residual += std::abs(RES[i][j].value);
        RES[i][j].value = LPR[i][j].value * (RES[i][j].value - LS[i][j].value * RES[i][j - 1].value - LW[i][j].value * RES[i - 1][j].value);
      }
    }
    double small = 1e-20;

    if (L == 0)
    {
      RESOR = Residual;
    }
    RSM = Residual / (RESOR + small);

    // if (L == iterations - 1)
    // {
    //   std::cout << EqnName << " Inner it " << L << " and Residual --> " << Residual << " RSM " << std::endl; // RSM << std::endl;
    // }

    // Back Sustitution dn Correction
    for (unsigned int i = NI - 2; i >= 1; i--)
    {
      for (unsigned int j = NJ - 2; j >= 1; j--)
      {
        RES[i][j].value = RES[i][j].value - (UN[i][j].value * RES[i][j + 1].value + UE[i][j].value * RES[i + 1][j].value);

        phi.value[i + j * NI] += RES[i][j].value;
      }
    }
  }
  return Residual;
}

void Equation::SetWallShearTX(Field &vec, int iStr, int iEnd,
                              int Ex1, int Ex2, int myEx1, int myEx2,
                              Field &massFluxEast,
                              Field &massFluxNorth,
                              Field::Direction side)
{
  PROFILE_FUNCTION();

  if (iEnd < Ex1 || iStr >= Ex2)
  {
    noWallShearXBoundaryConditions(vec, 1, vec.NI, side);
  }
  else if (iStr < Ex1 && iEnd < Ex2)
  {
    noWallShearXBoundaryConditions(vec, 1, myEx1, side);
  }
  else if (iStr > Ex1 && iEnd >= Ex2)
  {
    noWallShearXBoundaryConditions(vec, myEx2, vec.NI, side);
  }
  else if (iStr < Ex1 && iEnd >= Ex2)
  {
    noWallShearXBoundaryConditions(vec, 1, myEx1, side);
    noWallShearXBoundaryConditions(vec, myEx2, vec.NI, side);
  }
  SetDirichlet(vec, side, massFluxEast, massFluxNorth);
}

void Equation::SetWallShearX(Field &vec, Field::Direction side)
{
  PROFILE_FUNCTION();
  noWallShearXBoundaryConditions(vec, 1, vec.NI, side);
}

void Equation::SetWallShearY(Field &vec, Field::Direction side)
{
  PROFILE_FUNCTION();
  noWallShearYBoundaryConditions(vec, 1, vec.NJ, side);
}

void Equation::SetDirichlet(Field &vec, Field::Direction side, Field &massFluxEast,
                            Field &massFluxNorth)
{
  PROFILE_FUNCTION();
  switch (side)
  {
  case Field::west:
    for (int j = 1; j < NJ - 1; j++)
    {
      int i = 1;

      sourceB[i][j].value -= vec.value[i - 1 + j * NI] * AW[i][j].value;
      APU[i][j].value -= AW[i][j].value;
      AW[i][j].value = 0;
    }
    break;

  case Field::east:
    for (int j = 1; j < NJ - 1; j++)
    {
      int i = NI - 2;

      sourceB[i][j].value -= vec.value[i + 1 + j * NI] * AE[i][j].value;
      APU[i][j].value -= AE[i][j].value;
      AE[i][j].value = 0;
    }
    break;

  case Field::south:
    for (int i = 1; i < NI - 1; i++)
    {
      int j = 1;

      sourceB[i][j].value -= vec.value[i + (j - 1) * NI] * AS[i][j].value;
      APU[i][j].value -= AS[i][j].value;
      AS[i][j].value = 0;
    }
    break;

  case Field::north:
    for (int i = 1; i < NI - 1; i++)
    {
      int j = vec.NJ - 2;

      sourceB[i][j].value -= vec.value[i + (j + 1) * NI] * AN[i][j].value;
      APU[i][j].value -= AN[i][j].value;
      AN[i][j].value = 0;
    }
    break;
  default:
    std::cout << "SetDirichlet is messup" << std::endl;
    break;
  }
}

double Equation::solveGaussSeidel(Field &phi, FiniteMatrix::finiteMat &sourceinput, double &alpha, int &iterations)
{
  PROFILE_FUNCTION();
  double error;
  URF = 1.0;
  for (int k = 0; k < iterations; k++)
  {
    error = 0;
    forAllInterior(NI, NJ)
    {
      double newvalue = (sourceinput[i][j].value - AW[i][j].value * phi.value[i - 1 + j * NI] - AE[i][j].value * phi.value[i + 1 + j * NI] - AS[i][j].value * phi.value[i + (j - 1) * NI] - AN[i][j].value * phi.value[i + (j + 1) * NI]) / AP[i][j].value;
      error += std::abs(newvalue - phi.value[i + j * NI]);
      phi.value[i + j * NI] = newvalue * (1.0 - alpha) + alpha * phi.value[i + j * NI];
    }
  }
  return error;
}

double Equation::solveExplicit(Field &phi, FiniteMatrix::finiteMat &sourceinput, double &alpha, int &iterations)
{
  PROFILE_FUNCTION();
  double error;

  forAllInterior(NI, NJ)
  {
    double newvalue = DT * (sourceinput[i][j].value -
                            AW[i][j].value * phi.value[i - 1 + j * NI] -
                            AE[i][j].value * phi.value[i + 1 + j * NI] -
                            AS[i][j].value * phi.value[i + (j - 1) * NI] -
                            AN[i][j].value * phi.value[i + (j + 1) * NI] -
                            AP[i][j].value * phi.value[i + j * NI]) +
                      phi.value[i + j * NI];
    error += std::abs(newvalue - phi.value[i + j * NI]);
    phi.value[i + j * NI] = newvalue;
  }
  return error;
}
