#include "Equation.hpp"

Equation::Equation(const FiniteMatrix::finiteMat &Fmatrix) : value(0.0), Residual(0.0), RSM(0.0), RESOR(4), URF(0.8), // 0.8),
                                                             EqnName("No Name"), SOR(0.2),
                                                             A(Fmatrix),
                                                             DT(10e-6),
                                                             UE(Fmatrix.size(), vector<double>(Fmatrix[0].size(),0)),
                                                             UN(Fmatrix.size(), vector<double>(Fmatrix[0].size(),0)),
                                                             LW(Fmatrix.size(), vector<double>(Fmatrix[0].size(),0)),
                                                             LS(Fmatrix.size(), vector<double>(Fmatrix[0].size(),0)),
                                                             LPR(Fmatrix.size(), vector<double>(Fmatrix[0].size(),0)),
                                                             RES(Fmatrix.size(), vector<double>(Fmatrix[0].size(),0)),
                                                             NI(A.size()),
                                                             NJ(A[0].size()),
                                                             NIM(NI - 1),
                                                             NJM(NJ - 1),
                                                             Literations(5)
{
}

Equation::~Equation()
{
}

void Equation::assembleEquation()
{
  PROFILE_FUNCTION();

  forAllInternal(A)
  {
    A[i][j].ap += -A[i][j].aw - A[i][j].ae - A[i][j].as - A[i][j].an;
    // Source Initial does not contain boundaries
  }
}

void Equation::relax(Field &vec)
{
  PROFILE_FUNCTION();
  if (URF != 1.0)
    forAllInternal(A)
    {
      A[i][j].ap = (1.0 / URF) * A[i][j].ap;
      A[i][j].svalue += (1.0 - URF) * (A[i][j].ap * vec.value[i + j * NI]);
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
      A[i][j].ap -= A[i][j].as;
      A[i][j].svalue -= A[i][j].as * vec.value[i + j * NI];
      A[i][j].as = 0;
    }
    break;
  case Field::north:
    j = NJ - 2;
    for (int i = start; i < end; i++)
    {
      A[i][j].ap -= A[i][j].an;
      A[i][j].svalue -= A[i][j].an * vec.value[i + j * NI];
      A[i][j].an = 0;
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
      A[i][j].ap -= A[i][j].aw;
      A[i][j].svalue -= A[i][j].aw * vec.value[i + j * NI];
      A[i][j].aw = 0;
    }
    break;

  case Field::east:

    i = NI - 2;

    for (int j = start; j < end; j++)
    {
      A[i][j].ap -= A[i][j].ae;
      A[i][j].svalue -= A[i][j].ae * vec.value[i + j * NI];
      A[i][j].ae = 0;
    }
    break;
  default:
    std::cout << "NoWallShearYBoundary received a bad direction" << std::endl;
  }
}

double Equation::solve(Field &phi, double &alpha, int &niter, int &iterations, int iterChange)
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
    return solveGaussSeidel(phi, alpha, iterations);
  }
  iterations = 4;
  // return solveExplicit(phi, alpha, iterations);
  for (int i = 1; i < NI - 1; i++)
  {
    for (int j = 1; j < NJ - 1; j++)
    {
      LW[i][j] = A[i][j].aw / (1.0 + alpha * UN[i - 1][j]);
      LS[i][j] = A[i][j].as / (1.0 + alpha * UE[i][j - 1]);

      double P1 = alpha * LW[i][j] * UN[i - 1][j];
      double P2 = alpha * LS[i][j] * UE[i][j - 1];

      LPR[i][j] = 1.0 / (A[i][j].ap + P1 + P2 - LW[i][j] * UE[i - 1][j] - LS[i][j] * UN[i][j - 1]);

      UN[i][j] = (A[i][j].an - P1) * LPR[i][j];
      UE[i][j] = (A[i][j].ae - P2) * LPR[i][j];
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
        RES[i][j] = A[i][j].svalue - (A[i][j].an * phi.value[i + (j + 1) * NI] + A[i][j].as * phi.value[i + (j - 1) * NI] + A[i][j].ae * phi.value[i + 1 + j * NI] + A[i][j].aw * phi.value[i - 1 + j * NI] + A[i][j].ap * phi.value[i + j * NI]);
        Residual += std::abs(RES[i][j]);
        RES[i][j] = LPR[i][j] * (RES[i][j]- LS[i][j] * RES[i][j - 1]- LW[i][j] * RES[i - 1][j]);
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
        RES[i][j] = RES[i][j] - (UN[i][j] * RES[i][j + 1] + UE[i][j] * RES[i + 1][j]);

        phi.value[i + j * NI] += RES[i][j];
      }
    }
  }
  return Residual;
}

void Equation::SetWallShearTX(Field &vec, int iStr, int iEnd,
                              int Ex1, int Ex2, int myEx1, int myEx2,
                              Field::Direction side)
{
  PROFILE_FUNCTION();

  if (iEnd <= Ex1 || iStr >= Ex2)
    noWallShearXBoundaryConditions(vec, 1, vec.NI, side);
  else if (iStr < Ex1 && iEnd <= Ex2)
    noWallShearXBoundaryConditions(vec, 1, myEx1, side);
  else if (iStr >= Ex1 && iEnd > Ex2)
    noWallShearXBoundaryConditions(vec, myEx2, vec.NI, side);
  else if (iStr < Ex1 && iEnd >= Ex2)
  {
    noWallShearXBoundaryConditions(vec, 1, myEx1, side);
    noWallShearXBoundaryConditions(vec, myEx2, vec.NI, side);
  }

  SetDirichlet(vec, side);
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

void Equation::SetDirichlet(Field &vec, Field::Direction side)
{
  PROFILE_FUNCTION();
  switch (side)
  {
  case Field::west:
    for (int j = 1; j < NJ - 1; j++)
    {
      int i = 1;

      A[i][j].svalue -= vec.value[i - 1 + j * NI] * A[i][j].aw;
      A[i][j].ap -= A[i][j].aw;
      A[i][j].aw = 0;
    }
    break;

  case Field::east:
    for (int j = 1; j < NJ - 1; j++)
    {
      int i = NI - 2;

      A[i][j].svalue -= vec.value[i + 1 + j * NI] * A[i][j].ae;
      A[i][j].ap -= A[i][j].ae;
      A[i][j].ae = 0;
    }
    break;

  case Field::south:
    for (int i = 1; i < NI - 1; i++)
    {
      int j = 1;

      A[i][j].svalue -= vec.value[i + (j - 1) * NI] * A[i][j].as;
      A[i][j].ap -= A[i][j].as;
      A[i][j].as = 0;
    }
    break;

  case Field::north:
    for (int i = 1; i < NI - 1; i++)
    {
      int j = vec.NJ - 2;

      A[i][j].svalue -= vec.value[i + (j + 1) * NI] * A[i][j].an;
      A[i][j].ap -= A[i][j].an;
      A[i][j].an = 0;
    }
    break;
  default:
    std::cout << "SetDirichlet is messup" << std::endl;
    break;
  }
}

double Equation::solveGaussSeidel(Field &phi, double &alpha, int &iterations)
{
  PROFILE_FUNCTION();
  double error;
  URF = 1.0;
  for (int k = 0; k < iterations; k++)
  {
    error = 0;
    forAllInterior(NI, NJ)
    {
      double newvalue = (A[i][j].svalue - A[i][j].aw * phi.value[i - 1 + j * NI] - A[i][j].ae * phi.value[i + 1 + j * NI] - A[i][j].as * phi.value[i + (j - 1) * NI] - A[i][j].an * phi.value[i + (j + 1) * NI]) / A[i][j].ap;
      error += std::abs(newvalue - phi.value[i + j * NI]);
      phi.value[i + j * NI] = newvalue * (1.0 - alpha) + alpha * phi.value[i + j * NI];
    }
  }
  return error;
}

double Equation::solveExplicit(Field &phi)
{
  PROFILE_FUNCTION();
  double error;

  forAllInterior(NI, NJ)
  {
    double newvalue = DT * (A[i][j].svalue -
                            A[i][j].aw * phi.value[i - 1 + j * NI] -
                            A[i][j].ae * phi.value[i + 1 + j * NI] -
                            A[i][j].as * phi.value[i + (j - 1) * NI] -
                            A[i][j].an * phi.value[i + (j + 1) * NI] -
                            A[i][j].ap * phi.value[i + j * NI]) +
                      phi.value[i + j * NI];
    error += std::abs(newvalue - phi.value[i + j * NI]);
    phi.value[i + j * NI] = newvalue;
  }
  return error;
}
