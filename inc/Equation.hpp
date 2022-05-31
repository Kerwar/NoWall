#ifndef EQUATION_H
#define EQUATION_H

#include <string>

#include "FiniteMatrix.hpp"

class Equation {
 public:
  explicit Equation(const FiniteMatrix::finiteMat &);
  virtual ~Equation();

  typedef vector<vector<double>> Svector;

  double value;
  double Residual, RSM, RESOR;

  double URF;
  string EqnName;

  FiniteMatrix::finiteMat A;

  double DT;
  void assembleEquation();
  void relax(const Field &);

  void noWallShearXBoundaryConditions(const Field &vec, const int &start,
                                      const int &end,
                                      const Direction &side);
  void noWallShearYBoundaryConditions(const Field &vec, const int &start,
                                      const int &end,
                                      const Direction &side);

  void SetWallShearTX(const Field &vec, const int &iStr, const int &iEnd,
                      const int &Ex1, const int &Ex2, const int &myEx1,
                      const int &myEx2, const Direction &side);
  void SetWallShearX(const Field &vec, const Direction &side);
  void SetWallShearY(const Field &vec, const Direction &side);

  void SetDirichlet(const Field &vec, const Direction &side);
  double solve(Field &phi, const double &alpha, const int &niter,
               int &iterations, const int &iterChange);
  double solveGaussSeidel(Field &phi, const double &alpha,
                          const int &iterations);
  double solveExplicit(Field &phi);

  void inline updateEquation(const FiniteMatrix::finiteMat &fvm) { A = fvm; };

 private:
  Svector UE, UN, LW, LS, LPR, RES;
  int NI, NJ, Literations;
};
#endif