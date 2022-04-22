#ifndef EQUATION_H
#define EQUATION_H

#include "FiniteMatrix.h"
#include <string>

class Equation{
public:
  Equation(const FiniteMatrix::finiteMat&);
  virtual ~Equation();

  typedef vector<FiniteMatrix> Svector1d;
  typedef FiniteMatrix::finiteMat Svector;

  FiniteMatrix::finiteMat APinitial;

  double value;
  double Residual, RSM, RESOR;

  double URF, rURF; // rUrf removable?
	string EqnName;
	double SOR;

	FiniteMatrix::finiteMat AP;
	FiniteMatrix::finiteMat AW;
	FiniteMatrix::finiteMat AE;
	FiniteMatrix::finiteMat AS;
	FiniteMatrix::finiteMat AN;
	FiniteMatrix::finiteMat rAP;
	FiniteMatrix::finiteMat APU;
  FiniteMatrix::finiteMat APReaction;
	FiniteMatrix::finiteMat sourceInitial;
	FiniteMatrix::finiteMat sourceB;
	FiniteMatrix::finiteMat sourceFinal;
	FiniteMatrix::finiteMat sourceRelaxed;

  double DT;
  void assembleEquation();
  void relax(Field&);

  void noWallShearXBoundaryConditions(Field &vec, int start, int end, Field::Direction side);
  void noWallShearYBoundaryConditions(Field &vec, int start, int end, Field::Direction side);

  void SetWallShearTX(Field &vec, int iStr, int iEnd, int Ex1, int Ex2, int myEx1, int myEx2, Field &massFluxEast, Field &massFluxNorth, Field::Direction side);
  void SetWallShearX(Field &vec, Field::Direction side);
  void SetWallShearY(Field &vec, Field::Direction side);

  void SetDirichlet(Field &vec, Field::Direction side, Field &massFluxEast, Field &massFluxNorth);
  double solve(Field &phi, FiniteMatrix::finiteMat &sourceinput, double &alpha, int &niter, int &iterations, int iterChange);
  double solveGaussSeidel(Field &phi, FiniteMatrix::finiteMat &sourceinput, double &alpha, int &iterations);
  double solveExplicit(Field &phi, FiniteMatrix::finiteMat &sourceinput, double &alpha, int &iterations);
private:

	Svector UE, UN, LW, LS, LPR, RES;
	int NI, NJ, NIM, NJM, Literations;

};
#endif