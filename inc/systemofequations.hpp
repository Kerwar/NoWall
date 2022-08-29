#ifndef _SYSTEMOFEQUATIONS_HPP_
#define _SYSTEMOFEQUATIONS_HPP_

#include "equation.hpp"

struct SystemOfEquations {
  Equation *T, *F, *Z;
};

#endif