#ifndef _SYSTEMOFEQUATIONS_HPP_
#define _SYSTEMOFEQUATIONS_HPP_

#include "equation.hpp"

struct SystemOfEquations {
  Equation *T, *F, *Z;
  ~SystemOfEquations(){
    delete T;
    delete F;
    delete Z;
  }
};

struct Variables {
  Variables(int size_x, int size_y)
      : U(size_x, size_y),
        V(size_x, size_y),
        T(size_x, size_y),
        F(size_x, size_y),
        Z(size_x, size_y) {}

  Field U, V;
  Field T, F, Z;
};

#endif
