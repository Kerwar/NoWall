#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_approx.hpp>

#include "../inc/grid.hpp"
#include "../inc/field.hpp"
#include "../inc/paralel.hpp"
#include <iostream>

TEST_CASE("Initializing Grid", "Grid")
{
  int N = 10;
  int M = 2;
  double xMin = -60, xMax = 60.0, yWall = 0.5, yChannel = 0.5, xExMin = -20, xExMax = 20;
  double yMin = 0, yMax = 1.5;

  Grid mainGrid(9, 2, 0, 10, 0, 1);

  CHECK(mainGrid.N == 9);
  CHECK(mainGrid.M == 2);

  CHECK(mainGrid.NI == 11);
  CHECK(mainGrid.NJ == 4);
}