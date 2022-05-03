#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../inc/Grid.hpp"
#include "../inc/Field.hpp"
#include <iostream>

TEST_CASE("Factors of weight of points", "[field]")
{
  int N = 10, M = 10;
  double xMin = -10, xMax = 10.0, yMin = 0.0, yMax = 1.0;
  double viscX = 0.0, viscY = 0.0;

  Grid myGrid(N, M, xMin, xMax, yMin, yMax);

  Field fieldOper;
  Field::vectorField field(N, Field::vec1dfield(M));
  fieldOper.getGridInfoPassed(field, myGrid, viscX, viscY);

  forAllInternalUCVs(field)
  {
    INFO("i: " << i << " j: " << j);
    REQUIRE_THAT(field[i][j].FXE,
                 Catch::Matchers::WithinRel(0.5, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
    REQUIRE_THAT(field[i][j].FXP,
                 Catch::Matchers::WithinRel(0.5, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
  }
  forAllInternalVCVs(field)
  {
    INFO("i: " << i << " j: " << j);
    REQUIRE_THAT(field[i][j].FYN,
                 Catch::Matchers::WithinRel(0.5, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
    REQUIRE_THAT(field[i][j].FYP,
                 Catch::Matchers::WithinRel(0.5, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
  }

  forWestBoundary(field)
  {
    INFO("i: " << i << " j: " << j);
    REQUIRE_THAT(field[i][j].FXE,
                 Catch::Matchers::WithinRel(0.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
    REQUIRE_THAT(field[i][j].FXP,
                 Catch::Matchers::WithinRel(1.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
  }
  forEastBoundary(field)
  {
    INFO("i: " << i << " j: " << j);
    REQUIRE_THAT(field[i][j].FXE,
                 Catch::Matchers::WithinRel(0.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
    REQUIRE_THAT(field[i][j].FXP,
                 Catch::Matchers::WithinRel(1.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
  }
  forSouthBoundary(field)
  {
    INFO("i: " << i << " j: " << j);
    REQUIRE_THAT(field[i][j].FYN,
                 Catch::Matchers::WithinRel(0.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
    REQUIRE_THAT(field[i][j].FYP,
                 Catch::Matchers::WithinRel(1.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
  }
  forNorthBoundary(field)
  {
    INFO("i: " << i << " j: " << j);
    REQUIRE_THAT(field[i][j].FYN,
                 Catch::Matchers::WithinRel(0.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
    REQUIRE_THAT(field[i][j].FYP,
                 Catch::Matchers::WithinRel(1.0, 0.001) || Catch::Matchers::WithinAbs(0.5, 0.000001));
  }
};

TEST_CASE("Initialize Internal Field", "[field]")
{
  int N = 100, M = 50;

  Field fieldOper;
  Field::vectorField field(N, Field::vec1dfield(M));

  fieldOper.initializeInternalField(field, 1.0);

  forAllInternal(field)
  {
    REQUIRE(field[i][j].value == 1.0);
  }
};

TEST_CASE("Extrapolate Condition", "[field]")
{
  int N = 100, M = 50;
  double xMin = -10, xMax = 10.0, yMin = 0.0, yMax = 1.0;
  double viscX = 0.0, viscY = 0.0;

  Grid myGrid(N, M, xMin, xMax, yMin, yMax);

  Field fieldOper;
  Field::vectorField field(N, Field::vec1dfield(M));
  fieldOper.getGridInfoPassed(field, myGrid, viscX, viscY);

  fieldOper.initializeInternalField(field, 1.0);
  fieldOper.linearExtrapolateCondition(field, Field::west);

  fieldOper.initializeInternalField(field, 2.0);
  fieldOper.linearExtrapolateCondition(field, Field::east);

  fieldOper.initializeInternalField(field, 3.0);
  fieldOper.linearExtrapolateCondition(field, Field::south);

  fieldOper.initializeInternalField(field, 4.0);
  fieldOper.linearExtrapolateCondition(field, Field::north);

  fieldOper.initializeInternalField(field, 0.0);
  SECTION("West Boundary")
  {
    for (unsigned int j = 1; j < field[0].size() - 1; j++)
    {
      for (unsigned int i = 0; i < 1; i++)
      {
        INFO("i: " << i << " j: " << j);
        double expected = 1.0;
        REQUIRE_THAT(field[i][j].value,
                     Catch::Matchers::WithinRel(expected, 0.001) || Catch::Matchers::WithinAbs(expected, 0.000001));
      }
    }
  }
  SECTION("East Boundary")
  {
  for(unsigned int j = 1; j <field[0].size()-1; j++) 
  {
    for(unsigned int i = field.size()-1; i < field.size(); i++) 
    {
      INFO("i: " << i << " j: " << j);
      double expected = 2.0;
      REQUIRE_THAT(field[i][j].value,
                   Catch::Matchers::WithinRel(expected, 0.001) || Catch::Matchers::WithinAbs(expected, 0.000001));
    }
  }
  }

  SECTION("South Boundary")
  {
    forSouthBoundary(field)
    {
      INFO("i: " << i << " j: " << j);
      double expected = 3.0;
      REQUIRE_THAT(field[i][j].value,
                   Catch::Matchers::WithinRel(expected, 0.001) || Catch::Matchers::WithinAbs(expected, 0.000001));
    }
  }

  SECTION("North Boundary")
  {
    forNorthBoundary(field)
    {
      INFO("i: " << i << " j: " << j);
      double expected = 4.0;
      REQUIRE_THAT(field[i][j].value,
                   Catch::Matchers::WithinRel(expected, 0.001) || Catch::Matchers::WithinAbs(expected, 0.000001));
    }
  }
};