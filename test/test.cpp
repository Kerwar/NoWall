#define CATCH_CONFIG_MAIN

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_approx.hpp>

#include <mpi.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>

#include "../inc/Paralel.hpp"
#include "../inc/Field.hpp"

TEST_CASE("Processors are not set Up as intended", "[basic]")
{
  int N = 50;
  int M = 50;
  double xMin = -60, xMax = 60.0, xExMin = -20.0, xExMax = 20.0;
  double yMin = 0.0, yMax = 1.5, yWall = 0.5, yChannel = 0.5;

  int rank;

  Paralel paralel;

  paralel.setUpComm(N, M, xMax, yWall, xExMin, xExMax);
  Grid mainGrid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);

  paralel.setUpMesh(N, M, yChannel, yWall, mainGrid.exI1, mainGrid.exI2);

  int USize = paralel.cUSize;
  int DSize = paralel.cDSize;
  int WSize = paralel.wallSize;

  SECTION("Channel and wall distribution")
  {
    CHECK(paralel.worldNProcs == USize + DSize + WSize);
    CHECK(USize >= DSize);
  }

  SECTION("UP CHANNEL NEIGHBOURS", "[UC]")
  {
    if (paralel.loc == Paralel::up)
    {
      for (int i = 0; i < paralel.nProcs; i++)
      {
        // Left Neighbour
        int myLeftWannaBe;

        if (paralel.myLeft != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myLeft, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myLeftWannaBe);
        }
        else
          myLeftWannaBe = MPI_PROC_NULL;

        INFO("Left neightbour in the up channel \n");
        CHECK(paralel.myLeft == myLeftWannaBe);

        int myRightWannaBe;

        // Left Neighbour
        if (paralel.myRight != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myRight, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myRightWannaBe);
        }
        else
          myRightWannaBe = MPI_PROC_NULL;

        INFO("Right neighbour in the up channel \n");
        CHECK(paralel.myRight == myRightWannaBe);

        // Bot neighbour
        int myBotWannaBe;

        if (paralel.myBot != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myBot, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myBotWannaBe);
        }
        else
          myBotWannaBe = MPI_PROC_NULL;
        INFO("Bot neighbour in the up channel \n");
        CHECK(paralel.myBot == myBotWannaBe);

        // Top neighbour
        int myTopWannaBe;

        if (paralel.myTop != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myTop, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myTopWannaBe);
        }
        else
          myTopWannaBe = MPI_PROC_NULL;

        INFO("Top neighbour in the up channel \n");
        CHECK(paralel.myTop == myTopWannaBe);
      }
    }
  }

  SECTION("Down CHANNEL NEIGHBOURS", "[DC]")
  {
    if (paralel.loc == Paralel::down)
    {
      for (int i = 0; i < paralel.nProcs; i++)
      {
        // Left Neighbour
        int myLeftWannaBe;

        if (paralel.myLeft != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myLeft, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myLeftWannaBe);
        }
        else
          myLeftWannaBe = MPI_PROC_NULL;

        INFO("Left neightbour in the up channel \n");
        CHECK(paralel.myLeft == myLeftWannaBe);

        int myRightWannaBe;

        // Left Neighbour
        if (paralel.myRight != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myRight, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myRightWannaBe);
        }
        else
          myRightWannaBe = MPI_PROC_NULL;

        INFO("Right neighbour in the up channel \n");
        CHECK(paralel.myRight == myRightWannaBe);

        // Bot neighbour
        int myBotWannaBe;

        if (paralel.myBot != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myBot, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myBotWannaBe);
        }
        else
          myBotWannaBe = MPI_PROC_NULL;
        INFO("Bot neighbour in the up channel \n");
        CHECK(paralel.myBot == myBotWannaBe);

        // Top neighbour
        int myTopWannaBe;

        if (paralel.myTop != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myTop, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myTopWannaBe);
        }
        else
          myTopWannaBe = MPI_PROC_NULL;

        INFO("Top neighbour in the up channel \n");
        CHECK(paralel.myTop == myTopWannaBe);
      }
    }
  }

  SECTION("Wall NEIGHBOURS", "[W]")
  {
    if (paralel.loc == Paralel::wall)
    {
      for (int i = 0; i < paralel.nProcs; i++)
      {
        // Left Neighbour
        int myLeftWannaBe;

        if (paralel.myLeft != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myLeft, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myLeftWannaBe);
        }
        else
          myLeftWannaBe = MPI_PROC_NULL;

        INFO("Left neightbour in the up channel \n");
        CHECK(paralel.myLeft == myLeftWannaBe);

        int myRightWannaBe;

        // Left Neighbour
        if (paralel.myRight != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myRight, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myRightWannaBe);
        }
        else
          myRightWannaBe = MPI_PROC_NULL;

        INFO("Right neighbour in the up channel \n");
        CHECK(paralel.myRight == myRightWannaBe);

        // Bot neighbour
        int myBotWannaBe;

        if (paralel.myBot != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myBot, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myBotWannaBe);
        }
        else
          myBotWannaBe = MPI_PROC_NULL;
        INFO("Bot neighbour in the up channel \n");
        CHECK(paralel.myBot == myBotWannaBe);

        // Top neighbour
        int myTopWannaBe;

        if (paralel.myTop != MPI_PROC_NULL)
        {
          int myCoords[2];
          MPI_Cart_coords(paralel.myComm, paralel.myTop, 2, myCoords);
          MPI_Cart_rank(paralel.myComm, myCoords, &myTopWannaBe);
        }
        else
          myTopWannaBe = MPI_PROC_NULL;

        INFO("Top neighbour in the up channel \n");
        CHECK(paralel.myTop == myTopWannaBe);
      }
    }
  }
  paralel.freeComm();
};

TEST_CASE("ID is not covering all the processors", "[basic]")
{
  int N = GENERATE(50, 200);
  int M = GENERATE(50, 100);
  double xMin = -60, xMax = 60.0, xExMin = -20.0, xExMax = 20.0;
  double yMin = 0.0, yMax = 1.5, yWall = 0.5, yChannel = 0.5;

  int rank;

  Paralel paralel;

  paralel.setUpComm(N, M, xMax, yWall, xExMin, xExMax);
  Grid mainGrid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);

  paralel.setUpMesh(N, M, yChannel, yWall, mainGrid.exI1, mainGrid.exI2);

  int rowS = paralel.nProcsInRow;
  int colS = paralel.nProcsInCol;

  vector<int> ids(paralel.nProcs, -1);

  for (int i = 0; i < rowS; i++)
  {
    for (int j = 0; j < colS; j++)
    {
      int id = paralel.id(i, j);
      ids[id] = id;
    }
  }

  for (int i = 0; i < ids.size(); i++)
  {
    CHECK(ids[i] >= 0);
  }
  paralel.freeComm();
};

TEST_CASE("Points are not well distributed", "[basic]")
{
  int N = GENERATE(50, 200);
  int M = GENERATE(50, 100);
  double xMin = -60, xMax = 60.0, xExMin = -20.0, xExMax = 20.0;
  double yMin = 0.0, yMax = 1.5, yWall = 0.5, yChannel = 0.5;

  int rank;

  Paralel paralel;

  paralel.setUpComm(N, M, xMax, yWall, xExMin, xExMax);
  Grid mainGrid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);

  paralel.setUpMesh(N, M, yChannel, yWall, mainGrid.exI1, mainGrid.exI2);

  int xProcs = paralel.nProcsInRow;
  int yProcs = paralel.nProcsInCol;

  int NYChannel = floor(M * (2 * yChannel / (2.0 * yChannel + yWall)) / 2.0);
  int NYWall = M - 2 * NYChannel;

  vector<int> iStart = paralel.iStrAllProc;
  vector<int> iEnd = paralel.iEndAllProc;

  vector<int> jStart = paralel.jStrAllProc;
  vector<int> jEnd = paralel.jEndAllProc;

  for (int i = 0; i < paralel.nProcs; i++)
  {
    CHECK(iStart[i] == paralel.iStrAllProc[i]);
    CHECK(iEnd[i] == paralel.iEndAllProc[i]);
    CHECK(jStart[i] == paralel.jStrAllProc[i]);
    CHECK(jEnd[i] == paralel.jEndAllProc[i]);
  }

  SECTION("Points in the Top channel", "[UC]")
  {
    if (paralel.loc == Paralel::up)
    {
      int NY0 = NYChannel + NYWall;
      int NYMax = M;

      CHECK(iEnd[paralel.myProc] - iStart[paralel.myProc] == paralel.myNx - 1);
      CHECK(jEnd[paralel.myProc] - jStart[paralel.myProc] == paralel.myNy - 1);

      for (int j = 0; j < yProcs; j++)
      {
        int id;

        id = paralel.id(0, j);

        INFO("xPos: " << 0 << "/" << xProcs - 1 << " yPos: " << j << "/" << yProcs - 1
                      << std::right << " N: " << N << " M: " << M << " id " << id);
        CHECK(0 == iStart[id]);

        id = paralel.id(xProcs - 1, j);
        CHECK(N - 1 == iEnd[id]);
      }

      int id = paralel.id(paralel.myRowId, paralel.myColId);
      int rightId = paralel.myRight;

      if (paralel.myRight == MPI_PROC_NULL)
      {
        CHECK(N - 1 == iEnd[id]);
      }
      else
      {
        INFO("My proc " << paralel.myProc << "/" << paralel.nProcs - 1
                        << " My Right:" << paralel.myRight << std::right << " Goes from " << iStart[id] << " to " << iEnd[id] << " Right goes from " << iStart[rightId] << " to " << iEnd[rightId] << " of a X size of " << N << " points and " << paralel.myRowId << "/" << paralel.nProcsInRow - 1 << " procs.");
        CHECK(iStart[rightId] - 1 == iEnd[id]);
      }

      for (int i = 0; i < xProcs; i++)
      {
        id = paralel.id(i, 0);

        INFO("xPos: " << xProcs - 1 << "/" << xProcs - 1 << " yPos: " << i << "/" << yProcs - 1
                      << std::right << " N: " << N << " M: " << M << " id " << id);
        CHECK(NY0 == jStart[id]);

        id = paralel.id(i, yProcs - 1);
        CHECK(NYMax - 1 == jEnd[id]);
      }

      id = paralel.id(paralel.myRowId, paralel.myColId);
      int topId = paralel.myTop;

      if (topId == MPI_PROC_NULL)
      {
        CHECK(M - 1 == jEnd[id]);
      }
      else
      {
        CHECK(jStart[topId] - 1 == jEnd[id]);
      }
    }
  }
  SECTION("Points in the Bot channel", "[DC]")
  {
    if (paralel.loc == Paralel::down)
    {
      int NY0 = 0;
      int NYMax = NYChannel;

      CHECK(iEnd[paralel.myProc] - iStart[paralel.myProc] == paralel.myNx - 1);
      CHECK(jEnd[paralel.myProc] - jStart[paralel.myProc] == paralel.myNy - 1);

      for (int j = 0; j < yProcs; j++)
      {
        int id;

        id = paralel.id(0, j);

        INFO("xPos: " << 0 << "/" << xProcs - 1 << " yPos: " << j << "/" << yProcs - 1
                      << std::right << " N: " << N << " M: " << M << " id " << id);
        CHECK(0 == iStart[id]);

        id = paralel.id(xProcs - 1, j);
        CHECK(N - 1 == iEnd[id]);
      }

      for (int i = 1; i < xProcs - 1; i++)
      {
        for (int j = 1; j < yProcs - 1; j++)
        {
          int id = paralel.id(i, j);
          int rightId = paralel.myRight;

          CHECK(iStart[rightId] - 1 == iEnd[id]);
        }
      }

      for (int i = 0; i < xProcs; i++)
      {
        int id = paralel.id(i, 0);

        INFO("xPos: " << xProcs - 1 << "/" << xProcs - 1 << " yPos: " << i << "/" << yProcs - 1
                      << std::right << " N: " << N << " M: " << M << " id " << id);
        CHECK(NY0 == jStart[id]);

        id = paralel.id(i, yProcs - 1);
        CHECK(NYMax - 1 == jEnd[id]);
      }

      for (int i = 1; i < xProcs - 1; i++)
      {
        for (int j = 1; j < yProcs - 1; j++)
        {
          int id = paralel.id(i, j);
          int topId = paralel.myTop;

          CHECK(jStart[topId] - 1 == jEnd[id]);
        }
      }
    }
  }
  SECTION("Points in the Wall", "[W]")
  {

    if (paralel.loc == Paralel::wall)
    {
      int NY0 = NYChannel;
      int NYMax = NYChannel + NYWall;

      CHECK(iEnd[paralel.myProc] - iStart[paralel.myProc] == paralel.myNx - 1);
      CHECK(jEnd[paralel.myProc] - jStart[paralel.myProc] == paralel.myNy - 1);

      double DX = (2.0 * xMax) / N;
      int NXWall = 0;
      while (-xMax + NXWall * DX < xExMin)
      {
        NXWall++;
      }

      NXWall = N - 2 * NXWall;

      CHECK(iEnd[paralel.myProc] - iStart[paralel.myProc] == paralel.myNx - 1);
      CHECK(jEnd[paralel.myProc] - jStart[paralel.myProc] == paralel.myNy - 1);

      for (int j = 0; j < yProcs; j++)
      {
        int id;
        int NtoWall = (N - NXWall) / 2;

        id = paralel.id(0, j);

        INFO("xPos: " << 0 << "/" << xProcs - 1 << " yPos: " << j << "/" << yProcs - 1
                      << std::right << " N: " << N << " M: " << M << " id " << id);
        CHECK(NtoWall == iStart[id]);

        id = paralel.id(xProcs - 1, j);
        CHECK(NXWall + NtoWall - 1 == iEnd[id]);
      }

      for (int i = 1; i < xProcs - 1; i++)
      {
        for (int j = 1; j < yProcs - 1; j++)
        {
          int id = paralel.id(i, j);
          int rightId = paralel.myRight;

          CHECK(iStart[rightId] - 1 == iEnd[id]);
        }
      }

      for (int i = 0; i < xProcs; i++)
      {
        int id = paralel.id(i, 0);

        INFO("xPos: " << xProcs - 1 << "/" << xProcs - 1 << " yPos: " << i << "/" << yProcs - 1
                      << std::right << " N: " << N << " M: " << M << " id " << id);
        CHECK(NY0 == jStart[id]);

        id = paralel.id(i, yProcs - 1);
        CHECK(NYMax - 1 == jEnd[id]);
      }

      for (int i = 1; i < xProcs - 1; i++)
      {
        for (int j = 1; j < yProcs - 1; j++)
        {
          int id = paralel.id(i, j);
          int topId = paralel.myTop;

          CHECK(jStart[topId] - 1 == jEnd[id]);
        }
      }
    }
  }
};

TEST_CASE("Neighbours are not receiving the right info", "[info]"
                                                         "[neigh]")
{
  int N = 100;
  int M = 50;
  double xMin = -60, xMax = 60.0, xExMin = -20.0, xExMax = 20.0;
  double yMin = 0.0, yMax = 1.5, yWall = 0.5, yChannel = 0.5;

  Paralel paralel;

  paralel.setUpComm(N, M, xMax, yWall, xExMin, xExMax);
  Grid mainGrid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);

  paralel.setUpMesh(N, M, yChannel, yWall, mainGrid.exI1, mainGrid.exI2);

  double myRank = paralel.myProc;

    double myXMin = mainGrid.XC[paralel.iStr][0];
  double myXMax = mainGrid.XC[paralel.iEnd + 2][0];
  double myYMin = 0.0;
  double myYMax = yChannel;
    if (paralel.loc == Paralel::wall)
  {
    myXMin = mainGrid.X[paralel.exI1][0];
    myXMax = mainGrid.X[paralel.exI2][0];
    myYMin += yChannel;
    myYMax += yWall;
  }
  else if (paralel.loc == Paralel::up)
  {
    myYMin += yChannel + yWall;
    myYMax += yChannel + yWall;
  }

  Grid myGrid(paralel.myNx, paralel.myNy, myXMin, myXMax, myYMin, myYMax);
  myGrid.SetIEx(xExMin, xExMax);

  Field fieldOper;
  Field::vectorField field(paralel.myNx, Field::vec1dfield(paralel.myNy));

  double viscX = 0.0, viscY = 0.0;
  fieldOper.getGridInfoPassed(field, myGrid, viscX, viscY);
  // std::cout << "SITN " << paralel.myNx << " " << paralel.myNy << "\n";
  fieldOper.initializeInternalField(field, myRank);

  // for(unsigned int i = 0; i < field.size(); i++)
  //   {
  //     for(unsigned int j = 0; j < field[i].size(); j++)
  //   {
  //     std::cout << field[i][j].value  << " ";
  //   }
  //   std::cout<<  std::endl;
  //   }
  fieldOper.linearExtrapolateCondition(field, Field::west);
  fieldOper.linearExtrapolateCondition(field, Field::east);
  fieldOper.linearExtrapolateCondition(field, Field::south);
  fieldOper.linearExtrapolateCondition(field, Field::north);

  // for(unsigned int i = 0; i < field.size(); i++)
  //   {
  //     for(unsigned int j = 0; j < field[i].size(); j++)
  //   {
  //     std::cout << field[i][j].value  << " ";
  //   }
  //   std::cout<<  std::endl;
  //   }
  paralel.SendInfoToNeighbours(field);

  SECTION(std::string("Info Sent Left ") + std::to_string(paralel.loc))
  {
    if (paralel.myRight >= 0)
    {
      SECTION("Right is not NULL")
      {
        for (unsigned int j = 1; j < field[0].size() - 1; j++)
        {
          for (unsigned int i = field.size() - 1; i < field.size(); i++)
          {
            double expected = (myRank + paralel.myRight) / 2.0;
            double value = field[i][j].value;

            CHECK(value == Catch::Approx(expected));
          }
        }
      }
    }
    else
    {
      SECTION("Right is NULL")
      {
        for (unsigned int j = 1; j < field[0].size() - 1; j++)
        {
          for (unsigned int i = field.size() - 1; i < field.size(); i++)
          {
            double expected = myRank;
            double value = field[i][j].value;

            CHECK(value == Catch::Approx(expected));
          }
        }
      }
    }
  }

  SECTION(std::string("Info Sent Right ") + std::to_string(paralel.loc))
  {
    if (paralel.myLeft >= 0)
    {
      SECTION("Left is not NULL")
      {
        for (unsigned int j = 1; j < field[0].size() - 1; j++)
        {
          for (unsigned int i = 0; i < 1; i++)
          {
            double expected = (myRank + paralel.myLeft) / 2.0;
            double value = field[i][j].value;

            CHECK(value == Catch::Approx(expected));
          }
        }
      }
    }
    else
    {
      SECTION("Left is NULL")
      {
        for (unsigned int j = 1; j < field[0].size() - 1; j++)
        {
          for (unsigned int i = 0; i < 1; i++)
          {
            double expected = myRank;
            double value = field[i][j].value;

            CHECK(value == Catch::Approx(expected));
          }
        }
      }
    }
  }

  SECTION(std::string("Info Sent Bot ") + std::to_string(paralel.loc))
  {
    if (paralel.myTop >= 0)
    {
      SECTION("Top is not NULL")
      {
        forNorthBoundary(field)
        {
          double expected = (myRank + paralel.myTop) / 2.0;
          double value = field[i][j].value;

          CHECK(value == Catch::Approx(expected));
        }
      }
    }
    else
    {
      SECTION("Top is NULL")
      {
        forNorthBoundary(field)
        {
          double expected = myRank;
          double value = field[i][j].value;

          CHECK(value == Catch::Approx(expected));
        }
      }
    }
  }

  SECTION(std::string("Info Sent Top ") + std::to_string(paralel.loc))
  {
    if (paralel.myBot >= 0)
    {
      SECTION("Bot is not NULL")
      {
        forSouthBoundary(field)
        {
          double expected = (myRank + paralel.myBot) / 2.0;
          double value = field[i][j].value;

          INFO("xPos: " << i + paralel.iStrAllProc[paralel.myProc] << "/" << paralel.iEndAllProc[paralel.myProc] << " yPos: " << j + paralel.jStrAllProc[paralel.myProc] << "/" << paralel.jEndAllProc[paralel.myProc] << std::right << " Me: " << paralel.myProc << " My bot: " << paralel.myBot);
          CHECK(value == Catch::Approx(expected));
        }
      }
    }
    else
    {
      SECTION("Bot is NULL")
      {
        forSouthBoundary(field)
        {
          double expected = myRank;
          double value = field[i][j].value;

          CHECK(value == Catch::Approx(expected));
        }
      }
    }
  }
  paralel.freeComm();
};

TEST_CASE("Main Proc is not receiving the info properly", "[info]")
{
  int N = 100;
  int M = 50;
  double xMin = -60, xMax = 60.0, xExMin = -20.0, xExMax = 20.0;
  double yMin = 0.0, yMax = 1.5, yWall = 0.5, yChannel = 0.5;

  Paralel paralel;

  paralel.setUpComm(N, M, xMax, yWall, xExMin, xExMax);
  Grid mainGrid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);

  paralel.setUpMesh(N, M, yChannel, yWall, mainGrid.exI1, mainGrid.exI2);

  int myRank = paralel.myProc;

  Field fieldOper;
  Field::vectorField field(paralel.myNx, Field::vec1dfield(paralel.myNy));
  Field::vectorField sol(N, Field::vec1dfield(M));

  // std::cout << "SITN " << paralel.myNx << " " << paralel.myNy << "\n";
  fieldOper.initializeInternalField(field, myRank);

  fieldOper.linearExtrapolateCondition(field, Field::west);
  fieldOper.linearExtrapolateCondition(field, Field::east);
  fieldOper.linearExtrapolateCondition(field, Field::south);
  fieldOper.linearExtrapolateCondition(field, Field::north);

  paralel.SendInfoToCommMainProc(field, sol, paralel.myComm);

  if (myRank == 0)
  {
    vector<int> iStart = paralel.iStrAllProc;
    vector<int> iEnd = paralel.iEndAllProc;

    vector<int> jStart = paralel.jStrAllProc;
    vector<int> jEnd = paralel.jEndAllProc;

    for (unsigned int k = 0; k < paralel.nProcs; k++)
    {
      for (unsigned int i = iStart[k]; i <= iEnd[k]; i++)
      {
        for (unsigned int j = jStart[k]; j <= jEnd[k]; j++)
        {
          if (j == jStart[k] || j == jEnd[k])
            if (i == iStart[k] || i == iEnd[k])
              continue;
          INFO("xPos: " << i << "/" << iEnd[k] << " yPos: " << j << "/" << jEnd[k]);
          CHECK(sol[i][j].value == Catch::Approx(k));
        }
      }
    }
  }
  paralel.freeComm();
};

TEST_CASE("The Main Proc of each zones have trouble communicating", 
                                                                    "[exT]")
{
  // Values that are initialized in the test case and the solution to it
  // Index NN    1.5
  // Index  N    1.0
  // Wall BD      1,936750562
  // Index  O    2.0
  // Index  S    2.5
  //
  //
  // Index NN   2.5
  // Index  N   2.0
  // Wall BD    2,936750562
  // Index  O   3.0
  // Index  S   3.5

  int N = 100;
  int M = 50;
  double xMin = -60, xMax = 60.0, xExMin = -20.0, xExMax = 20.0;
  double yMin = 0.0, yMax = 1.5, yWall = 0.5, yChannel = 0.5;

  Paralel paralel;

  paralel.setUpComm(N, M, xMax, yWall, xExMin, xExMax);
  Grid mainGrid(N, M, xMin, xMax, yMin, yMax);
  mainGrid.SetIEx(xExMin, xExMax);

  paralel.setUpMesh(N, M, yChannel, yWall, mainGrid.exI1, mainGrid.exI2);

  int myRank = paralel.myProc;
  Paralel::Loc loc = paralel.loc;

  Field fieldOper;

  Field::vectorField T(N, Field::vec1dfield(M));
  Field::vec1dfield TWallUp(N), TWallDown(N);

  switch (loc) // Initialize T
  {
  case Paralel::up:
    fieldOper.initializeInternalField(T, 1.0);
    fieldOper.linearExtrapolateCondition(T, Field::south);
    fieldOper.initializeInternalField(T, 1.5);
    break;

  case Paralel::wall:
    fieldOper.initializeInternalField(T, 2.0);
    fieldOper.linearExtrapolateCondition(T, Field::north);
    fieldOper.linearExtrapolateCondition(T, Field::south);
    fieldOper.initializeInternalField(T, 2.5);
    // fieldOper.linearExtrapolateCondition(T, Field::east);
    // fieldOper.linearExtrapolateCondition(T, Field::west);
    break;

  case Paralel::down:
    fieldOper.initializeInternalField(T, 3.0);
    fieldOper.linearExtrapolateCondition(T, Field::north);
    fieldOper.initializeInternalField(T, 3.5);
    break;
  }

  double bK = 0.15, a2 = 0.01;

  if (myRank == 0)
  {
    // std::cout << myRank << " " << paralel.UCMainProc << " " << paralel.WMainProc << " " << paralel.DCMainProc << std::endl;
    paralel.ExchangeWallTemperature(T, TWallUp, TWallDown, bK, a2);

    int exI1 = paralel.exI1, exI2 = paralel.exI2;

    switch (loc)
    {
    case Paralel::up:
      for (int i = exI1 + 1; i < exI2; i++)
      {
        CHECK(T[i][0].value == Catch::Approx(1.936750562));
      }
      break;

    case Paralel::wall:
      for (int i = 1; i < N - 1; i++)
      {
        CHECK(T[i][paralel.NYWall - 1].value == Catch::Approx(1.936750562));
        CHECK(T[i][0].value == Catch::Approx(2.936750562).epsilon(0.01));
      }
      break;

    case Paralel::down:
      for (int i = exI1 + 1; i < exI2; i++)
      {
        CHECK(T[i][paralel.NYChannel - 1].value == Catch::Approx(2.936750562).epsilon(0.01));
      }
      break;
    }
  }
}

/*
TEST_CASE("The information of the Temperature is not distributed right", "[info]"
                                                                         "[exT]")
{
  // Values that are initialized in the test case and the solution to it
  // Index NN    1.5
  // Index  N    1.0
  // Wall BD      1,936750562
  // Index  O    2.0
  // Index  S    2.5
  //
  //
  // Index NN   2.5
  // Index  N   2.0
  // Wall BD    2,936750562
  // Index  O   3.0
  // Index  S   3.5

  int N = 100;
  int M = 50;
  double xMax = 60.0, yWall = 0.5, yChannel = 0.5, xEx = 20.0;

  Paralel paralel;

  paralel.setUpComm(N, M, xMax, yWall, xEx);
  paralel.setUpMesh(N, M, xMax, yWall, xEx);

  int myRank = paralel.myProc;
  Paralel::Loc loc = paralel.loc;

  Field fieldOper;

  int NX = paralel.myNx, NY = paralel.myNy;
  Field::vectorField Tsol(N, Field::vec1dfield(M));
  Field::vec1dfield TWallUp(N), TWallDown(M);
  Field::vectorField T(NX, Field::vec1dfield(NY));

  fieldOper.initializeInternalField(T, 0.0);

  switch (loc) // Initialize Tsol
  {
  case Paralel::up:
    fieldOper.initializeInternalField(Tsol, 1.0);
    fieldOper.linearExtrapolateCondition(Tsol, Field::south);
    fieldOper.initializeInternalField(Tsol, 1.5);
    break;

  case Paralel::wall:
    fieldOper.initializeInternalField(Tsol, 2.0);
    fieldOper.linearExtrapolateCondition(Tsol, Field::north);
    fieldOper.linearExtrapolateCondition(Tsol, Field::south);
    fieldOper.initializeInternalField(Tsol, 2.5);
    break;

  case Paralel::down:
    fieldOper.initializeInternalField(Tsol, 3.0);
    fieldOper.linearExtrapolateCondition(Tsol, Field::north);
    fieldOper.initializeInternalField(Tsol, 3.5);
    break;
  }

  double bK = 0.15, a2 = 0.01;

  if (myRank == 0)
  {
    paralel.ExchangeWallTemperature(Tsol, TWallUp, TWallDown, bK, a2);
  }
  int exI1 = paralel.exI1, exI2 = paralel.exI2;

  paralel.ShareWallTemperatureInfo(TWallUp, TWallDown, T);

  if ((paralel.iStr >= exI1 && paralel.iStr <= exI2) || (paralel.iEnd >= exI1 && paralel.iEnd <= exI2) || (paralel.iStr <= exI1 && paralel.iEnd >= exI2))
  {
    }
}*/