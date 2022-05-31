#include <unistd.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <thread>

#include "Instrumentor.hpp"
#include "Problem.hpp"
#include "mpi.h"

using std::cout;
using std::endl;

string showTime(std::chrono::duration<double> time) {
  string result =
      std::to_string(
          std::chrono::duration_cast<std::chrono::hours>(time).count()) +
      ":" +
      std::to_string(
          std::chrono::duration_cast<std::chrono::minutes>(time).count() % 60) +
      ":" +
      std::to_string(
          std::chrono::duration_cast<std::chrono::seconds>(time).count() % 60);

  return result;
}

void PrintCurrentStep(const std::chrono::duration<double> &tStart,
                      const std::chrono::duration<double> &tIter,
                      const int &iter, const double &error, const double &m) {
  cout << " |Time from start: " << std::setw(7) << showTime(tStart)
       << " |Time of last step: " << std::setw(7) << showTime(tIter)
       << " |Iteration number: " << std::setw(10) << iter 
       << " |Error: " << std::setw(12) << std::setprecision(7) << error 
       << " |M : " << std::setw(8) << std::setprecision(6) << m << "|\n";
       
}

int main(int argc, char **argv) {
  PMPI_Init(&argc, &argv);

  int worldRank, worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  Instrumentor::Get().BeginSession("../Profiling/Profile-" +
                                   std::to_string(worldRank) + ".json");

  auto startTime = std::chrono::high_resolution_clock::now();

  constexpr int n = 1200;
  constexpr int M = 10;
  constexpr int NPROCS = 2;
  constexpr int N = n + (NPROCS / 2 - n % (NPROCS / 2));

  if (worldSize == NPROCS) {
    double mprevious = 2;
    Problem<N, M, NPROCS> problem(mprevious);

    problem.setUpProblem();

    problem.initializeVariables();

    string filename = "../Results/";
    problem.writeSolution(filename, 0);

    auto fromStart = std::chrono::high_resolution_clock::now();
    auto fromPrev = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationFromStart = fromStart - startTime;
    std::chrono::duration<double> durationFromPrev = fromPrev - startTime;

    if (worldRank == 0) {
      cout << "The size of the Mesh is " << N << "x" << M << endl;
      cout << " " << std::setfill('_') << std::setw(120) << "\n" << std::setfill(' ');
      PrintCurrentStep(durationFromStart, durationFromPrev, 0, 0, problem.m);
      ;
    }
    for (int i = 1; i < problem.maxit; i++) {
      double error = 1;
      error = problem.mainIter(i);

      if (problem.isErrorSmallEnough(error) || i % problem.iShow == 0 ||
          i == problem.maxit - 1) {
        fromPrev = std::chrono::high_resolution_clock::now();
        durationFromPrev = fromPrev - fromStart;
        fromStart = std::chrono::high_resolution_clock::now();
        durationFromStart = fromStart - startTime;
        if (worldRank == 0)
          PrintCurrentStep(durationFromStart, durationFromPrev, i, error,
                           problem.m);
        problem.writeSolution(filename, i);
        if (problem.isErrorSmallEnough(error)) break;
      }
    }

    string previous = "";
    problem.writeSolution(previous, -1);

    if (worldRank == 0) cout << "We have Computed Everything!" << endl;
    Instrumentor::Get().EndSession();
  }
  PMPI_Finalize();
  return 0;
}
