#include <unistd.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>

#include "Instrumentor.hpp"
#include "communicator.hpp"
#include "environment.hpp"
#include "mpi.h"
#include "parameters.hpp"
#include "problem.hpp"

using std::cout;
using std::endl;

using parameters::M;
using parameters::N;
using parameters::NPROCS;

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
  Environment env(argc, argv);
  Communicator world(MPI_COMM_WORLD);

  Instrumentor::Get().BeginSession("../Profiling/Profile-" +
                                   std::to_string(world.rank()) + ".json");

  auto startTime = std::chrono::high_resolution_clock::now();
  if (world.size() == NPROCS) {
    double mprevious = 2;
    Problem problem(mprevious);
    if (argc > 1) problem.xfix = std::atof(argv[1]);

    problem.setUpProblem();

    problem.initializeVariables();

    string filename = "../Results/";
    problem.writeSolution(filename, 0);

    auto fromStart = std::chrono::high_resolution_clock::now();
    auto fromPrev = std::chrono::high_resolution_clock::now();

    auto durationFromStart = fromStart - startTime;
    auto durationFromPrev = fromPrev - startTime;

    if (world.rank() == 0) {
      cout << "The size of the Mesh is " << NTOTAL << "x" << MINPUT << endl;
      cout << " " << std::setfill('_') << std::setw(120) << "\n"
           << std::setfill(' ');
      PrintCurrentStep(durationFromStart, durationFromPrev, 0, 0, problem.q);
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
        if (world.rank() == 0)
          PrintCurrentStep(durationFromStart, durationFromPrev, i, error, problem.q);
        if (!problem.isErrorSmallEnough(error))
          problem.writeSolution(filename, i);
        if (problem.isErrorSmallEnough(error)) break;
      }
    }

    string previous = "";
    problem.writeSolution(filename, -1);

    if (world.rank() == 0) cout << "We have Computed Everything!" << endl;
    Instrumentor::Get().EndSession();
  }

  return 0;
}
