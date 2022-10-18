#include <unistd.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>

#include "communicator.hpp"
#include "environment.hpp"
#include "instrumentor.hpp"
#include "mpi.h"
#include "messages.hpp"
#include "parameters.hpp"
#include "problem.hpp"

using std::cout;
using std::endl;

using parameters::NPROCS;

int main(int argc, char **argv) {
  Environment env(argc, argv);
  Communicator world(MPI_COMM_WORLD);
  
  Instrumentor::Get().BeginSession("../Profiling/Profile-" +
                                   std::to_string(world.rank()) + ".json");

  auto startTime = std::chrono::high_resolution_clock::now();
  if (world.size() == NPROCS) {
    double mprevious = 0.8;
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
      cout << " " << std::setfill('_') << std::setw(130) << "\n"
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
          PrintCurrentStep(durationFromStart, durationFromPrev, i, error,
                           problem.q);
        if (!problem.isErrorSmallEnough(error))
          problem.writeSolution(filename, i);
        else break;
      }
    }

    problem.writeSolution(filename, -1);

    if (world.rank() == 0) cout << "We have Computed Everything!" << endl;
    Instrumentor::Get().EndSession();
  }
  return 0;
}
