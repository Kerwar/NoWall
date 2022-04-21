#include <iostream>
#include <cmath>
#include <unistd.h>
#include <thread>
#include <chrono>

#include "Problem.h"
#include "Instrumentor.h"
#include "mpi.h"

using std::cout;
using std::endl;

string showTime(std::chrono::duration<double> time)
{
  string result = std::to_string(std::chrono::duration_cast<std::chrono::hours>(time).count()) +
                  ":" + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(time).count() % 60) +
                  ":" + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(time).count() % 60);

  return result;
}

int main(int argc, char **argv)
{
  PMPI_Init(&argc, &argv);

  int worldRank, worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  Instrumentor::Get().BeginSession("../Profiling/Profile-"+std::to_string(worldRank)+".json");
  auto startTime = std::chrono::high_resolution_clock::now();
  bool solutionNotGoodEnough = true;
  int n = 0, m = 0;
  double mprevious;

  while (solutionNotGoodEnough)
  {
    Problem *problem;
    
    if(n != 0)
      problem = new Problem(worldSize, n, m, mprevious);
    else 
      problem = new Problem(worldSize);

    problem->setUpProblem();

    problem->initializeVariables();
    string filename = "../Results/";
    problem->writeSolution(filename, 0);

    auto fromStart = std::chrono::high_resolution_clock::now();
    auto fromPrev = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationFromStart = fromStart - startTime;
    std::chrono::duration<double> durationFromPrev = fromPrev - startTime;

    if (worldRank == 0)
      cout << "Time from start: " << showTime(durationFromStart) << " Time of last step: " << showTime(durationFromPrev) << " Iteration number: 0 Error: 1" << " M : " << problem->m << endl;;

    double error = 1;
    for (int i = 1; i < problem->maxit; i++)
    {
      error = problem->mainIter(i);


      if (problem->isErrorSmallEnough(error) || i % problem->iShow == 0 || i == problem->maxit - 1)
      {
        fromPrev = std::chrono::high_resolution_clock::now();
        durationFromPrev = fromPrev - fromStart;
        fromStart = std::chrono::high_resolution_clock::now();
        durationFromStart = fromStart - startTime;
        if (worldRank == 0)
          cout << "Time from start: " << showTime(durationFromStart) << " Time of last step: " << showTime(durationFromPrev) << " Iteration number: " << i << " Error: " << error << " M : " << problem->m << endl;
        problem->writeSolution(filename, i);
      if (problem->isErrorSmallEnough(error))
        break;
        // problem.increaseInnerIter(1);
      }
    }
    solutionNotGoodEnough = problem->isSolutionNotGoodEnough();
      // solutionNotGoodEnough = false;//true;

    if (solutionNotGoodEnough)
    {
      problem->retrieveNandM(n, m, mprevious);
      string previous = "";
      problem->writeSolution(previous, -1);
      if (worldRank == 0)
      cout << "Refining the mesh to " << n << " " << m << endl;
    }
    else
    if (worldRank == 0)
	    cout << "We have Computed Everything!" << endl;

    delete problem;
  }
  Instrumentor::Get().EndSession();
  PMPI_Finalize();
  return 0;
}
