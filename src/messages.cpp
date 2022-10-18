#include "messages.hpp"

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
  std::cout << " |Time from start: " << std::scientific <<  std::setw(7) << showTime(tStart)
       << " |Time of last step: " << std::setw(7) << showTime(tIter)
       << " |Iteration number: " << std::setw(10) << iter
       << " |Error: " << std::setw(12) << std::setprecision(9) << error
       << " |M : " << std::setw(10) << std::setprecision(8) << m << "|\n";
}
