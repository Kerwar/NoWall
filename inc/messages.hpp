#ifndef __MESSAGES_HPP__
#define __MESSAGES_HPP__

#include <iostream> 
#include <chrono>
#include <string>
#include <iomanip>

using std::string;


string showTime(std::chrono::duration<double> time);

void PrintCurrentStep(const std::chrono::duration<double> &tStart,
                      const std::chrono::duration<double> &tIter,
                      const int &iter, const double &error, const double &m);
#endif
