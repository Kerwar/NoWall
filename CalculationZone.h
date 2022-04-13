#ifndef CALCULATIONZONE_H
#define CALCULATIONZONE_H

#include <iostream>
#include <vector>

using std::vector;

class CalculationZone
{
public:

private:
  int NX_Zone, NY_Zone;
  int nProcsInMyZone;

  vector<int> iStrAllProc, iEndAllProc;
  vector<int> jStrAllProc, jEndAllProc;
};

#endif