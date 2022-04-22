#ifndef FORALLOPERATORS_H
#define FORALLOPERATORS_H
#pragma once
#include <vector>
#include "Instrumentor.h"

using std::vector;

template <class T>
class ForAllOperators
{
public:
  ForAllOperators<T>();
  virtual ~ForAllOperators<T>();

  vector<T> temp1dvector;
  vector<vector<T>> tempvector;
};

#define forAll(tempvector)                    \
  for (int i = 0; i < tempvector.size(); i++) \
    for (int j = 0; j < tempvector[i].size(); j++)

#define forAllInternal(tempvector)                \
  for (unsigned int i = 1; i < tempvector.size() - 1; i++) \
    for (unsigned int j = 1; j < tempvector[i].size() - 1; j++)

#define forAllInternalUCVs(tempvector)            \
  for (int i = 1; i < tempvector.size() - 2; i++) \
    for (int j = 1; j < tempvector[i].size() - 1; j++)

#define forAllInternalVCVs(tempvector)            \
  for (int i = 1; i < tempvector.size() - 1; i++) \
    for (int j = 1; j < tempvector[i].size() - 2; j++)

#define forWestBoundary(tempvector)              \
  for (int j = 0; j < tempvector[0].size(); j++) \
    for (int i = 0; i < 1; i++)

#define forEastBoundary(tempvector)              \
  for (int j = 0; j < tempvector[0].size(); j++) \
    for (int i = tempvector.size() - 1; i < tempvector.size(); i++)

#define forSouthBoundary(tempvector)              \
  for (int i = 1; i < tempvector.size() - 1; i++) \
    for (int j = 0; j < 1; j++)

#define forNorthBoundary(tempvector)              \
  for (int i = 1; i < tempvector.size() - 1; i++) \
    for (int j = tempvector[i].size() - 1; j < tempvector[i].size(); j++)

#define forAllN(NI, NJ)        \
  for (int j = 0; j < NJ; j++) \
    for (int i = 0; i < NI; i++)

#define forAllInterior(NI, NJ) \
  for (int j = 1; j < NJ - 1; j++) \
    for (int i = 1; i < NI - 1; i++)

#define forAllInteriorUCVs(NI, NJ) \
  for (int j = 1; j < NJ - 1; j++) \
    for (int i = 1; i < NI - 2; i++)

#define forAllInteriorVCVs(NI, NJ) \
  for (int j = 1; j < NJ - 2; j++) \
    for (int i = 1; i < NI - 1; i++)

#define forWBoundary(NI, NJ)  \
    for (int j = 0; j < NJ; j++) \
  for (int i = 0; i < 1; i++) 

#define forEBoundary(NI, NJ)   \
  for (int j = 0; j < NJ; j++) \
    for (int i = NI - 1; i < NI; i++)

#define forSBoundary(NI, NJ)  \
  for (int j = 0; j < 1; j++) \
    for (int i = 1; i < NI - 1; i++)

#define forNBoundary(NI, NJ)        \
  for (int j = NJ - 1; j < NJ; j++) \
    for (int i = 1; i < NI - 1; i++)

#endif // FORALLOPERATORS_H
