#ifndef FORALLOPERATORS_H
#define FORALLOPERATORS_H
#pragma once
#include <vector>
#include "Instrumentor.h"

using std::vector;

template<class T>
class ForAllOperators
{
    public:
        ForAllOperators<T>();
        virtual ~ForAllOperators<T>();

    vector<T> temp1dvector;
    vector<vector<T>> tempvector;
};

#define forAll(tempvector) \
for(int i = 0; i <tempvector.size(); i++) \
  for(int j = 0; j <tempvector[i].size(); j++)

#define forAllInternal(tempvector) \
for(int i = 1; i <tempvector.size()-1; i++) \
  for(int j = 1; j <tempvector[i].size()-1; j++)

#define forAllInternalUCVs(tempvector) \
for(int i = 1; i <tempvector.size()-2; i++) \
  for(int j = 1; j <tempvector[i].size()-1; j++)

#define forAllInternalVCVs(tempvector) \
for(int i = 1; i <tempvector.size()-1; i++) \
  for(int j = 1; j <tempvector[i].size()-2; j++)

#define forWestBoundary(tempvector) \
for(int j = 0; j <tempvector[0].size(); j++) \
  for(int i = 0; i < 1; i++) 
  
#define forEastBoundary(tempvector) \
for(int j = 0; j <tempvector[0].size(); j++) \
  for(int i = tempvector.size()-1; i < tempvector.size(); i++) 

#define forSouthBoundary(tempvector) \
for(int i = 1; i <tempvector.size()-1; i++) \
  for(int j = 0; j <1; j++)

#define forNorthBoundary(tempvector) \
for(int i = 1; i <tempvector.size()-1; i++) \
  for(int j = tempvector[i].size()-1; j < tempvector[i].size(); j++)

#endif // FORALLOPERATORS_H
