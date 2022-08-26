#ifndef READDATA
#define READDATA

#include <filesystem>
#include <iostream>
#include <string>

#include "filereader.hpp"

namespace fs = std::filesystem;

int indexofT1(const Field &vec, const int &NI, const int &NJ) {
  int result = 0;

  while (vec.value[result + (NJ - 1) * NI] < 1) result++;

  if (result >= NI) {
    result = 0;

    for (int i = 0; i < NI; i++)
      if (vec[result + (NJ - 1) * NI] < vec[i + (NJ - 1) * NI]) {
        result = i;
      }
    // std::cout << i << " " << vec[result + (NJ-1)*NI] << std::endl;}
  }
  return result;
}

// Sol_NxM-2432x20_Lxa-120x1_Ex-40_q-1.2_m-2_beta-10_LeFxZ-1x0.3_-1.q
void readdata(const std::string &path, std::vector<double> &a,
              std::vector<double> &T_before1, std::vector<double> &T_after1,
              std::vector<double> &XC, std::vector<double> &YC,
              std::vector<double> &q, std::vector<double> &m,
              std::vector<int> &indexes) {
  for (const auto &entry : fs::directory_iterator(path)) {
    std::string filename = entry.path();

    if (filename.find("Sol_NxM-2432x20") != -1 &&
        filename.find("_-1.f") != -1) {
      filename = filename.substr(filename.find("Sol"), filename.length());

      std::string dimensions =
          filename.substr(filename.find("NxM-") + 4,
                          filename.find("_Lxa") - filename.find("NxM-") - 4);

      int N, M;
      N = std::stoi(dimensions.substr(0, dimensions.find("x")));
      M = std::stoi(
          dimensions.substr(dimensions.find("x") + 1, dimensions.length() - 1));

      std::string adimensions =
          filename.substr(filename.find("_Lxa-") + 5,
                          filename.find("_Ex") - filename.find("_Lxa-") - 3);
      double a_value;
      a_value = std::stod(adimensions.substr(adimensions.find("x") + 1,
                                             adimensions.length() - 1));

      std::string qstring =
          filename.substr(filename.find("_q-") + 3,
                          filename.find("_m-") - filename.find("_q-") - 3);
      double q_value = std::stod(qstring);

      std::string mstring =
          filename.substr(filename.find("_m-") + 3,
                          filename.find("_beta-") - filename.find("_m-") - 3);
      double m_value = std::stod(mstring);

      FileReader filerader;

      Field vec(N + 2, M + 2);
      filerader.readField(entry.path(), 2, 1, vec, 0, N + 2);

      indexes.push_back(indexofT1(vec, N + 2, M + 2));
      T_before1.push_back(vec[indexes.back() - 1 + (M + 2 - 1) * (N + 2)]);
      T_after1.push_back(vec[indexes.back() + (M + 2 - 1) * (N + 2)]);
      a.push_back(a_value);
      q.push_back(q_value);
      m.push_back(m_value);
    }

    if ((int)filename.find("Grid_NxM-2432x20") != -1 &&
        (int)filename.find("_Lxa-120x10") != -1) {
      FileReader readfile;
      filename = filename.substr(filename.find("Grid"), filename.length());
      readfile.readGrid(entry.path(), 2, XC, YC, 2432, 20);
    }
  }
}
#endif