#include <filesystem>
#include <iostream>
#include <string>

#include "Field.hpp"
#include "readdata.hpp"

int main() {
  std::vector<double> a, q, m;
  std::vector<double> T_before1, T_after1;
  std::vector<int> indexes;
  std::vector<double> XC(2434, 0), YC(12, 0);

  // std::string path = "/home/javier/Documentos/Resultados/CF2D/A-vs-X_f";
  // readdata(path, a, T_before1, T_after1, XC, YC, q, m, indexes);

  // std::ofstream outfile;
  // outfile.open(
  //     "/home/javier/Documentos/Resultados/CF2D/A-vs-X_f/a~xf_m-2_beta_10.txt");

  // int n_files = a.size();

  // for (int i = 0; i < n_files; i++) {
  //   int min_i = 0;
  //   for (int k = 1; k < (int)a.size(); k++)
  //     if (a[k] < a[min_i]) min_i = k;

  //   double X = (XC[indexes[min_i]] - XC[indexes[min_i] - 1] -
  //               T_before1[min_i] * XC[indexes[min_i]] +
  //               T_after1[min_i] * XC[indexes[min_i]-1]) /
  //              (T_after1[min_i] - T_before1[min_i]);
  //   std::cout << a[min_i] << " " << X << " " << XC[indexes[min_i]] <<
  //   std::endl; outfile << a[min_i] << " " << X << std::endl;
  //   a.erase(a.begin() + min_i);
  //   T_before1.erase(T_before1.begin() + min_i);
  // T_after1.erase(T_after1.begin() + min_i);
  //   indexes.erase(indexes.begin() + min_i);
  // }
  // outfile.close();

  std::string path = "/home/javier/Documentos/Resultados/CF2D/Q-vs-xf";
  readdata(path, a, T_before1, T_after1, XC, YC, q, m, indexes);

  while (a.size() > 0) {
    // Choose a value of a to process
    double a_active = a[0];
    std::cout << a.size() << "\n";

    // Set the file where the info will go
    std::ofstream outfile;
    std::ostringstream filename;
    filename << "/home/javier/Documentos/Resultados/CF2D/Q-vs-xf/q~xf_a-"
             << std::fixed << std::setprecision(0) << a_active
             << "_m-2_beta-10.txt";
    outfile.open(filename.str());
    std::cout << filename.str() << std::endl;

    // Take the indexes of the files that have the right a value
    std::vector<int> index_files_a_active;
    std::vector<int> index_files_a_remove;

    for (int i = 0; i < (int)a.size(); i++)
      if (a[i] == a_active) index_files_a_active.push_back(i);

    // Loop over the right files to get the one with maximum q
    int n_files = index_files_a_active.size();

    for (int i = 0; i < n_files; i++) {
      int K = 0;
      int add_i = index_files_a_active[K];

      for (int k = 1; k < (int) index_files_a_active.size(); k++)
        if (q[index_files_a_active[k]] > q[add_i]) {
          K = k;
          add_i = index_files_a_active[K];
        }

      double X = (XC[indexes[add_i]] - XC[indexes[add_i] - 1] -
                  T_before1[add_i] * XC[indexes[add_i]] +
                  T_after1[add_i] * XC[indexes[add_i] - 1]) /
                 (T_after1[add_i] - T_before1[add_i]);

      if (T_after1[add_i] < 1.0) X = XC[indexes[add_i]];
      outfile << q[add_i] << " " << X << " " << a[add_i] << " " << m[add_i]
              << std::endl;

      index_files_a_remove.push_back(add_i);
      index_files_a_active.erase(index_files_a_active.begin() + K);
    }
    for (int i = 0; i < n_files; i++) {
      int del_i = index_files_a_remove[i];
      a.erase(a.begin() + del_i);
      q.erase(q.begin() + del_i);
      m.erase(m.begin() + del_i);
      T_before1.erase(T_before1.begin() + del_i);
      T_after1.erase(T_after1.begin() + del_i);
      indexes.erase(indexes.begin() + del_i);
      for (int k = i; k < n_files; k++)
        if(index_files_a_remove[i] < index_files_a_remove[k])
          index_files_a_remove[k]--;
    }
    outfile.close();
  }
}
