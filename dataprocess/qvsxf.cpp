#include <filesystem>
#include <iostream>
#include <string>

#include "datamanager.hpp"
#include "field.hpp"

int main() {
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

  DataManager data_processor(path);

  for (const auto &entry : std::filesystem::directory_iterator(path)) {
    string filename = entry.path();

    if (filename.find("Sol_NxM") != string::npos &&
        filename.find("_-1.f") != string::npos) {
      data_processor.add_member();
      data_processor.get_dimensions(filename);
      data_processor.get_parameters(filename);
      std::cout << filename << std::endl;
      data_processor.get_flame_position(filename);
    }
  }

  std::cout << "ALL DATA READ\n";

  while (data_processor.size() > 0) {
    // Choose a value of a to process
    double a_active = data_processor.database[0].a;
    std::cout << data_processor.size() << "\n";

    // Set the file where the info will go
    std::ostringstream filename;
    std::streamsize ss = filename.precision();
    filename << "/home/javier/Documentos/Resultados/CF2D/Q-vs-xf/q~xf_a-"
             << std::fixed << std::setprecision(0) << a_active
             << "_m-2_beta-10.txt";
    filename << std::setprecision(ss);

    std::ofstream outfile;
    outfile.open(filename.str());
    std::cout << filename.str() << std::endl;

    // Take the indexes of the files that have the right a value
    std::vector<int> index_files_active;
    std::vector<int> index_files_a_remove;

    index_files_active = data_processor.filter_by_a(a_active);

    std::cout << "Files filter!\n";
    // Loop over the right files to get the one with maximum q

    data_processor.write_data(outfile, index_files_active);
    std::cout << "Data Written!\n";
    data_processor.remove_data(index_files_active);
    std::cout << "Data deleted!\n";
    outfile.close();
  }
}
