#include <filesystem>
#include <iostream>
#include <string>

#include "datamanager.hpp"
#include "field.hpp"

int main() {
  std::string path = "/home/javier/Documentos/Resultados/CF2D/A-vs-X_f";

  DataManager data_processor(path);

  for (const auto &entry : std::filesystem::directory_iterator(path)) {
    string filename = entry.path();
    if (filename.find("Sol_NxM") != string::npos &&
        filename.find("_-1.f") != string::npos) {
      data_processor.add_member();
      data_processor.get_dimensions(filename);
      data_processor.get_parameters(filename);
      data_processor.get_flame_position(filename);
    }
  }

  std::ofstream outfile;
  
  outfile.open(
      "/home/javier/Documentos/Resultados/CF2D/A-vs-X_f/a~xf_m-2_beta_10.txt");


  vector<int> index_files_active = data_processor.filter_by_m(2.0);
  data_processor.write_data(outfile, index_files_active);
  data_processor.remove_data(index_files_active);

  outfile.close();
}
