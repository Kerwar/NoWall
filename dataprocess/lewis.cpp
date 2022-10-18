#include <filesystem>
#include <iostream>
#include <string>

#include "datamanager.hpp"
#include "field.hpp"

int main() {
  std::string path = "/home/javier/Documentos/Resultados/CF2D/Lewis";

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
    filename << "/home/javier/Documentos/Resultados/CF2D/Lewis/lewis_m-2_beta-10.txt";
    filename << std::setprecision(ss);

    std::ofstream outfile;
    outfile.open(filename.str());
    std::cout << filename.str() << std::endl;

    // Take the indexes of the files that have the right a value
    std::vector<int> index_files_active;
    std::vector<int> index_files_a_remove;

    index_files_active = data_processor.filter_by_a(a_active);

    // Loop over the right files to get the one with maximum q

    data_processor.write_data(outfile, index_files_active);
    data_processor.remove_data(index_files_active);
    outfile.close();
  }
}
