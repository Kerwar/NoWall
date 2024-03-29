#include "filereader.hpp"

FileReader::FileReader() {}

FileReader::~FileReader() {}

void FileReader::read_field(const string &name, int blockWanted,
                            int variableWanted, Field &vec, int locIStr,
                            int locIEnd) {
  std::ifstream infile(name, std::ios::binary);

  vector<int> NX, NY, nVar;

  read_file_header(infile, NX, NY, nVar);
  int nBlocs = NX.size();
  int NI = vec.NI;
  int NJ = vec.NJ / 2;

  for (int k = 0; k < nBlocs; k++)
    for (int l = 0; l < nVar[k]; l++) {
      std::vector<double> inputField(NX[k] * NY[k]);

      inputField = read_variable(infile, NX[0], NY[0]);

      if (blockWanted == k + 1 && variableWanted == l + 1) {
        double iRel = (double)NX[k] / (locIEnd - locIStr);
        double jRel = (double)NY[k] / NJ;
        for (int i = locIStr; i < locIEnd; i++)
          for (int j = NJ * k; j < NJ * blockWanted; j++) {
            const int I = i - locIStr;
            const int J = j - NJ * k;

            int inputi = (int)std::floor(I * iRel);
            int inputj = (int)std::floor(J * jRel);

            double ifactor = 1.0 - (I * iRel - inputi);
            double jfactor = 1.0 - (J * jRel - inputj);

            int nextinputi = std::min(inputi + 1, NX[k] - 1);
            int nextinputj = std::min(inputj + 1, NY[k] - 1);

            vec.value[i + j * NI] =
                inputField[inputi + inputj * NX[k]] * ifactor * jfactor +
                inputField[nextinputi + inputj * NX[k]] * (1.0 - ifactor) *
                    jfactor +
                inputField[inputi + nextinputj * NX[k]] * (1.0 - jfactor) *
                    ifactor +
                inputField[nextinputi + nextinputj * NX[k]] * (1.0 - ifactor) *
                    (1.0 - jfactor);
          }
      }
    }
  infile.close();
}

void FileReader::read_exact_field(const string &name, int blockWanted,
                                  int variableWanted, Field &vec) {
  std::ifstream infile(name, std::ios::binary);

  vector<int> NX, NY, nVar;

  read_file_header(infile, NX, NY, nVar);
  int nBlocs = NX.size();
  for (int k = 0; k < nBlocs; k++)
    for (int l = 0; l < nVar[k]; l++) {
      std::vector<double> inputField(NX[k] * NY[k]);

      inputField = read_variable(infile, NX[k], NY[k]);

      if (blockWanted == k && variableWanted == l) {
        if (vec.NI != NX[k] || vec.NJ != NY[k])
        std::cout << "Not Equal " << vec.NI << " " << NX[k] << " " << vec.NJ
                  << " " << NY[k] << "\n";
        
        vec.value = inputField;
        
      }
    }
  infile.close();
}

void FileReader::read_grid(const string &name, vector<double> &XC,
                           vector<double> &YC) {
  std::ifstream gridFile(name, std::ios::binary);

  vector<int> NX, NY;

  read_grid_header(gridFile, NX, NY);

  XC = read_variable(gridFile, NX[0], NY[0]);
  XC.resize(NX[0]);

  vector<double> all_YC = read_variable(gridFile, NX[0], NY[0]);
  YC.resize(NY[0]);
  for (int j = 0; j < NY[0]; j++) YC[j] = all_YC[NX[0] * j];

  gridFile.close();
}

void FileReader::read_file_header(std::ifstream &file, vector<int> &NX,
                                  vector<int> &NY, vector<int> &nVar) {
  int nBlocs;
  file.read((char *)&nBlocs, sizeof(nBlocs));

  NX.resize(nBlocs);
  NY.resize(nBlocs);
  nVar.resize(nBlocs);

  for (int k = 0; k < nBlocs; k++) {
    file.read((char *)&NX[k], sizeof(NX[k]));
    file.read((char *)&NY[k], sizeof(NY[k]));
    file.read((char *)&nVar[k], sizeof(nVar[k]));
  }
}

void FileReader::read_grid_header(std::ifstream &file, vector<int> &NX,
                                  vector<int> &NY) {
  int nBlocs;
  file.read((char *)&nBlocs, sizeof(nBlocs));

  NX.resize(nBlocs);
  NY.resize(nBlocs);

  for (int k = 0; k < nBlocs; k++) {
    file.read((char *)&NX[k], sizeof(NX[k]));
    file.read((char *)&NY[k], sizeof(NY[k]));
  }
}

vector<double> FileReader::read_variable(std::ifstream &file, int NX, int NY) {
  vector<double> result(NX * NY);

  for (int j = 0; j < NY; j++)
    for (int i = 0; i < NX; i++) {
      double inputvalue;
      file.read((char *)&inputvalue, sizeof(double));
      result[i + j * NX] = inputvalue;
    }

  return result;
}