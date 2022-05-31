#include "FileReader.hpp"

FileReader::FileReader() {}

FileReader::~FileReader() {}

void FileReader::readField(string &name, int blockWanted, int variableWanted,
                           Field &vec, int locIStr, int locIEnd) {
  char cstr[name.size() + 2];
  strcpy(cstr, name.c_str());

  std::ifstream infile(cstr, std::ios::binary);
  // std::ifstream readGrid("Grid.xyz", std::ios::binary);

  // if (!infile.is_open())
  // {
  //   cout << "The file " << cstr << " could not be open." << endl;
  // }

  // readGrid.close();
  int nBlocs;
  infile.read((char *)&nBlocs, sizeof(nBlocs));

  vector<int> NX(nBlocs);
  vector<int> NY(nBlocs);
  vector<int> nVar(nBlocs);

  int NI = vec.NI;
  int NJ = vec.NJ / 2;

  for (int k = 0; k < nBlocs; k++) {
    infile.read((char *)&NX[k], sizeof(NX[k]));
    infile.read((char *)&NY[k], sizeof(NY[k]));
    infile.read((char *)&nVar[k], sizeof(nVar[k]));
  }
  for (int k = 0; k < nBlocs; k++)
    for (int l = 0; l < nVar[k]; l++) {
      double inputField[NX[k] * NY[k]];

      for (int j = 0; j < NY[k]; j++)
        for (int i = 0; i < NX[k]; i++) {
          double inputvalue;
          infile.read((char *)&inputvalue, sizeof(double));
          inputField[i + j * NX[k]] = inputvalue;
        }

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
