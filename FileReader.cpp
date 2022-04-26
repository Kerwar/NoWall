#include "FileReader.hpp"

FileReader::FileReader()
{
}

FileReader::~FileReader()
{
}

void FileReader::readField(string &name, int blockWanted, int variableWanted, Field &vec, int locIStr, int locIEnd)
{

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

  for (int k = 0; k < nBlocs; k++)
  {
    infile.read((char *)&NX[k], sizeof(NX[k]));
    infile.read((char *)&NY[k], sizeof(NY[k]));
    infile.read((char *)&nVar[k], sizeof(nVar[k]));
  }
  for (int k = 0; k < nBlocs; k++)
    for (int l = 0; l < nVar[k]; l++)
    {
      double inputField[NX[k] * NY[k]];

      for (int j = 0; j < NY[k]; j++)
        for (int i = 0; i < NX[k]; i++)
        {
          double inputvalue;
          infile.read((char *)&inputvalue, sizeof(double));
          inputField[i + j * NX[k]] = inputvalue;
        }

      if (blockWanted == k + 1 && variableWanted == l + 1)
      {
        double iRel = (double)NX[k] / (locIEnd - locIStr);
        double jRel = (double)NY[k] / NJ;
        double ifactor = 0;
        double jfactor = 0;
        
        double checking[vec.NI * vec.NJ];

        for(int i = 0; i < NI; i++)
          for(int j = 0; j < NJ; j++)
            checking[i + j*NI] = i + j *NI;

        for (int i = locIStr; i < locIEnd; i++)
          for (int j = NJ * k; j < NJ * blockWanted; j++)
          {
            int inputi = std::min(NX[k] - 1, (int)std::floor((i - locIStr) * iRel));
            int inputj = std::min(NY[k] - 1, (int)std::floor((j - NJ * k) * jRel));

            int previnputi = std::max((int)std::floor((i - 1 - locIStr) * iRel), 0);
            int previnputj = std::max((int)std::floor((j - 1 - NJ * k) * jRel), 0);

            ifactor = previnputi == inputi ? 1.0 : 0.5;
            jfactor = previnputj == inputj ? 1.0 : 0.5;
            
            checking[i + (j-NJ*k)*NI] = 0;
            vec.value[i + j * NI] = inputField[inputi + inputj * NX[k]] * ifactor * jfactor +
                                    inputField[previnputi + inputj * NX[k]] * (1.0 - ifactor) * jfactor +
                                    inputField[inputi + previnputj * NX[k]] * (1.0 - jfactor) * ifactor +
                                    inputField[previnputi + previnputj * NX[k]] * (1.0 - ifactor) * (1.0 - jfactor);
          }
        
        for(int i = 0; i < NI; i++)
          for(int j = 0; j < NJ; j++)
            if (checking[i + j*NI] !=0)
            std::cout << i << " " << j << " " << checking[i+j*NI];
      }
    }
  infile.close();
}
