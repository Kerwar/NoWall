#include "FileReader.h"

FileReader::FileReader()
{
}

FileReader::~FileReader()
{
}

void FileReader::readField(string &name, int blockWanted, int variableWanted, Field::vectorField &vec, int locIStr, int locIEnd)
{

  char cstr[name.size() + 2];
  strcpy(cstr, name.c_str());

  std::ifstream infile(cstr, std::ios::binary);
  std::ifstream readGrid("Grid.xyz", std::ios::binary);

  if (!infile.is_open())
  {
    cout << "The file " << cstr << " could not be open." << endl;
  }

  readGrid.close();
  int nBlocs;
  infile.read((char *)&nBlocs, sizeof(nBlocs));

  vector<int> NX(nBlocs);
  vector<int> NY(nBlocs);
  vector<int> nVar(nBlocs);

  int N = vec.size();
  int M = vec[0].size() / 2;

  for (int k = 0; k < nBlocs; k++)
  {
    infile.read((char *)&NX[k], sizeof(NX[k]));
    infile.read((char *)&NY[k], sizeof(NY[k]));
    infile.read((char *)&nVar[k], sizeof(nVar[k]));
  }

  for (int k = 0; k < nBlocs; k++)
    for (int l = 0; l < nVar[k]; l++)
    {
      Field::vectorField inputField(NX[k], Field::vec1dfield(NY[k]));

      for (int j = 0; j < NY[k]; j++)
        for (int i = 0; i < NX[k]; i++)
        {
          double inputvalue;
          infile.read((char *)&inputvalue, sizeof(inputvalue));
          inputField[i][j].value = inputvalue;
        }

      if (blockWanted == k + 1 && variableWanted == l + 1)
      {
        double iRel = (double)NX[k] / (locIEnd - locIStr);
        double jRel = (double)NY[k] / M;
        double ifactor = 0;
        double jfactor = 0;
        for (int i = locIStr; i < locIEnd; i++)
          for (int j = M * k; j < M * blockWanted; j++)
          {
            int inputi = std::min(NX[k] - 1, (int)  std::floor((i - locIStr)* iRel));
            int inputj = std::min(NY[k] - 1, (int)  std::floor((j - M * k) * jRel));

            int previnputi = std::max((int)std::floor((i - 1 - locIStr)* iRel),0);
            int previnputj = std::max((int)std::floor((j - 1 - M * k) * jRel),0);
            
            ifactor = previnputi == inputi ? 1.0 : 0.5;
            jfactor = previnputj == inputj ? 1.0 : 0.5;
            
            vec[i][j].value = inputField[inputi][inputj].value * ifactor * jfactor 
            + inputField[previnputi][inputj].value * (1.0 - ifactor) * jfactor
            + inputField[inputi][previnputj].value * (1.0 - jfactor) * ifactor
            + inputField[previnputi][previnputj].value * (1.0 - ifactor) * (1.0 - jfactor);
          }
      }
    }
  infile.close();
}
