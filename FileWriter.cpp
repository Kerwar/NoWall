#include "FileWriter.h"

FileWriter::FileWriter()
{
  // ctor
}

FileWriter::~FileWriter()
{
  // dtor
}

// void FileWriter::writeTFZ(string &name, int time, Grid &myGrid,
//                           Field::vectorField &Ttemp, Field::vectorField &Ftemp,
//                           Field::vectorField &Ztemp, int iStr, int iEnd,
//                           int jStr, int jEnd, Paralel::Loc loc)
// {
//   string myFileType = ".dat";

//   int ttime = time;

//   string timename;
//   string location;

//   switch (loc)
//   {
//   case Paralel::up:
//     location = "Up";
//     break;

//   case Paralel::down:
//     location = "Down";
//     break;

//     // case Paralel::wall:
//     //   location = "Wall";
//     //   break;
//   }

//   std::ostringstream temp;
//   temp << ttime;

//   timename = temp.str();

//   string newname = name;

//   string newfilename = newname.append(timename);
//   newfilename = newfilename.append(location);
//   string name2 = newfilename.append(myFileType);

//   char cstr[name2.size() + 2];
//   strcpy(cstr, newfilename.c_str());

//   std::ofstream outfile;

//   outfile.open(cstr, std::ios::out | std::ios::binary);

//   int NXtemp = iEnd - iStr;
//   int NYtemp = jEnd - jStr;

//   outfile << "VARIABLES=\"X\", \"Y\", \"T\", \"F\", \"Z\" \n";
//   outfile << "ZONE T=\"" << location << "\" ,I=" << NXtemp + 1 << ", J=" << NYtemp + 1 << ", DATAPACKING=BLOCK, VARLOCATION=([3-5]=CELLCENTERED)\n";

//   for (int j = jStr; j < jEnd + 1; j++)
//   {
//     for (int i = iStr; i < iEnd + 1; i++)
//     {

//       double xpos = myGrid.X[i][j];
//       outfile << prd(xpos, 2, 7);
//     }
//   }
//   outfile << std::endl;
//   for (int j = jStr; j < jEnd + 1; j++)
//   {
//     for (int i = iStr; i < iEnd + 1; i++)
//     {

//       double ypos = myGrid.Y[i][j];
//       outfile << prd(ypos, 2, 7);
//     }
//   }
//   outfile << std::endl;

//   for (int j = jStr; j < jEnd; j++)
//   {
//     for (int i = iStr; i < iEnd; i++)
//     {

//       double TT = Ttemp[i][j].value;
//       outfile << prd(TT, 4, 10);
//     }
//   }
//   outfile << std::endl;
//   for (int j = jStr; j < jEnd; j++)
//   {
//     for (int i = iStr; i < iEnd; i++)
//     {

//       double FF = Ftemp[i][j].value;
//       outfile << prd(FF, 4, 10);
//     }
//   }
//   outfile << std::endl;
//   for (int j = jStr; j < jEnd; j++)
//   {
//     for (int i = iStr; i < iEnd; i++)
//     {

//       double ZZ = Ztemp[i][j].value;
//       outfile << prd(ZZ, 4, 10);
//     }
//   }
//   outfile << std::endl;

//   outfile.close();
// }

// std::string FileWriter::prd(const double x, const int decDigits, const int width)
// {
//   std::stringstream ss;
//   ss << std::fixed << std::right;
//   ss.fill(' ');            // fill space around displayed #
//   ss.width(width);         // set  width around displayed #
//   ss.precision(decDigits); // set # places after decimal
//   ss << x;
//   return ss.str();
// }

void FileWriter::WriteBin(string prefix, string sufix, int time, Grid &myGrid, Field &Utemp, Field &Vtemp, Field &Ttemp, Field &Ftemp, Field &Ztemp, int iStr, int iEnd, int jStr, int jEnd, Paralel::Loc loc)
{

  int solNI = Utemp.NI;
  int NI = myGrid.NI;

  // Needed for q file
  double mach = 1;
  double alpha = 1;
  double reyn = 1;
  double qtime = 0;
  int ttime = time;

  string timename;

  if (time == 0)
  {
    if (loc == Paralel::down)
    {
      string myGridType = ".xyz";

      string newname = "Grid";

      string name2 = prefix + newname;
      
      name2.append(sufix);
      name2.append(myGridType);

      char cstr[name2.size() + 2];
      strcpy(cstr, name2.c_str());

      std::ofstream outfile;

      outfile.open(cstr, std::ios::binary);

      int Blocks = 2;
      outfile.write((char *)&Blocks, sizeof(Blocks));

      int NXtemp = myGrid.N;
      int NYtemp = jEnd - jStr + 1;

      outfile.write((char *)&NXtemp, sizeof(NXtemp));
      outfile.write((char *)&NYtemp, sizeof(NYtemp));

      NXtemp = myGrid.N;
      NYtemp = jEnd - jStr + 1;

      outfile.write((char *)&NXtemp, sizeof(NXtemp));
      outfile.write((char *)&NYtemp, sizeof(NYtemp));

      for (int j = jStr; j <= jEnd; j++)
        for (int i = 0; i < iEnd; i++)
          outfile.write(reinterpret_cast<char *>(&myGrid.XC[i + NI * j]), sizeof(myGrid.XC[i + NI * j]));
      
      for (int j = jStr; j <= jEnd; j++)
        for (int i = 0; i < iEnd; i++)
          outfile.write(reinterpret_cast<char *>(&myGrid.YC[i + NI * j]), sizeof(myGrid.YC[i + NI * j]));
      
      outfile.close();
    }
    else
    {
      string myGridType = ".xyz";

      string newname = "Grid";

      string name2 = prefix + newname;
      name2.append(sufix);
      name2.append(myGridType);

      char cstr[name2.size() + 2];
      strcpy(cstr, name2.c_str());

      std::ofstream outfile;
      outfile.open(cstr, std::ios::app | std::ios::binary);


      for (int j = jStr; j <= jEnd; j++)
        for (int i = iStr; i < iEnd; i++)
          outfile.write(reinterpret_cast<char *>(&myGrid.XC[i + NI * j]), sizeof(myGrid.XC[i + NI * j]));
      
      for (int j = jStr; j <= jEnd; j++)
        for (int i = iStr; i < iEnd; i++)
          outfile.write(reinterpret_cast<char *>(&myGrid.YC[i + NI * j]), sizeof(myGrid.YC[i + NI * j]));
      
      outfile.close();
    }
  }

  std::ostringstream temp;
  temp << ttime;

  timename = temp.str();

  string newname = "Sol";

  string name = prefix + newname;
  name.append(sufix);

  string newfilename = name.append(timename);

  std::ofstream outfile;

  string myType = ".q";

  string name2 = newfilename.append(myType);
  char qstr[name2.size() + 2];
  strcpy(qstr, name2.c_str());

  int Blocks = 2;
  int NXtemp = myGrid.N;
  int NYtemp = jEnd - jStr + 1;
  if (loc == Paralel::down)
  {

    outfile.open(qstr, std::ios::binary);

    outfile.write((char *)&Blocks, sizeof(Blocks));

    NXtemp = myGrid.N;
    NYtemp = jEnd - jStr + 1;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));

    NXtemp = myGrid.N;
    NYtemp = jEnd - jStr + 1;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
  }
  else
  {
    outfile.open(qstr, std::ios::app | std::ios::binary);
  }
  // In q file there are 4 variables:
  outfile.write((char *)&mach, sizeof(mach));
  outfile.write((char *)&alpha, sizeof(alpha));
  outfile.write((char *)&reyn, sizeof(reyn));
  outfile.write((char *)&qtime, sizeof(qtime));

  // 1) Density
  // 2) U
  // 3) V
  // 4) Energy
  double rho = 1;
  for (int j = jStr; j <= jEnd; j++)
  {
    for (int i = iStr; i < iEnd; i++)
    {
      outfile.write(reinterpret_cast<char *>(&rho), sizeof(rho));
    }
  }
  for (int j = jStr; j <= jEnd; j++)
  {
    for (int i = iStr; i < iEnd; i++)
    {
      outfile.write(reinterpret_cast<char *>(&Utemp.value[i + j * solNI]), sizeof(Utemp.value[i + j * solNI]));
    }
  }
  for (int j = jStr; j <= jEnd; j++)
  {
    for (int i = iStr; i < iEnd; i++)
    {
      outfile.write(reinterpret_cast<char *>(&Vtemp.value[i + j * solNI]), sizeof(Vtemp.value[i + j * solNI]));
    }
  }
  for (int j = jStr; j <= jEnd; j++)
  {
    for (int i = iStr; i < iEnd; i++)
    {
      outfile.write(reinterpret_cast<char *>(&Ttemp.value[i + j * solNI]), sizeof(Ttemp.value[i + j * solNI]));
      
    }
  }
  outfile.close();

  newname = "Sol";

  name = prefix + newname;
  name.append(sufix);

  newfilename = name.append(timename);

  string myFileType = ".f";

  string name3 = newfilename.append(myFileType);
  char cstr[name3.size() + 2];
  strcpy(cstr, newfilename.c_str());
  if (loc == Paralel::down)
  {

    outfile.open(cstr, std::ios::out | std::ios::binary);

    int NVar = 3;

    outfile.write((char *)&Blocks, sizeof(Blocks));
    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
    outfile.write((char *)&NVar, sizeof(NVar));

    NXtemp = myGrid.exI2 - myGrid.exI1;
    NYtemp = jEnd - jStr + 1;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
    outfile.write((char *)&NVar, sizeof(NVar));

    NXtemp = myGrid.N;
    NYtemp = jEnd - jStr + 1;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
    outfile.write((char *)&NVar, sizeof(NVar));
  }
  else
  {
    outfile.open(cstr, std::ios::app | std::ios::binary);
  }
  for (int j = jStr; j <= jEnd; j++)
  {
    for (int i = iStr; i < iEnd; i++)
    {
      outfile.write(reinterpret_cast<char *>(&Ttemp.value[i + j * solNI]), sizeof(Ttemp.value[i + j * solNI]));
    }
  }
  for (int j = jStr; j <= jEnd; j++)
  {
    for (int i = iStr; i < iEnd; i++)
    {
      outfile.write(reinterpret_cast<char *>(&Ftemp.value[i + j * solNI]), sizeof(Ftemp.value[i + j * solNI]));
    }
  }
  for (int j = jStr; j <= jEnd; j++)
  {
    for (int i = iStr; i < iEnd; i++)
    {
      outfile.write(reinterpret_cast<char *>(&Ztemp.value[i + j * solNI]), sizeof(Ztemp.value[i + j * solNI]));
    }
  }
  outfile.close();
}

void FileWriter::WriteInter(string prefix, string sufix, int time, Grid &mainGrid, Grid &myGrid, Field &Utemp, Field &Vtemp, Field &Ttemp, Field &Ftemp, Field &Ztemp, int iStr, int iEnd, int jStr, int jEnd, Paralel::Loc loc)
{
  // Needed for q file
  double mach = 1;
  double alpha = 1;
  double reyn = 1;
  double qtime = 0;
  int ttime = time;

  int mainNI = mainGrid.NI;
  int solNI = Utemp.NI;
  int NJ = myGrid.NJ;

  string timename;

  if (time == 0)
  {
    if (loc == Paralel::down)
    {
      string myGridType = ".xyz";

      string newname = "Grid";

      string name2 = prefix + newname;
      
      name2.append(sufix);
      name2.append(myGridType);

      char cstr[name2.size() + 2];
      strcpy(cstr, name2.c_str());

      std::ofstream outfile;

      outfile.open(cstr, std::ios::binary);

      int Blocks = 2;
      outfile.write((char *)&Blocks, sizeof(Blocks));

      int NXtemp = mainNI - 1;
      int NYtemp = NJ - 1;

      outfile.write((char *)&NXtemp, sizeof(NXtemp));
      outfile.write((char *)&NYtemp, sizeof(NYtemp));

      
      NXtemp = mainNI - 1;
      NYtemp = NJ - 1;

      outfile.write((char *)&NXtemp, sizeof(NXtemp));
      outfile.write((char *)&NYtemp, sizeof(NYtemp));

      for (int j = jStr; j < NJ + jStr - 1; j++)
        for (int i = 0; i < mainNI - 1; i++)
          outfile.write(reinterpret_cast<char *>(&mainGrid.X[i + mainNI * j]), sizeof(mainGrid.X[i + mainNI * j]));
      
      for (int j = jStr; j < NJ + jStr - 1; j++)
        for (int i = 0; i < mainNI - 1; i++)
          outfile.write(reinterpret_cast<char *>(&mainGrid.Y[i + mainNI * j]), sizeof(mainGrid.Y[i + mainNI * j]));
      
      outfile.close();
    }
    else if (loc == Paralel::up)
    {
      string myGridType = ".xyz";

      string newname = "Grid";

      string name2 = prefix + newname;
      name2.append(sufix);
      name2.append(myGridType);

      char cstr[name2.size() + 2];
      strcpy(cstr, name2.c_str());

      std::ofstream outfile;
      outfile.open(cstr, std::ios::app | std::ios::binary);

      int jMax = myGrid.NJ + jStr - 1;

      for (int j = jStr; j < jMax; j++)
        for (int i = 0; i < mainNI - 1; i++)
          outfile.write(reinterpret_cast<char *>(&mainGrid.X[i + mainNI * j]), sizeof(mainGrid.X[i + mainNI * j]));
      
      for (int j = jStr; j < jMax; j++)
        for (int i = 0; i < mainNI - 1; i++)
          outfile.write(reinterpret_cast<char *>(&mainGrid.Y[i + mainNI * j]), sizeof(mainGrid.Y[i + mainNI * j]));
      
      outfile.close();
    }
  }

  std::ostringstream temp;
  temp << ttime;

  timename = temp.str();

  string newname = "Sol";

  string name = prefix + newname;
  name.append(sufix);

  string newfilename = name.append(timename);

  std::ofstream outfile;

  string myType = ".q";

  string name2 = newfilename.append(myType);
  char qstr[name2.size() + 2];
  strcpy(qstr, name2.c_str());
  
  int Blocks = 2;
  int NXtemp = mainGrid.NI - 1;
  int NYtemp = myGrid.NJ - 1;
  if (loc == Paralel::down)
  {

    outfile.open(qstr, std::ios::binary);

    outfile.write((char *)&Blocks, sizeof(Blocks));

    NXtemp = mainGrid.NI - 1;
    NYtemp = myGrid.NJ - 1;
    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));

    NXtemp = mainGrid.NI - 1;
    NYtemp = myGrid.NJ - 1;
    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
  }
  else
  {
    outfile.open(qstr, std::ios::app | std::ios::binary);
  }
  // In q file there are 4 variables:
  outfile.write((char *)&mach, sizeof(mach));
  outfile.write((char *)&alpha, sizeof(alpha));
  outfile.write((char *)&reyn, sizeof(reyn));
  outfile.write((char *)&qtime, sizeof(qtime));

  // 1) Density
  // 2) U
  // 3) V
  // 4) Energy
  double rho = 1;
  int jMax = myGrid.NJ + jStr - 1;
  for (int j = jStr; j < jMax; j++)
  {
    for (int i = iStr; i < iEnd + 1; i++)
    {
      outfile.write(reinterpret_cast<char *>(&rho), sizeof(rho));
    }
  }
  for (int j = jStr; j < jMax; j++)
  {
    for (int i = iStr; i < iEnd + 1; i++)
    {
      double value = Utemp.value[std::min(i, iEnd - 1) + std::min(j, jMax - 2) * solNI];
      if (i != iStr && i != iEnd)
        value = 0.5 * (value + Utemp.value[i - 1 + std::min(j, jMax - 2) * solNI]);
      if (j != jStr && j != jMax - 1)
        value = 0.5 * (value + Utemp.value[std::min(i, iEnd - 1) + (j - 1) * solNI]);
      
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  }
  for (int j = jStr; j < jMax; j++)
  {
    for (int i = iStr; i < iEnd + 1; i++)
    {
      double value = Vtemp.value[std::min(i, iEnd - 1) + std::min(j, jMax - 2) * solNI];
      if (i != iStr && i != iEnd)
        value = 0.5 * (value + Vtemp.value[i - 1 + std::min(j, jMax - 2) * solNI]);
      if (j != jStr && j != jMax - 1)
        value = 0.5 * (value + Vtemp.value[std::min(i, iEnd - 1) + (j - 1) * solNI]);
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  }

  for (int j = jStr; j < jMax; j++)
  {
    for (int i = iStr; i < iEnd + 1; i++)
    {
      double value = Ttemp.value[std::min(i, iEnd - 1) + std::min(j, jMax - 2) * solNI];
      if (i != iStr && i != iEnd)
        value = 0.5 * (value + Ttemp.value[i - 1 + std::min(j, jMax - 2) * solNI]);
      if (j != jStr && j != jMax - 1)
        value = 0.5 * (value + Ttemp.value[std::min(i, iEnd - 1) + (j - 1) * solNI]);
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  }
  outfile.close();

  newname = "Sol";

  name = prefix + newname;
  name.append(sufix);

  newfilename = name.append(timename);

  string myFileType = ".f";

  string name3 = newfilename.append(myFileType);
  char cstr[name3.size() + 2];
  strcpy(cstr, newfilename.c_str());
  if (loc == Paralel::down)
  {

    outfile.open(cstr, std::ios::out | std::ios::binary);

    int NVar = 3;

    outfile.write((char *)&Blocks, sizeof(Blocks));

    NXtemp = mainGrid.NI - 1;
    NYtemp = myGrid.NJ - 1;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
    outfile.write((char *)&NVar, sizeof(NVar));

    NXtemp = mainGrid.NI - 1;
    NYtemp = myGrid.NJ - 1;

    outfile.write((char *)&NXtemp, sizeof(NXtemp));
    outfile.write((char *)&NYtemp, sizeof(NYtemp));
    outfile.write((char *)&NVar, sizeof(NVar));
  }
  else
  {
    outfile.open(cstr, std::ios::app | std::ios::binary);
  }
  for (int j = jStr; j < jMax; j++)
  {
    for (int i = iStr; i < iEnd + 1; i++)
    {
      double value = Ttemp.value[std::min(i, iEnd - 1) + std::min(j, jMax - 2) * solNI];
      if (i != iStr && i != iEnd)
        value = 0.5 * (value + Ttemp.value[i - 1 + std::min(j, jMax - 2) * solNI]);
      if (j != jStr && j != jMax - 1)
        value = 0.5 * (value + Ttemp.value[std::min(i, iEnd - 1) + (j - 1) * solNI]);
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  }
  for (int j = jStr; j < jMax; j++)
  {
    for (int i = iStr; i < iEnd + 1; i++)
    {
      double value = Ftemp.value[std::min(i, iEnd - 1) + std::min(j, jMax - 2) * solNI];
      if (i != iStr && i != iEnd)
        value = 0.5 * (value + Ftemp.value[i - 1 + std::min(j, jMax - 2) * solNI]);
      if (j != jStr && j != jMax - 1)
        value = 0.5 * (value + Ftemp.value[std::min(i, iEnd - 1) + (j - 1) * solNI]);
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  }
  for (int j = jStr; j < jMax; j++)
  {
    for (int i = iStr; i < iEnd + 1; i++)
    {
      double value = Ztemp.value[std::min(i, iEnd - 1) + std::min(j, jMax - 2) * solNI];
      if (i != iStr && i != iEnd)
        value = 0.5 * (value + Ztemp.value[i - 1 + std::min(j, jMax - 2) * solNI]);
      if (j != jStr && j != jMax - 1)
        value = 0.5 * (value + Ztemp.value[std::min(i, iEnd - 1) + (j - 1) * solNI]);
      outfile.write(reinterpret_cast<char *>(&value), sizeof(value));
    }
  }
  outfile.close();
}