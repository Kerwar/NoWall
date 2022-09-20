#include "datamanager.hpp"

DataManager::DataManager(string path) : path_(path) {}

int DataManager::indexofT1(const Field &vec, const int &ni, const int &nj) {
  int result = ni * (nj/2 - 1);

  while (vec.value[result] < 1 && result < nj * ni) {
    // std::cout << result << " " << vec[result];
    result++;
  }

  if (result >= ni * nj) {
    result = ni * (nj - 1);
    int best_index = 0;

    for (int i = 0; i < ni; i++) {
      // std::cout <<  i << " " << vec.value[i + (ni-1)*nj] << " " << best_index
      // << std::endl;
      if (std::abs(1.0 - vec[result + i]) <
          std::abs(1.0 - vec[result + best_index]))
        best_index = i;
    }
    result = best_index;
    // std::cout << "One problem" << std::endl;
  }
  return result;
}

int DataManager::max_Z(const Field &vec, const int &ni, const int &nj) {
  int result = ni * (nj/2 - 1);

  int best_index = 0;

  for (int i = 0; i < ni; i++) {
    if (std::abs(1.0 - vec[result + i]) <
        std::abs(1.0 - vec[result + best_index]))
      best_index = i;
  }
  result = best_index;
  return result;
}

double DataManager::interpolate_quadratic_maximum(
    const double &x1, const double &x2, const double &x3, const double &y1,
    const double &y2, const double &y3) {
  return 0.5 *
         (pow(x1, 2) * y2 - pow(x1, 2) * y3 - pow(x2, 2) * y1 +
          pow(x2, 2) * y3 + pow(x3, 2) * y1 - pow(x3, 2) * y2) /
         (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
}

string DataManager::find_parameter_value(const string fname,
                                         const string &preString,
                                         const string &postString) {
  return fname.substr(
      fname.find(preString) + preString.length(),
      fname.find(postString) - fname.find(preString) - preString.length());
}

// Sol_NxM-2432x20_Lxa-120x1_Ex-40_q-1.2_m-2_beta-10_LeFxZ-1x0.3_-1.q
void DataManager::get_dimensions(string filename) {
  string dims = find_parameter_value(filename, "NxM-", "_Lxa");
  if (XC.find(dims) == XC.end()) add_grid(filename, dims);

  string adims = find_parameter_value(filename, "_Lxa-", "_Ex");

  database.back().a =
      stod(adims.substr(adims.find("x") + 1, adims.length() - 1));
}

void DataManager::add_grid(string filename, string dims) {
  string gridFile = "/" + find_parameter_value(filename, "/", "Sol") + "Grid_" +
                    dims + ".xyz";
  FileReader readfile;

  vector<double> xc, yc;
  readfile.read_grid(gridFile, xc, yc);

  XC[dims] = xc;
  YC[dims] = yc;
}

void DataManager::get_parameters(string filename) {
  string qstring = find_parameter_value(filename, "_q-", "_m-");
  database.back().q = stod(qstring);

  string mstring = find_parameter_value(filename, "_m-", "_beta-");
  database.back().m = stod(mstring);
}

void DataManager::get_flame_position(string filename) {
  string dims = find_parameter_value(filename, "NxM-", "_Lxa");

  Data &active = database.back();
  int n = XC[dims].size();
  int m = YC[dims].size();

  FileReader filereader;

  Field vec(n, m);

  filereader.read_exact_field(filename, 1, 0, vec);

  int index = indexofT1(vec, n, m);
  active.indexes = index;

  active.T_before1 = vec[index - 1];
  active.T_after1 = vec[index];
  active.T = (XC[dims][index % n] - XC[dims][index % n - 1] -
              active.T_before1 * XC[dims][index%n] +
              active.T_after1 * XC[dims][index%n - 1]) /
             (active.T_after1 - active.T_before1);
  if(active.T > XC[dims][index%n]) active.T = XC[dims][index%n];
  if(active.T < XC[dims][index%n - 1]) active.T = XC[dims][index%n];

  filereader.read_exact_field(filename, 1, 2, vec);

  index = max_Z(vec, n, m);
  active.indexes = index;

  active.Z_Max_prev = vec[index - 1];
  active.Z_Max = vec[index];
  active.Z_Max_post = vec[index + 1];

  active.Z = interpolate_quadratic_maximum(
      XC[dims][index % n - 1], XC[dims][index % n], XC[dims][index % n + 1],
      active.Z_Max_prev, active.Z_Max, active.Z_Max_post);
}

vector<int> DataManager::filter_by_a(double a) {
  vector<int> result;
  for (int i = 0; i < size(); i++)
    if (std::abs(database[i].a - a) < 10e-5) result.push_back(i);
  
  return result;
}

vector<int> DataManager::filter_by_m(double m) {
  vector<int> result;
  for (int i = 0; i < size(); i++)
    if (std::abs(database[i].m - m) < 10e-5) result.push_back(i);
  
  return result;
}

void DataManager::write_data(std::ofstream &file, vector<int> data_index) {
  
  while (data_index.size() > 0) {
    int add_i = data_index[0];
    vector<int> copyhelp = data_index;

    for (int k : data_index) 
      if (database[k].T < database[add_i].T) 
        add_i = k;
    
    file << database[add_i].q << " " << database[add_i].T << " " << database[add_i].Z << " " 
         << database[add_i].a << " " << database[add_i].m << std::endl;

    data_index.clear();
    remove_copy(copyhelp.begin(), copyhelp.end(), std::back_inserter(data_index), add_i);
  }
}

void DataManager::remove_data(vector<int> data_index) {
  for (int i = 0; i < (int)data_index.size(); i++) {
    database.erase(database.begin() + data_index[i]);
    for (int k = i; k < (int)data_index.size(); k++)
      if (data_index[i] < data_index[k]) data_index[k]--;
  }
}