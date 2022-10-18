#ifndef _DATAMANAGER_HPP_
#define _DATAMANAGER_HPP_

#include <filesystem>
#include <iostream>
#include <map>
#include <iterator>
#include <string>

#include "data.hpp"
#include "filereader.hpp"

using std::stod, std::stoi;
using std::vector;

namespace fs = std::filesystem;

class DataManager {
 public:
  DataManager(string path);

  vector<Data> database;
  std::map<std::string, vector<double>> XC, YC;

  string path() { return path_; };

  void get_dimensions(string filename);
  void get_parameters(string filename);
  void get_flame_position(string filename);
  void add_member() { database.push_back(Data()); };
  int size() { return database.size(); };

  vector<int> filter_by_a(double a);
  vector<int> filter_by_m(double m);
	vector<int> filter_by_LeF(double LeF);
  
  void write_data(std::ofstream &file, vector<int> data_index);
  void remove_data(vector<int> data_index);

 private:
  string path_;
  int n_solutions_ = 0;

  int indexofT1(const Field &vec, const int &NI, const int &NJ);
  int max_Z(const Field &vec, const int &ni, const int &nj);
  string find_parameter_value(const string fname, const string &preString,
                              const string &postString);
  double interpolate_quadratic_maximum(const double &x1, const double &x2,
                                       const double &x3, const double &y1,
                                       const double &y2, const double &y3);
  void add_grid(string filename, string dims);
};

#endif
