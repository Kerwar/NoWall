#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <algorithm>

namespace parameters {
// NUMERICAL
constexpr int NPROCS = 64;

constexpr int NINPUT = 1200;
constexpr int MINPUT = 80;

constexpr double TOL = 10E-8;

// DOMAIN
constexpr double channel_xmin = 0;
constexpr double channel_xmax = 120;

constexpr double wall_xmin = 40;
constexpr double wall_xmax = 80;

constexpr double y_bot_min = 0.5;
constexpr double y_bot_max = 1.;

constexpr double y_wall_min = 1.;
constexpr double y_wall_max = 1.;

constexpr double y_top_min = 1.0;
constexpr double y_top_max = 1.5;

constexpr double a = 1;
constexpr double a2 = 1 / (a * a);

constexpr double BK = 0.15;

// REACTION
// static double q = 1.2;
constexpr double beta_reaction = 10;
constexpr double gamma_reaction = 0.7;

constexpr double LeF = 1.0;
constexpr double LeZ = 0.3;
constexpr double xHotSpot = 46;
constexpr double z0hs = 0.5;
constexpr double r0hs = 1;
// FROM HERE NO NEED TO CHANGE THEM,
// THEY ARE PARAMETERS THAT COME FROM THE OTHERS

constexpr int NPROCSINROW = NPROCS / 2;
constexpr int NPROCSINCOL = 1;

constexpr int NTOTAL = NINPUT + (NPROCS + (-NINPUT) % NPROCS) % NPROCS;
constexpr int MTOTAL = MINPUT;

constexpr int N = NTOTAL / (NPROCS / 2);
constexpr int M = MINPUT / 2;

constexpr int NI = N + 2;
constexpr int NJ = M + 2;

constexpr double DXX = (channel_xmax - channel_xmin) / NTOTAL;
constexpr double DYY = (y_top_max - y_bot_min) / MTOTAL;

constexpr int calc_wall_i_min() {
  int index = 0;
  for (int i = 0; i < NTOTAL + 2; i++)
    if (std::max(channel_xmin, channel_xmin + DXX * i - DXX / 2) >= wall_xmin) {
      index = i;
      break;
    }
  return index;
}

constexpr int calc_wall_i_max() {
  int index = 0;
  for (int i = 0; i < NTOTAL + 2; i++)
    if (std::max(channel_xmin, channel_xmin + DXX * i - DXX / 2) > wall_xmax) {
      index = i;
      break;
    }
  return index;
}

constexpr int IWALLMIN = calc_wall_i_min();
constexpr int IWALLMAX = calc_wall_i_max();

constexpr int NTOTALWALL = IWALLMAX - IWALLMIN;
constexpr int NWALL = NTOTALWALL / NPROCS;

constexpr double help_EXCTE = BK * a * a / 2.0;
constexpr double EXCTE = help_EXCTE / (1.0 + DYY * help_EXCTE);
}  // namespace parameters
#endif