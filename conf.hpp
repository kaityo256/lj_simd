#pragma once
//----------------------------------------------------------------------
#include <math.h>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>
#include <chrono>
#include <fstream>
//----------------------------------------------------------------------
const double density = 1.0;
const int N = 400000;
const int MAX_PAIRS = 30 * N;
double L = 50.0;
const double dt = 0.001;
enum {X, Y, Z, W, PX, PY, PZ, PW};
const double CUTOFF_LENGTH = 3.0;
const double SEARCH_LENGTH = 3.3;
const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
//----------------------------------------------------------------------
typedef double v4df __attribute__((vector_size(32)));
//----------------------------------------------------------------------
void
printv(v4df r) {
  double *a = (double*)(&r);
  printf("%.10f %.10f %.10f %.10f\n", a[0], a[1], a[2], a[3]);
}
//----------------------------------------------------------------------
// インテルコンパイラのループ交換最適化阻害のためのダミー変数
int sum = 0;
//----------------------------------------------------------------------
void
measure(void(*pfunc)(), const char *name, int particle_number) {
  const auto s = std::chrono::system_clock::now();
  const int LOOP = 100;
  for (int i = 0; i < LOOP; i++) {
    sum++; // ループ交換阻害
    pfunc();
  }
  const auto e = std::chrono::system_clock::now();
  const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();
  fprintf(stderr, "N=%d, %s %lld [ms]\n", particle_number, name, elapsed);
}
//----------------------------------------------------------------------
typedef void (* fp)(void);
//----------------------------------------------------------------------
void
loadpair(int &number_of_pairs, int number_of_partners[N], int i_particles[MAX_PAIRS], int j_particles[MAX_PAIRS]) {
  std::ifstream ifs("pair.dat", std::ios::binary);
  ifs.read((char*)&number_of_pairs, sizeof(int));
  ifs.read((char*)number_of_partners, sizeof(int)*N);
  ifs.read((char*)i_particles, sizeof(int)*MAX_PAIRS);
  ifs.read((char*)j_particles, sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
savepair(int number_of_pairs, int number_of_partners[N], int i_particles[MAX_PAIRS], int j_particles[MAX_PAIRS]) {
  std::ofstream ofs("pair.dat", std::ios::binary);
  ofs.write((char*)&number_of_pairs, sizeof(int));
  ofs.write((char*)number_of_partners, sizeof(int)*N);
  ofs.write((char*)i_particles, sizeof(int)*MAX_PAIRS);
  ofs.write((char*)j_particles, sizeof(int)*MAX_PAIRS);
}
//----------------------------------------------------------------------
void
register_pair(int index1, int index2, int &number_of_pairs, int number_of_partners[N], int i_particles[MAX_PAIRS], int j_particles[MAX_PAIRS]) {
  int i, j;
  if (index1 < index2) {
    i = index1;
    j = index2;
  } else {
    i = index2;
    j = index1;
  }
  i_particles[number_of_pairs] = i;
  j_particles[number_of_pairs] = j;
  number_of_partners[i]++;
  number_of_pairs++;
}
//----------------------------------------------------------------------
class DataManager {
public:
  virtual void add_particle(double x, double y, double z, int &particle_number) = 0;
  virtual double calc_distance(int i, int j) = 0;
};
//----------------------------------------------------------------------
void
makepair(int particle_number, int &number_of_pairs, int number_of_partners[N], int i_particles[MAX_PAIRS], int j_particles[MAX_PAIRS], DataManager *dm) {
  const double SL2 = SEARCH_LENGTH * SEARCH_LENGTH;
  for (int i = 0; i < particle_number; i++) {
    number_of_partners[i] = 0;
  }
  for (int i = 0; i < particle_number - 1; i++) {
    for (int j = i + 1; j < particle_number; j++) {
      /*
      const double dx = q[i][X] - q[j][X];
      const double dy = q[i][Y] - q[j][Y];
      const double dz = q[i][Z] - q[j][Z];
      const double r2 = dx * dx + dy * dy + dz * dz;
      */
      const double r2 = dm->calc_distance(i, j);
      if (r2 < SL2) {
        register_pair(i, j, number_of_pairs, number_of_partners, i_particles, j_particles);
      }
    }
  }
}
//----------------------------------------------------------------------
void
check_pairlist(int particle_number, int &number_of_pairs, int number_of_partners[N], int i_particles[MAX_PAIRS], int j_particles[MAX_PAIRS], DataManager *dm) {
  struct stat st;
  int ret = stat("pair.dat", &st);
  if (ret == 0) {
    std::cerr << "A pair-file is found. I use it." << std::endl;
    loadpair(number_of_pairs, number_of_partners, i_particles, j_particles);
  } else {
    std::cerr << "Make pairlist." << std::endl;
    //makepair();
    makepair(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, dm);
    savepair(number_of_pairs, number_of_partners, i_particles, j_particles);
  }
  std::cerr << "Number of pairs: " << number_of_pairs << std::endl;
}
//----------------------------------------------------------------------
void
sortpair(const int pn, int &number_of_pairs, int number_of_partners[N], int i_particles[MAX_PAIRS], int j_particles[MAX_PAIRS], int pointer[N], int sorted_list[MAX_PAIRS]) {
  static int pointer2[N];
  int pos = 0;
  pointer[0] = 0;
  for (int i = 0; i < pn - 1; i++) {
    pos += number_of_partners[i];
    pointer[i + 1] = pos;
  }
  for (int i = 0; i < pn; i++) {
    pointer2[i] = 0;
  }
  const int s = number_of_pairs;
  for (int k = 0; k < s; k++) {
    int i = i_particles[k];
    int j = j_particles[k];
    int index = pointer[i] + pointer2[i];
    sorted_list[index] = j;
    pointer2[i] ++;
  }
}
//----------------------------------------------------------------------
void
init(DataManager *dm, int &particle_number) {
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L / s);
  int sy = static_cast<int>(L / s);
  int sz = static_cast<int>(L / s);
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        double x = ix * s;
        double y = iy * s;
        double z = iz * s;
        dm->add_particle(x, y, z, particle_number);
        dm->add_particle(x, y + hs, z + hs, particle_number);
        dm->add_particle(x + hs, y, z + hs, particle_number);
        dm->add_particle(x + hs, y + hs, z, particle_number);
      }
    }
  }
}
//----------------------------------------------------------------------
