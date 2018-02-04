#include <cstdio>
#include <x86intrin.h>
#include <random>
#include "conf.hpp"
//----------------------------------------------------------------------
double q[D][N];
double p[D][N] = {};
int particle_number = 0;
int number_of_pairs = 0;
int number_of_partners[N];
int i_particles[MAX_PAIRS];
int j_particles[MAX_PAIRS];
int pointer[N], pointer2[N];
int sorted_list[MAX_PAIRS];
//----------------------------------------------------------------------
class SoADataManager : public DataManager {
private:
  double (*p)[N];
  double (*q)[N];
public:
  SoADataManager(double _p[D][N], double _q[D][N]) {
    p = _p;
    q = _q;
  }
  void add_particle(double x, double y, double z, int &particle_number) {
    static std::mt19937 mt(2);
    std::uniform_real_distribution<double> ud(0.0, 0.1);
    q[X][particle_number] = x + ud(mt);
    q[Y][particle_number] = y + ud(mt);
    q[Z][particle_number] = z + ud(mt);
    particle_number++;
  }
  double calc_distance(int i, int j) {
    const double dx = q[X][i] - q[X][j];
    const double dy = q[Y][i] - q[Y][j];
    const double dz = q[Z][i] - q[Z][j];
    return dx * dx + dy * dy + dz * dz;
  }

  void print_results(int particle_number) {
    for (int i = 0; i < 5; i++) {
      printf("%.10f %.10f %.10f\n", p[X][i], p[Y][i], p[Z][i]);
    }
    for (int i = particle_number - 5; i < particle_number; i++) {
      printf("%.10f %.10f %.10f\n", p[X][i], p[Y][i], p[Z][i]);
    }
  }
};
//----------------------------------------------------------------------
void
force_pair(void) {
  for (int k = 0; k < number_of_pairs; k++) {
    const int i = i_particles[k];
    const int j = j_particles[k];
    double dx = q[X][j] - q[X][i];
    double dy = q[Y][j] - q[Y][i];
    double dz = q[Z][j] - q[Z][i];
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    p[X][i] += df * dx;
    p[Y][i] += df * dy;
    p[Z][i] += df * dz;
    p[X][j] -= df * dx;
    p[Y][j] -= df * dy;
    p[Z][j] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
force_sorted(void) {
  const int pn = particle_number;
  const int * __restrict sorted_list2 = sorted_list;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[X][i];
    const double qy_key = q[Y][i];
    const double qz_key = q[Z][i];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list2[kp + k];
      double dx = q[X][j] - qx_key;
      double dy = q[Y][j] - qy_key;
      double dz = q[Z][j] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[X][j] -= df * dx;
      p[Y][j] -= df * dy;
      p[Z][j] -= df * dz;
    }
    p[X][i] += pfx;
    p[Y][i] += pfy;
    p[Z][i] += pfz;
  }
}
//----------------------------------------------------------------------
void
force_avx2(void) {
  const int pn = particle_number;
  const v4df vzero = _mm256_set_pd(0, 0, 0, 0);
  const v4df vcl2 = _mm256_set_pd(CL2, CL2, CL2, CL2);
  const v4df vc24 = _mm256_set_pd(24 * dt, 24 * dt, 24 * dt, 24 * dt);
  const v4df vc48 = _mm256_set_pd(48 * dt, 48 * dt, 48 * dt, 48 * dt);
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[X][i];
    const double qy_key = q[Y][i];
    const double qz_key = q[Z][i];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    const v4df vqix = _mm256_set_pd(q[X][i], q[X][i], q[X][i], q[X][i]);
    const v4df vqiy = _mm256_set_pd(q[Y][i], q[Y][i], q[Y][i], q[Y][i]);
    const v4df vqiz = _mm256_set_pd(q[Z][i], q[Z][i], q[Z][i], q[Z][i]);

    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_a = sorted_list[kp + k];
      const int j_b = sorted_list[kp + k + 1];
      const int j_c = sorted_list[kp + k + 2];
      const int j_d = sorted_list[kp + k + 3];

      const v4df vqjx = _mm256_set_pd(q[X][j_a], q[X][j_b], q[X][j_c], q[X][j_d]);
      const v4df vqjy = _mm256_set_pd(q[Y][j_a], q[Y][j_b], q[Y][j_c], q[Y][j_d]);
      const v4df vqjz = _mm256_set_pd(q[Z][j_a], q[Z][j_b], q[Z][j_c], q[Z][j_d]);

      const v4df vdx = vqjx - vqix;
      const v4df vdy = vqjy - vqiy;
      const v4df vdz = vqjz - vqiz;
      const v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      const v4df vr6 = vr2 * vr2 * vr2;
      v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      const v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);
      v4df vdfx = vdf * vdx;
      v4df vdfy = vdf * vdy;
      v4df vdfz = vdf * vdz;
      double *dfx = (double*)(&vdfx);
      double *dfy = (double*)(&vdfy);
      double *dfz = (double*)(&vdfz);

      pfx += dfx[3];
      pfy += dfy[3];
      pfz += dfz[3];

      p[X][j_a] -= dfx[3];
      p[Y][j_a] -= dfy[3];
      p[Z][j_a] -= dfz[3];

      pfx += dfx[2];
      pfy += dfy[2];
      pfz += dfz[2];
      p[X][j_b] -= dfx[2];
      p[Y][j_b] -= dfy[2];
      p[Z][j_b] -= dfz[2];

      pfx += dfx[1];
      pfy += dfy[1];
      pfz += dfz[1];
      p[X][j_c] -= dfx[1];
      p[Y][j_c] -= dfy[1];
      p[Z][j_c] -= dfz[1];

      pfx += dfx[0];
      pfy += dfy[0];
      pfz += dfz[0];
      p[X][j_d] -= dfx[0];
      p[Y][j_d] -= dfy[0];
      p[Z][j_d] -= dfz[0];
    }
    p[X][i] += pfx;
    p[Y][i] += pfy;
    p[Z][i] += pfz;
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[X][j] - qx_key;
      double dy = q[Y][j] - qy_key;
      double dz = q[Z][j] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      p[X][i] += df * dx;
      p[Y][i] += df * dy;
      p[Z][i] += df * dz;
      p[X][j] -= df * dx;
      p[Y][j] -= df * dy;
      p[Z][j] -= df * dz;
    }
  }
}
//----------------------------------------------------------------------
void
force_swp(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[X][i];
    const double qy_key = q[Y][i];
    const double qz_key = q[Z][i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    int ja = sorted_list[kp];
    double dxa = q[X][ja] - qx_key;
    double dya = q[Y][ja] - qy_key;
    double dza = q[Z][ja] - qz_key;
    double df = 0.0;
    double dxb = 0.0, dyb = 0.0, dzb = 0.0;
    int jb = 0;
    const int np = number_of_partners[i];
    for (int k = kp; k < np + kp; k++) {
      const double dx = dxa;
      const double dy = dya;
      const double dz = dza;
      double r2 = (dx * dx + dy * dy + dz * dz);
      const int j = ja;
      ja = sorted_list[k + 1];
      dxa = q[X][ja] - qx_key;
      dya = q[Y][ja] - qy_key;
      dza = q[Z][ja] - qz_key;
      if (r2 > CL2)continue;
      pfx += df * dxb;
      pfy += df * dyb;
      pfz += df * dzb;
      p[X][jb] -= df * dxb;
      p[Y][jb] -= df * dyb;
      p[Z][jb] -= df * dzb;
      const double r6 = r2 * r2 * r2;
      df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      jb = j;
      dxb = dx;
      dyb = dy;
      dzb = dz;
    }
    p[X][jb] -= df * dxb;
    p[Y][jb] -= df * dyb;
    p[Z][jb] -= df * dzb;
    p[X][i] += pfx + df * dxb;
    p[Y][i] += pfy + df * dyb;
    p[Z][i] += pfz + df * dzb;
  }
}
//----------------------------------------------------------------------
int
main(void) {
  SoADataManager soadm(p, q);
  init(&soadm, particle_number);
  check_pairlist(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, &soadm);
  sortpair(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, pointer, sorted_list);
#ifdef PAIR
  measure(&force_pair, "pair", particle_number);
  soadm.print_results(particle_number);
#elif SORTED
  measure(&force_sorted, "sorted", particle_number);
  soadm.print_results(particle_number);
#elif AVX2
  measure(&force_avx2, "avx2", particle_number);
  soadm.print_results(particle_number);
#else
  measure(&force_pair, "pair", particle_number);
  measure(&force_sorted, "sorted", particle_number);
  measure(&force_swp, "sorted_swp", particle_number);
  measure(&force_avx2, "avx2", particle_number);
#endif
}
//----------------------------------------------------------------------
