#include <cstdio>
#include <x86intrin.h>
#include <random>
#include "conf.hpp"
//----------------------------------------------------------------------
double q[N][D];
double p[N][D] = {};
__attribute__((aligned(64))) double z[N][8];
const double *zp = &(z[0][0]);
int particle_number = 0;
int number_of_pairs = 0;
int number_of_partners[N];
int i_particles[MAX_PAIRS];
int j_particles[MAX_PAIRS];
int pointer[N];
int sorted_list[MAX_PAIRS];
//----------------------------------------------------------------------
class AoSDataManager : public DataManager {
private:
  double (*p)[D];
  double (*q)[D];
public:
  AoSDataManager(double _p[N][D], double _q[N][D]) {
    p = _p;
    q = _q;
  }
  void add_particle(double x, double y, double z, int &particle_number) {
    static std::mt19937 mt(2);
    std::uniform_real_distribution<double> ud(0.0, 0.1);
    q[particle_number][X] = x + ud(mt);
    q[particle_number][Y] = y + ud(mt);
    q[particle_number][Z] = z + ud(mt);
    particle_number++;
  }
  double calc_distance(int i, int j) {
    const double dx = q[i][X] - q[j][X];
    const double dy = q[i][Y] - q[j][Y];
    const double dz = q[i][Z] - q[j][Z];
    return dx * dx + dy * dy + dz * dz;
  }
  void print_results(int particle_number) {
    for (int i = 0; i < 5; i++) {
      printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
    }
    for (int i = particle_number - 5; i < particle_number; i++) {
      printf("%.10f %.10f %.10f\n", p[i][X], p[i][Y], p[i][Z]);
    }
  }
};
//----------------------------------------------------------------------
void
copy_to_z(void) {
  for (int i = 0; i < particle_number; i++) {
    z[i][X] = q[i][X];
    z[i][Y] = q[i][Y];
    z[i][Z] = q[i][Z];
    z[i][PX] = p[i][X];
    z[i][PY] = p[i][Y];
    z[i][PZ] = p[i][Z];
  }
}
//----------------------------------------------------------------------
void
copy_from_z(void) {
  for (int i = 0; i < particle_number; i++) {
    q[i][X] = z[i][X];
    q[i][Y] = z[i][Y];
    q[i][Z] = z[i][Z];
    p[i][X] = z[i][PX];
    p[i][Y] = z[i][PY];
    p[i][Z] = z[i][PZ];
  }
}
//----------------------------------------------------------------------
void
force_pair(void) {
  for (int k = 0; k < number_of_pairs; k++) {
    const int i = i_particles[k];
    const int j = j_particles[k];
    double dx = q[j][X] - q[i][X];
    double dy = q[j][Y] - q[i][Y];
    double dz = q[j][Z] - q[i][Z];
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    p[i][X] += df * dx;
    p[i][Y] += df * dy;
    p[i][Z] += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
force_sorted(void) {
  const int pn = particle_number;
  const int * __restrict sorted_list2 = sorted_list;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list2[kp + k];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      //if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df=0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
void
force_sorted_z(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qix = z[i][X];
    const double qiy = z[i][Y];
    const double qiz = z[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = z[j][X] - qix;
      double dy = z[j][Y] - qiy;
      double dz = z[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      z[j][PX] -= df * dx;
      z[j][PY] -= df * dy;
      z[j][PZ] -= df * dz;
    }
    z[i][PX] += pfx;
    z[i][PY] += pfy;
    z[i][PZ] += pfz;
  }
}
//----------------------------------------------------------------------
void
force_swp(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    int ja = sorted_list[kp];
    double dxa = q[ja][X] - qx_key;
    double dya = q[ja][Y] - qy_key;
    double dza = q[ja][Z] - qz_key;
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
      dxa = q[ja][X] - qx_key;
      dya = q[ja][Y] - qy_key;
      dza = q[ja][Z] - qz_key;
      if (r2 > CL2)continue;
      pfx += df * dxb;
      pfy += df * dyb;
      pfz += df * dzb;
      p[jb][X] -= df * dxb;
      p[jb][Y] -= df * dyb;
      p[jb][Z] -= df * dzb;
      const double r6 = r2 * r2 * r2;
      df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      jb = j;
      dxb = dx;
      dyb = dy;
      dzb = dz;
    }
    p[jb][X] -= df * dxb;
    p[jb][Y] -= df * dyb;
    p[jb][Z] -= df * dzb;
    p[i][X] += pfx + df * dxb;
    p[i][Y] += pfy + df * dyb;
    p[i][Z] += pfz + df * dzb;
  }
}
//----------------------------------------------------------------------
void
force_avx2(void) {
  const v4df vzero = _mm256_set1_pd(0.0);
  const v4df vcl2 = _mm256_set1_pd(CL2);
  const v4df vc24 = _mm256_set1_pd(24 * dt);
  const v4df vc48 = _mm256_set1_pd(48 * dt);
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(q + i));
    v4df vpi = _mm256_load_pd((double*)(p + i));
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_a = sorted_list[kp + k];
      v4df vqj_a = _mm256_load_pd((double*)(q + j_a));
      v4df vdq_a = (vqj_a - vqi);

      const int j_b = sorted_list[kp + k + 1];
      v4df vqj_b = _mm256_load_pd((double*)(q + j_b));
      v4df vdq_b = (vqj_b - vqi);

      const int j_c = sorted_list[kp + k + 2];
      v4df vqj_c = _mm256_load_pd((double*)(q + j_c));
      v4df vdq_c = (vqj_c - vqi);

      const int j_d = sorted_list[kp + k + 3];
      v4df vqj_d = _mm256_load_pd((double*)(q + j_d));
      v4df vdq_d = (vqj_d - vqi);

      v4df tmp0 = _mm256_unpacklo_pd(vdq_a, vdq_b);
      v4df tmp1 = _mm256_unpackhi_pd(vdq_a, vdq_b);
      v4df tmp2 = _mm256_unpacklo_pd(vdq_c, vdq_d);
      v4df tmp3 = _mm256_unpackhi_pd(vdq_c, vdq_d);

      v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      v4df vr6 = vr2 * vr2 * vr2;
      v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      v4df vdf_a = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_b = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_c = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_d = _mm256_permute4x64_pd(vdf, 255);

      v4df vpj_a = _mm256_load_pd((double*)(p + j_a));
      vpi += vdq_a * vdf_a;
      vpj_a -= vdq_a * vdf_a;
      _mm256_store_pd((double*)(p + j_a), vpj_a);

      v4df vpj_b = _mm256_load_pd((double*)(p + j_b));
      vpi += vdq_b * vdf_b;
      vpj_b -= vdq_b * vdf_b;
      _mm256_store_pd((double*)(p + j_b), vpj_b);

      v4df vpj_c = _mm256_load_pd((double*)(p + j_c));
      vpi += vdq_c * vdf_c;
      vpj_c -= vdq_c * vdf_c;
      _mm256_store_pd((double*)(p + j_c), vpj_c);

      v4df vpj_d = _mm256_load_pd((double*)(p + j_d));
      vpi += vdq_d * vdf_d;
      vpj_d -= vdq_d * vdf_d;
      _mm256_store_pd((double*)(p + j_d), vpj_d);
    }
    _mm256_store_pd((double*)(p + i), vpi);
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - q[i][X];
      double dy = q[j][Y] - q[i][Y];
      double dz = q[j][Z] - q[i][Z];
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      p[i][X] += df * dx;
      p[i][Y] += df * dy;
      p[i][Z] += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
  }
}
//----------------------------------------------------------------------
void
force_sorted_z_avx2(void) {
  const v4df vzero = _mm256_set1_pd(0.0);
  const v4df vcl2 = _mm256_set1_pd(CL2);
  const v4df vc24 = _mm256_set1_pd(24 * dt);
  const v4df vc48 = _mm256_set1_pd(48 * dt);
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(z + i));
    v4df vpi = _mm256_load_pd((double*)(zp + i * 8 + 4));
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_a = sorted_list[kp + k];
      v4df vqj_a = _mm256_load_pd((double*)(z + j_a));
      v4df vdq_a = (vqj_a - vqi);

      const int j_b = sorted_list[kp + k + 1];
      v4df vqj_b = _mm256_load_pd((double*)(z + j_b));
      v4df vdq_b = (vqj_b - vqi);

      const int j_c = sorted_list[kp + k + 2];
      v4df vqj_c = _mm256_load_pd((double*)(z + j_c));
      v4df vdq_c = (vqj_c - vqi);

      const int j_d = sorted_list[kp + k + 3];
      v4df vqj_d = _mm256_load_pd((double*)(z + j_d));
      v4df vdq_d = (vqj_d - vqi);

      v4df tmp0 = _mm256_unpacklo_pd(vdq_a, vdq_b);
      v4df tmp1 = _mm256_unpackhi_pd(vdq_a, vdq_b);
      v4df tmp2 = _mm256_unpacklo_pd(vdq_c, vdq_d);
      v4df tmp3 = _mm256_unpackhi_pd(vdq_c, vdq_d);

      v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      v4df vr6 = vr2 * vr2 * vr2;
      v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      v4df vdf_a = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_b = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_c = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_d = _mm256_permute4x64_pd(vdf, 255);

      v4df vpj_a = _mm256_load_pd((double*)(zp + j_a * 8 + 4));
      vpi += vdq_a * vdf_a;
      vpj_a -= vdq_a * vdf_a;
      _mm256_store_pd((double*)(zp + j_a * 8 + 4), vpj_a);

      v4df vpj_b = _mm256_load_pd((double*)(zp + j_b * 8 + 4));
      vpi += vdq_b * vdf_b;
      vpj_b -= vdq_b * vdf_b;
      _mm256_store_pd((double*)(zp + j_b * 8 + 4), vpj_b);

      v4df vpj_c = _mm256_load_pd((double*)(zp + j_c * 8 + 4));
      vpi += vdq_c * vdf_c;
      vpj_c -= vdq_c * vdf_c;
      _mm256_store_pd((double*)(zp + j_c * 8 + 4), vpj_c);

      v4df vpj_d = _mm256_load_pd((double*)(zp + j_d * 8 + 4));
      vpi += vdq_d * vdf_d;
      vpj_d -= vdq_d * vdf_d;
      _mm256_store_pd((double*)(zp + j_d * 8 + 4), vpj_d);
    }
    _mm256_store_pd((double*)(zp + i * 8 + 4), vpi);
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = z[j][X] - z[i][X];
      double dy = z[j][Y] - z[i][Y];
      double dz = z[j][Z] - z[i][Z];
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      z[i][PX] += df * dx;
      z[i][PY] += df * dy;
      z[i][PZ] += df * dz;
      z[j][PX] -= df * dx;
      z[j][PY] -= df * dy;
      z[j][PZ] -= df * dz;
    }
  }
}
//----------------------------------------------------------------------
int
main(void) {
  AoSDataManager aosdm(p, q);
  init(&aosdm, particle_number);
  check_pairlist(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, &aosdm);
  sortpair(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, pointer, sorted_list);
#ifdef PAIR
  measure(&force_pair, "pair", particle_number);
  aosdm.print_results(particle_number);
#elif AVX2
  measure(&force_avx2, "avx2", particle_number);
  aosdm.print_results(particle_number);
#elif SORTED_Z
  copy_to_z();
  measure(&force_sorted_z, "sorted_z", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#elif SORTED_Z_AVX2
  copy_to_z();
  measure(&force_sorted_z_avx2, "sorted_z_avx2", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#else
  measure(&force_pair, "pair", particle_number);
  measure(&force_sorted, "sorted", particle_number);
  measure(&force_swp, "sorted_swp", particle_number);
  measure(&force_avx2, "avx2", particle_number);
  measure(&force_sorted_z, "sorted_z", particle_number);
  measure(&force_sorted_z_avx2, "sorted_z_avx2", particle_number);
#endif
}
//----------------------------------------------------------------------
