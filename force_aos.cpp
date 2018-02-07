#include <cstdio>
#include <x86intrin.h>
#include <random>
#include "conf.hpp"
//----------------------------------------------------------------------
double q[N][D];
double p[N][D] = {};
__attribute__((aligned(64))) double z[N][8];
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
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list2[kp + k];
      double dx = q[j][X] - qix;
      double dy = q[j][Y] - qiy;
      double dz = q[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      //if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
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
      //if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
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
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    int ja = sorted_list[kp];
    double dxa = q[ja][X] - qix;
    double dya = q[ja][Y] - qiy;
    double dza = q[ja][Z] - qiz;
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
      dxa = q[ja][X] - qix;
      dya = q[ja][Y] - qiy;
      dza = q[ja][Z] - qiz;
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
      const int j_1 = sorted_list[kp + k];
      const int j_2 = sorted_list[kp + k + 1];
      const int j_3 = sorted_list[kp + k + 2];
      const int j_4 = sorted_list[kp + k + 3];
      v4df vqj_1 = _mm256_load_pd((double*)(q + j_1));
      v4df vqj_2 = _mm256_load_pd((double*)(q + j_2));
      v4df vqj_3 = _mm256_load_pd((double*)(q + j_3));
      v4df vqj_4 = _mm256_load_pd((double*)(q + j_4));

      v4df vdq_1 = (vqj_1 - vqi);
      v4df vdq_2 = (vqj_2 - vqi);
      v4df vdq_3 = (vqj_3 - vqi);
      v4df vdq_4 = (vqj_4 - vqi);

      v4df tmp0 = _mm256_unpacklo_pd(vdq_1, vdq_2);
      v4df tmp1 = _mm256_unpackhi_pd(vdq_1, vdq_2);
      v4df tmp2 = _mm256_unpacklo_pd(vdq_3, vdq_4);
      v4df tmp3 = _mm256_unpackhi_pd(vdq_3, vdq_4);

      v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      v4df vr6 = vr2 * vr2 * vr2;
      v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

      v4df vpj_1 = _mm256_load_pd((double*)(p + j_1));
      v4df vpj_2 = _mm256_load_pd((double*)(p + j_2));
      v4df vpj_3 = _mm256_load_pd((double*)(p + j_3));
      v4df vpj_4 = _mm256_load_pd((double*)(p + j_4));

      vpj_1 -= vdq_1 * vdf_1;
      vpj_2 -= vdq_2 * vdf_2;
      vpj_3 -= vdq_3 * vdf_3;
      vpj_4 -= vdq_4 * vdf_4;

      vpi += vdq_1 * vdf_1;
      vpi += vdq_2 * vdf_2;
      vpi += vdq_3 * vdf_3;
      vpi += vdq_4 * vdf_4;

      _mm256_store_pd((double*)(p + j_1), vpj_1);
      _mm256_store_pd((double*)(p + j_2), vpj_2);
      _mm256_store_pd((double*)(p + j_3), vpj_3);
      _mm256_store_pd((double*)(p + j_4), vpj_4);
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
force_avx2_swp(void) {
  const v4df vzero = _mm256_set1_pd(0.0);
  const v4df vcl2 = _mm256_set1_pd(CL2);
  const v4df vc24 = _mm256_set1_pd(24 * dt);
  const v4df vc48 = _mm256_set1_pd(48 * dt);
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(&q[i][X]));
    v4df vpi = _mm256_load_pd((double*)(&p[i][X]));
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    int k = 0;
    // --- 8< ---
    int j_1 = sorted_list[kp + k];
    int j_2 = sorted_list[kp + k + 1];
    int j_3 = sorted_list[kp + k + 2];
    int j_4 = sorted_list[kp + k + 3];
    v4df vqj_1 = _mm256_load_pd((double*)(&q[j_1][X]));
    v4df vqj_2 = _mm256_load_pd((double*)(&q[j_2][X]));
    v4df vqj_3 = _mm256_load_pd((double*)(&q[j_3][X]));
    v4df vqj_4 = _mm256_load_pd((double*)(&q[j_4][X]));

    v4df vdq_1 = (vqj_1 - vqi);
    v4df vdq_2 = (vqj_2 - vqi);
    v4df vdq_3 = (vqj_3 - vqi);
    v4df vdq_4 = (vqj_4 - vqi);

    v4df tmp0 = _mm256_unpacklo_pd(vdq_1, vdq_2);
    v4df tmp1 = _mm256_unpackhi_pd(vdq_1, vdq_2);
    v4df tmp2 = _mm256_unpacklo_pd(vdq_3, vdq_4);
    v4df tmp3 = _mm256_unpackhi_pd(vdq_3, vdq_4);

    v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
    v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
    v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

    v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
    v4df vr6 = vr2 * vr2 * vr2;
    v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
    v4df mask = vcl2 - vr2;
    vdf = _mm256_blendv_pd(vdf, vzero, mask);

    if(np<4){
      vdf =  _mm256_setzero_pd();
    }

    v4df vdf_1, vdf_2, vdf_3, vdf_4;
    v4df vpj_1, vpj_2, vpj_3, vpj_4;
    for (k = 4; k < (np / 4) * 4; k += 4) {
      // --- 8< ---
      const int j_1_b = sorted_list[kp + k];
      const int j_2_b = sorted_list[kp + k + 1];
      const int j_3_b = sorted_list[kp + k + 2];
      const int j_4_b = sorted_list[kp + k + 3];
      vqj_1 = _mm256_load_pd((double*)(&q[j_1_b][X]));
      vqj_2 = _mm256_load_pd((double*)(&q[j_2_b][X]));
      vqj_3 = _mm256_load_pd((double*)(&q[j_3_b][X]));
      vqj_4 = _mm256_load_pd((double*)(&q[j_4_b][X]));
      v4df vdq_1_b = (vqj_1 - vqi);
      v4df vdq_2_b = (vqj_2 - vqi);
      v4df vdq_3_b = (vqj_3 - vqi);
      v4df vdq_4_b = (vqj_4 - vqi);
      tmp0 = _mm256_unpacklo_pd(vdq_1_b, vdq_2_b);
      tmp1 = _mm256_unpackhi_pd(vdq_1_b, vdq_2_b);
      tmp2 = _mm256_unpacklo_pd(vdq_3_b, vdq_4_b);
      tmp3 = _mm256_unpackhi_pd(vdq_3_b, vdq_4_b);
      vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      vr6 = vr2 * vr2 * vr2;

      // --- 8< ---
      vdf_1 = _mm256_permute4x64_pd(vdf, 0);
      vdf_2 = _mm256_permute4x64_pd(vdf, 85);
      vdf_3 = _mm256_permute4x64_pd(vdf, 170);
      vdf_4 = _mm256_permute4x64_pd(vdf, 255);

      vpj_1 = _mm256_load_pd((double*)(&p[j_1][X]));
      vpj_2 = _mm256_load_pd((double*)(&p[j_2][X]));
      vpj_3 = _mm256_load_pd((double*)(&p[j_3][X]));
      vpj_4 = _mm256_load_pd((double*)(&p[j_4][X]));

      vpj_1 -= vdq_1 * vdf_1;
      vpj_2 -= vdq_2 * vdf_2;
      vpj_3 -= vdq_3 * vdf_3;
      vpj_4 -= vdq_4 * vdf_4;

      vpi += vdq_1 * vdf_1;
      vpi += vdq_4 * vdf_4;
      vpi += vdq_2 * vdf_2;
      vpi += vdq_3 * vdf_3;

      _mm256_store_pd((double*)(&p[j_1][X]), vpj_1);
      _mm256_store_pd((double*)(&p[j_2][X]), vpj_2);
      _mm256_store_pd((double*)(&p[j_3][X]), vpj_3);
      _mm256_store_pd((double*)(&p[j_4][X]), vpj_4);

      // --- 8< ---
      j_1 = j_1_b; 
      j_2 = j_2_b; 
      j_3 = j_3_b; 
      j_4 = j_4_b; 
      vdq_1 = vdq_1_b;
      vdq_2 = vdq_2_b;
      vdq_3 = vdq_3_b;
      vdq_4 = vdq_4_b;
      vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);
    }
    // --- 8< ---
    vdf_1 = _mm256_permute4x64_pd(vdf, 0);
    vdf_2 = _mm256_permute4x64_pd(vdf, 85);
    vdf_3 = _mm256_permute4x64_pd(vdf, 170);
    vdf_4 = _mm256_permute4x64_pd(vdf, 255);

    vpj_1 = _mm256_load_pd((double*)(&p[j_1][X]));
    vpj_2 = _mm256_load_pd((double*)(&p[j_2][X]));
    vpj_3 = _mm256_load_pd((double*)(&p[j_3][X]));
    vpj_4 = _mm256_load_pd((double*)(&p[j_4][X]));

    vpi += vdq_1 * vdf_1;
    vpi += vdq_2 * vdf_2;
    vpi += vdq_3 * vdf_3;
    vpi += vdq_4 * vdf_4;
    vpj_1 -= vdq_1 * vdf_1;
    vpj_2 -= vdq_2 * vdf_2;
    vpj_3 -= vdq_3 * vdf_3;
    vpj_4 -= vdq_4 * vdf_4;

    _mm256_store_pd((double*)(&p[j_1][X]), vpj_1);
    _mm256_store_pd((double*)(&p[j_2][X]), vpj_2);
    _mm256_store_pd((double*)(&p[j_3][X]), vpj_3);
    _mm256_store_pd((double*)(&p[j_4][X]), vpj_4);

    _mm256_store_pd((double*)(&p[i][X]), vpi);

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
//------------------------------------------------------------------------
void
force_sorted_z_avx2(void) {
  const v4df vzero = _mm256_set1_pd(0.0);
  const v4df vcl2 = _mm256_set1_pd(CL2);
  const v4df vc24 = _mm256_set1_pd(24 * dt);
  const v4df vc48 = _mm256_set1_pd(48 * dt);
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(z + i));
    v4df vpi = _mm256_load_pd((double*)(&z[i][PX]));
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_1 = sorted_list[kp + k];
      const int j_2 = sorted_list[kp + k + 1];
      const int j_3 = sorted_list[kp + k + 2];
      const int j_4 = sorted_list[kp + k + 3];
      v4df vqj_1 = _mm256_load_pd((double*)(&z[j_1][X]));
      v4df vqj_2 = _mm256_load_pd((double*)(&z[j_2][X]));
      v4df vqj_3 = _mm256_load_pd((double*)(&z[j_3][X]));
      v4df vqj_4 = _mm256_load_pd((double*)(&z[j_4][X]));

      v4df vdq_1 = (vqj_1 - vqi);
      v4df vdq_2 = (vqj_2 - vqi);
      v4df vdq_3 = (vqj_3 - vqi);
      v4df vdq_4 = (vqj_4 - vqi);

      v4df tmp0 = _mm256_unpacklo_pd(vdq_1, vdq_2);
      v4df tmp1 = _mm256_unpackhi_pd(vdq_1, vdq_2);
      v4df tmp2 = _mm256_unpacklo_pd(vdq_3, vdq_4);
      v4df tmp3 = _mm256_unpackhi_pd(vdq_3, vdq_4);

      v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      v4df vr6 = vr2 * vr2 * vr2;
      v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

      v4df vpj_1 = _mm256_load_pd((double*)(&z[j_1][PX]));
      v4df vpj_2 = _mm256_load_pd((double*)(&z[j_2][PX]));
      v4df vpj_3 = _mm256_load_pd((double*)(&z[j_3][PX]));
      v4df vpj_4 = _mm256_load_pd((double*)(&z[j_4][PX]));

      vpi += vdq_1 * vdf_1;
      vpi += vdq_2 * vdf_2;
      vpi += vdq_3 * vdf_3;
      vpi += vdq_4 * vdf_4;

      vpj_1 -= vdq_1 * vdf_1;
      vpj_2 -= vdq_2 * vdf_2;
      vpj_3 -= vdq_3 * vdf_3;
      vpj_4 -= vdq_4 * vdf_4;

      _mm256_store_pd((double*)(&z[j_1][PX]), vpj_1);
      _mm256_store_pd((double*)(&z[j_2][PX]), vpj_2);
      _mm256_store_pd((double*)(&z[j_3][PX]), vpj_3);
      _mm256_store_pd((double*)(&z[j_4][PX]), vpj_4);
    }
    _mm256_store_pd((double*)(&z[i][PX]), vpi);
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = z[j][X] - z[i][X];
      double dy = z[j][Y] - z[i][Y];
      double dz = z[j][Z] - z[i][Z];
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df=0.0;
      z[i][PX] += df * dx;
      z[i][PY] += df * dy;
      z[i][PZ] += df * dz;
      z[j][PX] -= df * dx;
      z[j][PY] -= df * dy;
      z[j][PZ] -= df * dz;
    }
  }
}
//------------------------------------------------------------------------
void
force_sorted_z_avx2_swp(void) {
  const v4df vzero = _mm256_set1_pd(0.0);
  const v4df vcl2 = _mm256_set1_pd(CL2);
  const v4df vc24 = _mm256_set1_pd(24 * dt);
  const v4df vc48 = _mm256_set1_pd(48 * dt);
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const v4df vqi = _mm256_load_pd((double*)(&z[i][X]));
    v4df vpi = _mm256_load_pd((double*)(&z[i][PX]));
    const int np = number_of_partners[i];
    const int kp = pointer[i];
    int k = 0;
    // --- 8< ---
    int j_1 = sorted_list[kp + k];
    int j_2 = sorted_list[kp + k + 1];
    int j_3 = sorted_list[kp + k + 2];
    int j_4 = sorted_list[kp + k + 3];
    v4df vqj_1 = _mm256_load_pd((double*)(&z[j_1][X]));
    v4df vqj_2 = _mm256_load_pd((double*)(&z[j_2][X]));
    v4df vqj_3 = _mm256_load_pd((double*)(&z[j_3][X]));
    v4df vqj_4 = _mm256_load_pd((double*)(&z[j_4][X]));

    v4df vdq_1 = (vqj_1 - vqi);
    v4df vdq_2 = (vqj_2 - vqi);
    v4df vdq_3 = (vqj_3 - vqi);
    v4df vdq_4 = (vqj_4 - vqi);

    v4df tmp0 = _mm256_unpacklo_pd(vdq_1, vdq_2);
    v4df tmp1 = _mm256_unpackhi_pd(vdq_1, vdq_2);
    v4df tmp2 = _mm256_unpacklo_pd(vdq_3, vdq_4);
    v4df tmp3 = _mm256_unpackhi_pd(vdq_3, vdq_4);

    v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
    v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
    v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

    v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
    v4df vr6 = vr2 * vr2 * vr2;
    v4df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
    v4df mask = vcl2 - vr2;
    vdf = _mm256_blendv_pd(vdf, vzero, mask);

    if(np<4){
      vdf =  _mm256_setzero_pd();
    }

    v4df vdf_1, vdf_2, vdf_3, vdf_4;
    v4df vpj_1, vpj_2, vpj_3, vpj_4;
    for (k = 4; k < (np / 4) * 4; k += 4) {
      // --- 8< ---
      const int j_1_b = sorted_list[kp + k];
      const int j_2_b = sorted_list[kp + k + 1];
      const int j_3_b = sorted_list[kp + k + 2];
      const int j_4_b = sorted_list[kp + k + 3];
      vqj_1 = _mm256_load_pd((double*)(&z[j_1_b][X]));
      vqj_2 = _mm256_load_pd((double*)(&z[j_2_b][X]));
      vqj_3 = _mm256_load_pd((double*)(&z[j_3_b][X]));
      vqj_4 = _mm256_load_pd((double*)(&z[j_4_b][X]));
      v4df vdq_1_b = (vqj_1 - vqi);
      v4df vdq_2_b = (vqj_2 - vqi);
      v4df vdq_3_b = (vqj_3 - vqi);
      v4df vdq_4_b = (vqj_4 - vqi);
      tmp0 = _mm256_unpacklo_pd(vdq_1_b, vdq_2_b);
      tmp1 = _mm256_unpackhi_pd(vdq_1_b, vdq_2_b);
      tmp2 = _mm256_unpacklo_pd(vdq_3_b, vdq_4_b);
      tmp3 = _mm256_unpackhi_pd(vdq_3_b, vdq_4_b);
      vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      vr6 = vr2 * vr2 * vr2;

      // --- 8< ---
      vdf_1 = _mm256_permute4x64_pd(vdf, 0);
      vdf_2 = _mm256_permute4x64_pd(vdf, 85);
      vdf_3 = _mm256_permute4x64_pd(vdf, 170);
      vdf_4 = _mm256_permute4x64_pd(vdf, 255);

      vpj_1 = _mm256_load_pd((double*)(&z[j_1][PX]));
      vpj_2 = _mm256_load_pd((double*)(&z[j_2][PX]));
      vpj_3 = _mm256_load_pd((double*)(&z[j_3][PX]));
      vpj_4 = _mm256_load_pd((double*)(&z[j_4][PX]));

      vpj_1 -= vdq_1 * vdf_1;
      vpj_2 -= vdq_2 * vdf_2;
      vpj_3 -= vdq_3 * vdf_3;
      vpj_4 -= vdq_4 * vdf_4;

      vpi += vdq_1 * vdf_1;
      vpi += vdq_4 * vdf_4;
      vpi += vdq_2 * vdf_2;
      vpi += vdq_3 * vdf_3;

      _mm256_store_pd((double*)(&z[j_1][PX]), vpj_1);
      _mm256_store_pd((double*)(&z[j_2][PX]), vpj_2);
      _mm256_store_pd((double*)(&z[j_3][PX]), vpj_3);
      _mm256_store_pd((double*)(&z[j_4][PX]), vpj_4);

      // --- 8< ---
      j_1 = j_1_b; 
      j_2 = j_2_b; 
      j_3 = j_3_b; 
      j_4 = j_4_b; 
      vdq_1 = vdq_1_b;
      vdq_2 = vdq_2_b;
      vdq_3 = vdq_3_b;
      vdq_4 = vdq_4_b;
      vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);
    }
    // --- 8< ---
    vdf_1 = _mm256_permute4x64_pd(vdf, 0);
    vdf_2 = _mm256_permute4x64_pd(vdf, 85);
    vdf_3 = _mm256_permute4x64_pd(vdf, 170);
    vdf_4 = _mm256_permute4x64_pd(vdf, 255);

    vpj_1 = _mm256_load_pd((double*)(&z[j_1][PX]));
    vpj_2 = _mm256_load_pd((double*)(&z[j_2][PX]));
    vpj_3 = _mm256_load_pd((double*)(&z[j_3][PX]));
    vpj_4 = _mm256_load_pd((double*)(&z[j_4][PX]));

    vpi += vdq_1 * vdf_1;
    vpi += vdq_2 * vdf_2;
    vpi += vdq_3 * vdf_3;
    vpi += vdq_4 * vdf_4;
    vpj_1 -= vdq_1 * vdf_1;
    vpj_2 -= vdq_2 * vdf_2;
    vpj_3 -= vdq_3 * vdf_3;
    vpj_4 -= vdq_4 * vdf_4;

    _mm256_store_pd((double*)(&z[j_1][PX]), vpj_1);
    _mm256_store_pd((double*)(&z[j_2][PX]), vpj_2);
    _mm256_store_pd((double*)(&z[j_3][PX]), vpj_3);
    _mm256_store_pd((double*)(&z[j_4][PX]), vpj_4);

    _mm256_store_pd((double*)(&z[i][PX]), vpi);

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
//------------------------------------------------------------------------
#define putsz(x) puts(z[j_1][x],z[j_2][x],z[j_3][x],z[j_4][x],z[j_5][x],z[j_6][x],z[j_7][x],z[j_8][x])
#define putsval(x) puts(x##_1,x##_2,x##_3,x##_4,x##_5,x##_6,x##_7,x##_8);
//------------------------------------------------------------------------
// c.f. https://github.com/kohnakagawa/lj_knl
//------------------------------------------------------------------------
#ifdef AVX512
void
force_avx512(void) {
  const int pn = particle_number;
  const v8df vc24 = _mm512_set1_pd(24.0 * dt);
  const v8df vc48 = _mm512_set1_pd(48.0 * dt);
  const v8df vcl2  = _mm512_set1_pd(CL2);
  const v8df vzero = _mm512_setzero_pd();
  for (int i = 0; i < pn; i++) {
    const double qix = z[i][X];
    const double qiy = z[i][Y];
    const double qiz = z[i][Z];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    const v8df vqxi = _mm512_set1_pd(z[i][X]);
    const v8df vqyi = _mm512_set1_pd(z[i][Y]);
    const v8df vqzi = _mm512_set1_pd(z[i][Z]);
    const int kp = pointer[i];
    for (int k = 0; k < (np / 8) * 8; k += 8) {
      auto vindex = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
      vindex = _mm256_slli_epi32(vindex, 3);
      const v8df vqxj = _mm512_i32gather_pd(vindex, &(z[0][X]), 8);
      const v8df vqyj = _mm512_i32gather_pd(vindex, &(z[0][Y]), 8);
      const v8df vqzj = _mm512_i32gather_pd(vindex, &(z[0][Z]), 8);
      const v8df vdx = vqxj - vqxi;
      const v8df vdy = vqyj - vqyi;
      const v8df vdz = vqzj - vqzi;
      const v8df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      const v8df vr6 =  vr2 * vr2 * vr2;
      v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      const auto vmask = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      vdf = _mm512_mask_blend_pd(vmask, vzero, vdf);

      v8df vpxj = _mm512_i32gather_pd(vindex, &(z[0][PX]), 8);
      v8df vpyj = _mm512_i32gather_pd(vindex, &(z[0][PY]), 8);
      v8df vpzj = _mm512_i32gather_pd(vindex, &(z[0][PZ]), 8);
      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;

      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      _mm512_i32scatter_pd(&(z[0][PX]), vindex, vpxj, 8);
      _mm512_i32scatter_pd(&(z[0][PY]), vindex, vpyj, 8);
      _mm512_i32scatter_pd(&(z[0][PZ]), vindex, vpzj, 8);
    }
    pfx = _mm512_reduce_add_pd(vpxi);
    pfy = _mm512_reduce_add_pd(vpyi);
    pfz = _mm512_reduce_add_pd(vpzi);
    for (int k = (np / 8) * 8; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = z[j][X] - qix;
      double dy = z[j][Y] - qiy;
      double dz = z[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
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
// ループ端数最適化については以下のコードを参照した
// c.f. https://github.com/kohnakagawa/lj_knl
//----------------------------------------------------------------------
void
force_avx512_loopopt(void) {
  const int pn = particle_number;
  const v8df vc24 = _mm512_set1_pd(24.0 * dt);
  const v8df vc48 = _mm512_set1_pd(48.0 * dt);
  const v8df vcl2  = _mm512_set1_pd(CL2);
  const v8df vzero = _mm512_setzero_pd();
  const auto vpitch = _mm512_set1_epi64(8);
  for (int i = 0; i < pn; i++) {
    const int np = number_of_partners[i];
    const auto vnp = _mm512_set1_epi64(np);
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    const v8df vqxi = _mm512_set1_pd(z[i][X]);
    const v8df vqyi = _mm512_set1_pd(z[i][Y]);
    const v8df vqzi = _mm512_set1_pd(z[i][Z]);
    auto vk_idx = _mm512_set_epi64(7LL, 6LL, 5LL, 4LL, 3LL, 2LL, 1LL, 0LL);
    const int kp = pointer[i];
    for (int k = 0; k < np; k += 8) {
      const auto mask_loop = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
      auto vindex = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
      vindex = _mm256_slli_epi32(vindex, 3);
      const v8df vqxj = _mm512_i32gather_pd(vindex, &(z[0][X]), 8);
      const v8df vqyj = _mm512_i32gather_pd(vindex, &(z[0][Y]), 8);
      const v8df vqzj = _mm512_i32gather_pd(vindex, &(z[0][Z]), 8);
      v8df vpxj = _mm512_i32gather_pd(vindex, &(z[0][PX]), 8);
      v8df vpyj = _mm512_i32gather_pd(vindex, &(z[0][PY]), 8);
      v8df vpzj = _mm512_i32gather_pd(vindex, &(z[0][PZ]), 8);
      const v8df vdx = vqxj - vqxi;
      const v8df vdy = vqyj - vqyi;
      const v8df vdz = vqzj - vqzi;
      const v8df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      const v8df vr6 =  vr2 * vr2 * vr2;
      v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      const auto mask_cutoff = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      const auto mask = _mm512_kand(mask_cutoff, mask_loop);
      vdf = _mm512_mask_blend_pd(mask, vzero, vdf);

      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;

      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      _mm512_mask_i32scatter_pd(&(z[0][PX]), mask_loop, vindex, vpxj, 8);
      _mm512_mask_i32scatter_pd(&(z[0][PY]), mask_loop, vindex, vpyj, 8);
      _mm512_mask_i32scatter_pd(&(z[0][PZ]), mask_loop, vindex, vpzj, 8);

      vk_idx = _mm512_add_epi64(vk_idx, vpitch);
    }
    z[i][PX] += _mm512_reduce_add_pd(vpxi);
    z[i][PY] += _mm512_reduce_add_pd(vpyi);
    z[i][PZ] += _mm512_reduce_add_pd(vpzi);
  }
}
//----------------------------------------------------------------------
void
force_avx512_loopopt_swp(void) {
  const int pn = particle_number;
  const v8df vc24 = _mm512_set1_pd(24.0 * dt);
  const v8df vc48 = _mm512_set1_pd(48.0 * dt);
  const v8df vcl2  = _mm512_set1_pd(CL2);
  const v8df vzero = _mm512_setzero_pd();
  const auto vpitch = _mm512_set1_epi64(8);
  for (int i = 0; i < pn; i++) {
    const int np = number_of_partners[i];
    const auto vnp = _mm512_set1_epi64(np);
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    v8df vqxi = _mm512_set1_pd(z[i][X]);
    v8df vqyi = _mm512_set1_pd(z[i][Y]);
    v8df vqzi = _mm512_set1_pd(z[i][Z]);
    auto vk_idx = _mm512_set_epi64(7LL, 6LL, 5LL, 4LL, 3LL, 2LL, 1LL, 0LL);
    const int kp = pointer[i];
    int k = 0;
    auto mask_loop = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
    auto vindex = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
    vindex = _mm256_slli_epi32(vindex, 3);
    v8df vqxj = _mm512_i32gather_pd(vindex, &(z[0][X]), 8);
    v8df vqyj = _mm512_i32gather_pd(vindex, &(z[0][Y]), 8);
    v8df vqzj = _mm512_i32gather_pd(vindex, &(z[0][Z]), 8);
    v8df vdx = vqxj - vqxi;
    v8df vdy = vqyj - vqyi;
    v8df vdz = vqzj - vqzi;
    v8df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
    v8df vr6 =  vr2 * vr2 * vr2;
    v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
    auto mask_cutoff = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
    auto mask = _mm512_kand(mask_cutoff, mask_loop);
    vdf = _mm512_mask_blend_pd(mask, vzero, vdf);
    v8df vpxj, vpyj, vpzj;
    for (k = 8; k < np; k += 8) {
      vk_idx = _mm512_add_epi64(vk_idx, vpitch);
      auto mask_loop_b = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
      auto vindex_b = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
      vindex_b = _mm256_slli_epi32(vindex_b, 3);
      vqxj = _mm512_i32gather_pd(vindex_b, &(z[0][X]), 8);
      vqyj = _mm512_i32gather_pd(vindex_b, &(z[0][Y]), 8);
      vqzj = _mm512_i32gather_pd(vindex_b, &(z[0][Z]), 8);
      auto vdx_b = vqxj - vqxi;
      auto vdy_b = vqyj - vqyi;
      auto vdz_b = vqzj - vqzi;
      vr2 = vdx_b * vdx_b + vdy_b * vdy_b + vdz_b * vdz_b;
      vr6 =  vr2 * vr2 * vr2;
      mask_cutoff = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      mask = _mm512_kand(mask_cutoff, mask_loop_b);

      vpxj = _mm512_i32gather_pd(vindex, &(z[0][PX]), 8);
      vpyj = _mm512_i32gather_pd(vindex, &(z[0][PY]), 8);
      vpzj = _mm512_i32gather_pd(vindex, &(z[0][PZ]), 8);


      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;
      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      _mm512_mask_i32scatter_pd(&(z[0][PX]), mask_loop, vindex, vpxj, 8);
      _mm512_mask_i32scatter_pd(&(z[0][PY]), mask_loop, vindex, vpyj, 8);
      _mm512_mask_i32scatter_pd(&(z[0][PZ]), mask_loop, vindex, vpzj, 8);

      // --8<--
      mask_loop = mask_loop_b;
      vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      vdf = _mm512_mask_blend_pd(mask, vzero, vdf);
      vindex = vindex_b;
      vdx = vdx_b;
      vdy = vdy_b;
      vdz = vdz_b;
    }
    vpxj = _mm512_i32gather_pd(vindex, &(z[0][PX]), 8);
    vpyj = _mm512_i32gather_pd(vindex, &(z[0][PY]), 8);
    vpzj = _mm512_i32gather_pd(vindex, &(z[0][PZ]), 8);
    vpxj -= vdf * vdx;
    vpyj -= vdf * vdy;
    vpzj -= vdf * vdz;
    vpxi += vdf * vdx;
    vpyi += vdf * vdy;
    vpzi += vdf * vdz;

    _mm512_mask_i32scatter_pd(&(z[0][PX]), mask_loop, vindex, vpxj, 8);
    _mm512_mask_i32scatter_pd(&(z[0][PY]), mask_loop, vindex, vpyj, 8);
    _mm512_mask_i32scatter_pd(&(z[0][PZ]), mask_loop, vindex, vpzj, 8);

    z[i][PX] += _mm512_reduce_add_pd(vpxi);
    z[i][PY] += _mm512_reduce_add_pd(vpyi);
    z[i][PZ] += _mm512_reduce_add_pd(vpzi);
  }
}
//----------------------------------------------------------------------
void
force_avx512_gatheronly(void) {
  const int pn = particle_number;
  const v8df vc24 = _mm512_set1_pd(24.0 * dt);
  const v8df vc48 = _mm512_set1_pd(48.0 * dt);
  const v8df vcl2  = _mm512_set1_pd(CL2);
  const v8df vzero = _mm512_setzero_pd();
  const auto vpitch = _mm512_set1_epi64(8);
  for (int i = 0; i < pn; i++) {
    const int np = number_of_partners[i];
    const auto vnp = _mm512_set1_epi64(np);
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    const v8df vqxi = _mm512_set1_pd(z[i][X]);
    const v8df vqyi = _mm512_set1_pd(z[i][Y]);
    const v8df vqzi = _mm512_set1_pd(z[i][Z]);
    auto vk_idx = _mm512_set_epi64(7LL, 6LL, 5LL, 4LL, 3LL, 2LL, 1LL, 0LL);
    const int kp = pointer[i];
    for (int k = 0; k < np; k += 8) {

      const auto mask_loop = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
      auto vindex2 = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
      auto vindex = _mm256_slli_epi32(vindex2, 3);
      const int j_1 = sorted_list[kp + k];
      const int j_2 = sorted_list[kp + k + 1];
      const int j_3 = sorted_list[kp + k + 2];
      const int j_4 = sorted_list[kp + k + 3];
      const int j_5 = sorted_list[kp + k + 4];
      const int j_6 = sorted_list[kp + k + 5];
      const int j_7 = sorted_list[kp + k + 6];
      const int j_8 = sorted_list[kp + k + 7];

      const v8df vqxj = _mm512_i32gather_pd(vindex, &(z[0][X]), 8);
      const v8df vqyj = _mm512_i32gather_pd(vindex, &(z[0][Y]), 8);
      const v8df vqzj = _mm512_i32gather_pd(vindex, &(z[0][Z]), 8);
      v8df vpxj = _mm512_i32gather_pd(vindex, &(z[0][PX]), 8);
      v8df vpyj = _mm512_i32gather_pd(vindex, &(z[0][PY]), 8);
      v8df vpzj = _mm512_i32gather_pd(vindex, &(z[0][PZ]), 8);
      const v8df vdx = vqxj - vqxi;
      const v8df vdy = vqyj - vqyi;
      const v8df vdz = vqzj - vqzi;
      const v8df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      const v8df vr6 =  vr2 * vr2 * vr2;
      v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      const auto mask_cutoff = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      const auto mask = _mm512_kand(mask_cutoff, mask_loop);
      vdf = _mm512_mask_blend_pd(mask, vzero, vdf);

      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;

      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      v8df t1 = vpxj;
      v8df t2 = vpyj;
      v8df t3 = vpzj;
      v8df t4 = _mm512_setzero_pd();
      transpose_4x4(t1, t2, t3, t4);
      v4df pj_1 = _mm512_extractf64x4_pd(t1, 0);
      v4df pj_2 = _mm512_extractf64x4_pd(t2, 0);
      v4df pj_3 = _mm512_extractf64x4_pd(t3, 0);
      v4df pj_4 = _mm512_extractf64x4_pd(t4, 0);
      v4df pj_5 = _mm512_extractf64x4_pd(t1, 1);
      v4df pj_6 = _mm512_extractf64x4_pd(t2, 1);
      v4df pj_7 = _mm512_extractf64x4_pd(t3, 1);
      v4df pj_8 = _mm512_extractf64x4_pd(t4, 1);
      /*
            __mmask8 k1 = (k<np)? 255: 0;
            __mmask8 k2 = (k+1<np)? 255: 0;
            __mmask8 k3 = (k+2<np)? 255: 0;
            __mmask8 k4 = (k+3<np)? 255: 0;
            __mmask8 k5 = (k+4<np)? 255: 0;
            __mmask8 k6 = (k+5<np)? 255: 0;
            __mmask8 k7 = (k+6<np)? 255: 0;
            __mmask8 k8 = (k+7<np)? 255: 0;
            _mm256_mask_store_pd(&(z[j_1][PX]), k1, pj_1);
            _mm256_mask_store_pd(&(z[j_2][PX]), k2, pj_2);
            _mm256_mask_store_pd(&(z[j_3][PX]), k3, pj_3);
            _mm256_mask_store_pd(&(z[j_4][PX]), k4, pj_4);
            _mm256_mask_store_pd(&(z[j_5][PX]), k5, pj_5);
            _mm256_mask_store_pd(&(z[j_6][PX]), k6, pj_6);
            _mm256_mask_store_pd(&(z[j_7][PX]), k7, pj_7);
            _mm256_mask_store_pd(&(z[j_8][PX]), k8, pj_8);
      */
      if (k  < np)_mm256_store_pd(&(z[j_1][PX]), pj_1);
      if (k + 1 < np)_mm256_store_pd(&(z[j_2][PX]), pj_2);
      if (k + 2 < np)_mm256_store_pd(&(z[j_3][PX]), pj_3);
      if (k + 3 < np)_mm256_store_pd(&(z[j_4][PX]), pj_4);
      if (k + 4 < np)_mm256_store_pd(&(z[j_5][PX]), pj_5);
      if (k + 5 < np)_mm256_store_pd(&(z[j_6][PX]), pj_6);
      if (k + 6 < np)_mm256_store_pd(&(z[j_7][PX]), pj_7);
      if (k + 7 < np)_mm256_store_pd(&(z[j_8][PX]), pj_8);

      vk_idx = _mm512_add_epi64(vk_idx, vpitch);
    }
    z[i][PX] += _mm512_reduce_add_pd(vpxi);
    z[i][PY] += _mm512_reduce_add_pd(vpyi);
    z[i][PZ] += _mm512_reduce_add_pd(vpzi);
  }
}

//----------------------------------------------------------------------
void
force_avx512_transpose(void) {
  const int pn = particle_number;
  const v8df vc24 = _mm512_set1_pd(24.0 * dt);
  const v8df vc48 = _mm512_set1_pd(48.0 * dt);
  const v8df vcl2  = _mm512_set1_pd(CL2);
  const v8df vzero = _mm512_setzero_pd();
  const auto vpitch = _mm512_set1_epi64(8);
  for (int i = 0; i < pn; i++) {
    const int np = number_of_partners[i];
    const auto vnp = _mm512_set1_epi64(np);
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    const v8df vqxi = _mm512_set1_pd(z[i][X]);
    const v8df vqyi = _mm512_set1_pd(z[i][Y]);
    const v8df vqzi = _mm512_set1_pd(z[i][Z]);
    auto vk_idx = _mm512_set_epi64(7LL, 6LL, 5LL, 4LL, 3LL, 2LL, 1LL, 0LL);
    const int kp = pointer[i];
    for (int k = 0; k < np; k += 8) {
      const auto mask_loop = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
      const int j_1 = sorted_list[kp + k];
      const int j_2 = sorted_list[kp + k + 1];
      const int j_3 = sorted_list[kp + k + 2];
      const int j_4 = sorted_list[kp + k + 3];
      const int j_5 = sorted_list[kp + k + 4];
      const int j_6 = sorted_list[kp + k + 5];
      const int j_7 = sorted_list[kp + k + 6];
      const int j_8 = sorted_list[kp + k + 7];
      auto vindex = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
      vindex = _mm256_slli_epi32(vindex, 3);
      v8df vz_1 = _mm512_load_pd(&z[j_1][0]);
      v8df vz_2 = _mm512_load_pd(&z[j_2][0]);
      v8df vz_3 = _mm512_load_pd(&z[j_3][0]);
      v8df vz_4 = _mm512_load_pd(&z[j_4][0]);
      transpose_4x4(vz_1, vz_2, vz_3, vz_4);
      v8df vz_5 = _mm512_load_pd(&z[j_5][0]);
      v8df vz_6 = _mm512_load_pd(&z[j_6][0]);
      v8df vz_7 = _mm512_load_pd(&z[j_7][0]);
      v8df vz_8 = _mm512_load_pd(&z[j_8][0]);
      transpose_4x4(vz_5, vz_6, vz_7, vz_8);
      v8df vqxj = _mm512_insertf64x4(vz_1, _mm512_extractf64x4_pd(vz_5, 0), 1);
      v8df vqyj = _mm512_insertf64x4(vz_2, _mm512_extractf64x4_pd(vz_6, 0), 1);
      v8df vqzj = _mm512_insertf64x4(vz_3, _mm512_extractf64x4_pd(vz_7, 0), 1);
      v8df vpxj = _mm512_insertf64x4(vz_5, _mm512_extractf64x4_pd(vz_1, 1), 0);
      v8df vpyj = _mm512_insertf64x4(vz_6, _mm512_extractf64x4_pd(vz_2, 1), 0);
      v8df vpzj = _mm512_insertf64x4(vz_7, _mm512_extractf64x4_pd(vz_3, 1), 0);
      const v8df vdx = vqxj - vqxi;
      const v8df vdy = vqyj - vqyi;
      const v8df vdz = vqzj - vqzi;
      const v8df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      const v8df vr6 =  vr2 * vr2 * vr2;
      v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      const auto mask_cutoff = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      const auto mask = _mm512_kand(mask_cutoff, mask_loop);
      vdf = _mm512_mask_blend_pd(mask, vzero, vdf);

      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;

      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      _mm512_mask_i32scatter_pd(&(z[0][PX]), mask_loop, vindex, vpxj, 8);
      _mm512_mask_i32scatter_pd(&(z[0][PY]), mask_loop, vindex, vpyj, 8);
      _mm512_mask_i32scatter_pd(&(z[0][PZ]), mask_loop, vindex, vpzj, 8);

      vk_idx = _mm512_add_epi64(vk_idx, vpitch);
    }
    z[i][PX] += _mm512_reduce_add_pd(vpxi);
    z[i][PY] += _mm512_reduce_add_pd(vpyi);
    z[i][PZ] += _mm512_reduce_add_pd(vpzi);
  }
}
//----------------------------------------------------------------------
#endif //AVX512
int
main(void) {
  AoSDataManager aosdm(p, q);
  init(&aosdm, particle_number);
  check_pairlist(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, &aosdm);
  sortpair(particle_number, number_of_pairs, number_of_partners, i_particles, j_particles, pointer, sorted_list);
#ifdef PAIR
  measure(&force_pair, "pair", particle_number);
  aosdm.print_results(particle_number);
#elif SORTED
  measure(&force_sorted, "sorted", particle_number);
  aosdm.print_results(particle_number);
#elif AVX2
  measure(&force_avx2, "avx2", particle_number);
  aosdm.print_results(particle_number);
#elif AVX2_SWP
  measure(&force_avx2_swp, "avx2_swp", particle_number);
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
#elif SORTED_Z_AVX2_SWP
  copy_to_z();
  measure(&force_sorted_z_avx2_swp, "sorted_z_avx2_swp", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#elif AVX512_SIMPLE
  copy_to_z();
  measure(&force_avx512, "avx512", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#elif AVX512_LOOPOPT
  copy_to_z();
  measure(&force_avx512_loopopt, "avx512_loopopt", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#elif AVX512_LOOPOPT_SWP
  copy_to_z();
  measure(&force_avx512_loopopt_swp, "avx512_loopopt_swp", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#elif AVX512_GATHERONLY
  copy_to_z();
  measure(&force_avx512_gatheronly, "avx512_gatheronly", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#elif AVX512_TRANSPOSE
  copy_to_z();
  measure(&force_avx512_transpose, "avx512_transpose", particle_number);
  copy_from_z();
  aosdm.print_results(particle_number);
#else
  measure(&force_pair, "pair", particle_number);
  measure(&force_sorted, "sorted", particle_number);
  measure(&force_swp, "sorted_swp", particle_number);
  measure(&force_avx2, "avx2", particle_number);
  measure(&force_avx2_swp, "avx2_swp", particle_number);
  copy_to_z();
  measure(&force_sorted_z, "sorted_z", particle_number);
  measure(&force_sorted_z_avx2, "sorted_z_avx2", particle_number);
#ifdef AVX512
  measure(&force_avx512, "avx512", particle_number);
  measure(&force_avx512_loopopt, "avx512_loopopt", particle_number);
  measure(&force_avx512_loopopt_swp, "avx512_loopopt_swp", particle_number);
  measure(&force_avx512_transpose, "avx512_transpose", particle_number);
#endif //AVX512
#endif
}
//----------------------------------------------------------------------
