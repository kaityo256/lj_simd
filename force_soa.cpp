#include <cstdio>
#include <x86intrin.h>
#include <random>
#include "conf.hpp"
//----------------------------------------------------------------------
double q[D][N];
double p[D][N] = {};
int particle_number = 0;
int number_of_pairs = 0;
int number_of_partners[N] = {};
int i_particles[MAX_PAIRS] = {};
int j_particles[MAX_PAIRS] = {};
int pointer[N], pointer2[N];
int sorted_list[MAX_PAIRS] = {};
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
    //if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    if (r2 > CL2) df = 0.0;
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
    const double qix = q[X][i];
    const double qiy = q[Y][i];
    const double qiz = q[Z][i];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k = 0; k < np; k++) {
      const int j = sorted_list2[kp + k];
      double dx = q[X][j] - qix;
      double dy = q[Y][j] - qiy;
      double dz = q[Z][j] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      //if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      if (r2 > CL2) df = 0.0;
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
    const double qix = q[X][i];
    const double qiy = q[Y][i];
    const double qiz = q[Z][i];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    const v4df vqix = _mm256_set_pd(q[X][i], q[X][i], q[X][i], q[X][i]);
    const v4df vqiy = _mm256_set_pd(q[Y][i], q[Y][i], q[Y][i], q[Y][i]);
    const v4df vqiz = _mm256_set_pd(q[Z][i], q[Z][i], q[Z][i], q[Z][i]);

    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_1 = sorted_list[kp + k];
      const int j_2 = sorted_list[kp + k + 1];
      const int j_3 = sorted_list[kp + k + 2];
      const int j_4 = sorted_list[kp + k + 3];

      const v4df vqjx = _mm256_set_pd(q[X][j_1], q[X][j_2], q[X][j_3], q[X][j_4]);
      const v4df vqjy = _mm256_set_pd(q[Y][j_1], q[Y][j_2], q[Y][j_3], q[Y][j_4]);
      const v4df vqjz = _mm256_set_pd(q[Z][j_1], q[Z][j_2], q[Z][j_3], q[Z][j_4]);

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

      p[X][j_1] -= dfx[3];
      p[Y][j_1] -= dfy[3];
      p[Z][j_1] -= dfz[3];

      pfx += dfx[2];
      pfy += dfy[2];
      pfz += dfz[2];
      p[X][j_2] -= dfx[2];
      p[Y][j_2] -= dfy[2];
      p[Z][j_2] -= dfz[2];

      pfx += dfx[1];
      pfy += dfy[1];
      pfz += dfz[1];
      p[X][j_3] -= dfx[1];
      p[Y][j_3] -= dfy[1];
      p[Z][j_3] -= dfz[1];

      pfx += dfx[0];
      pfy += dfy[0];
      pfz += dfz[0];
      p[X][j_4] -= dfx[0];
      p[Y][j_4] -= dfy[0];
      p[Z][j_4] -= dfz[0];
    }
    p[X][i] += pfx;
    p[Y][i] += pfy;
    p[Z][i] += pfz;
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[X][j] - qix;
      double dy = q[Y][j] - qiy;
      double dz = q[Z][j] - qiz;
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
    const double qix = q[X][i];
    const double qiy = q[Y][i];
    const double qiz = q[Z][i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    int ja = sorted_list[kp];
    double dxa = q[X][ja] - qix;
    double dya = q[Y][ja] - qiy;
    double dza = q[Z][ja] - qiz;
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
      dxa = q[X][ja] - qix;
      dya = q[Y][ja] - qiy;
      dza = q[Z][ja] - qiz;
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
#define putsarray(x) puts(x[j##_1],x[j##_2],x[j##_3],x[j##_4],x[j##_5],x[j##_6],x[j##_7],x[j##_8]);
#define putsval(x) puts(x##_1,x##_2,x##_3,x##_4,x##_5,x##_6,x##_7,x##_8);
//----------------------------------------------------------------------
#ifdef AVX512
//----------------------------------------------------------------------
void
force_avx512(void) {
  const int pn = particle_number;
  const int * __restrict sorted_list2 = sorted_list;
  const v8df vc24 = _mm512_set1_pd(24.0 * dt);
  const v8df vc48 = _mm512_set1_pd(48.0 * dt);
  const v8df vcl2  = _mm512_set1_pd(CL2);
  const v8df vzero = _mm512_setzero_pd();
  for (int i = 0; i < pn; i++) {
    const double qix = q[X][i];
    const double qiy = q[Y][i];
    const double qiz = q[Z][i];
    const v8df vqxi = _mm512_set1_pd(q[X][i]);
    const v8df vqyi = _mm512_set1_pd(q[Y][i]);
    const v8df vqzi = _mm512_set1_pd(q[Z][i]);
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    const int kp = pointer[i];
    for (int k = 0; k < (np / 8) * 8; k += 8) {
      const auto vindex = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
      const v8df vqxj = _mm512_i32gather_pd(vindex, &(q[X][0]), 8);
      const v8df vqyj = _mm512_i32gather_pd(vindex, &(q[Y][0]), 8);
      const v8df vqzj = _mm512_i32gather_pd(vindex, &(q[Z][0]), 8);

      v8df vpxj = _mm512_i32gather_pd(vindex, &(p[X][0]), 8);
      v8df vpyj = _mm512_i32gather_pd(vindex, &(p[Y][0]), 8);
      v8df vpzj = _mm512_i32gather_pd(vindex, &(p[Z][0]), 8);

      const v8df vdx = vqxj - vqxi;
      const v8df vdy = vqyj - vqyi;
      const v8df vdz = vqzj - vqzi;
      const v8df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      const v8df vr6 =  vr2 * vr2 * vr2;
      v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      const auto vmask = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      vdf = _mm512_mask_blend_pd(vmask, vzero, vdf);

      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;

      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      _mm512_i32scatter_pd(&(p[X][0]), vindex, vpxj, 8);
      _mm512_i32scatter_pd(&(p[Y][0]), vindex, vpyj, 8);
      _mm512_i32scatter_pd(&(p[Z][0]), vindex, vpzj, 8);
    }
    pfx = _mm512_reduce_add_pd(vpxi);
    pfy = _mm512_reduce_add_pd(vpyi);
    pfz = _mm512_reduce_add_pd(vpzi);

    for (int k = (np / 8) * 8; k < np; k++) {
      const int j = sorted_list2[kp + k];
      double dx = q[X][j] - qix;
      double dy = q[Y][j] - qiy;
      double dz = q[Z][j] - qiz;
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
// AVX-512 with loop optimization
// c.f. https://github.com/kohnakagawa/lj_knl/blob/master/pd/force_soa.cpp
void
force_avx512_loopopt(void) {
  const int pn = particle_number;
  const v8df vc24 = _mm512_set1_pd(24.0 * dt);
  const v8df vc48 = _mm512_set1_pd(48.0 * dt);
  const v8df vcl2  = _mm512_set1_pd(CL2);
  const v8df vzero = _mm512_setzero_pd();
  const auto vpitch = _mm512_set1_epi64(8);
  for (int i = 0; i < pn; i++) {
    const v8df vqxi = _mm512_set1_pd(q[X][i]);
    const v8df vqyi = _mm512_set1_pd(q[Y][i]);
    const v8df vqzi = _mm512_set1_pd(q[Z][i]);
    const int np = number_of_partners[i];
    const auto vnp = _mm512_set1_epi64(np);
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    auto vk_idx = _mm512_set_epi64(7LL, 6LL, 5LL, 4LL, 3LL, 2LL, 1LL, 0LL);
    const int kp = pointer[i];
    for (int k = 0; k < np; k += 8) {
      const auto mask_loop = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
      const auto vindex = _mm256_lddqu_si256((const __m256i*)(&sorted_list[kp + k]));
      const v8df vqxj = _mm512_i32gather_pd(vindex, &(q[X][0]), 8);
      const v8df vqyj = _mm512_i32gather_pd(vindex, &(q[Y][0]), 8);
      const v8df vqzj = _mm512_i32gather_pd(vindex, &(q[Z][0]), 8);

      v8df vpxj = _mm512_i32gather_pd(vindex, &(p[X][0]), 8);
      v8df vpyj = _mm512_i32gather_pd(vindex, &(p[Y][0]), 8);
      v8df vpzj = _mm512_i32gather_pd(vindex, &(p[Z][0]), 8);

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

      _mm512_mask_i32scatter_pd(&(p[X][0]), mask_loop, vindex, vpxj, 8);
      _mm512_mask_i32scatter_pd(&(p[Y][0]), mask_loop, vindex, vpyj, 8);
      _mm512_mask_i32scatter_pd(&(p[Z][0]), mask_loop, vindex, vpzj, 8);
      vk_idx = _mm512_add_epi64(vk_idx, vpitch);
    }
    p[X][i] += _mm512_reduce_add_pd(vpxi);
    p[Y][i] += _mm512_reduce_add_pd(vpyi);
    p[Z][i] += _mm512_reduce_add_pd(vpzi);
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
    const v8df vqxi = _mm512_set1_pd(q[X][i]);
    const v8df vqyi = _mm512_set1_pd(q[Y][i]);
    const v8df vqzi = _mm512_set1_pd(q[Z][i]);
    const int np = number_of_partners[i];
    const auto vnp = _mm512_set1_epi64(np);
    v8df vpxi = _mm512_setzero_pd();
    v8df vpyi = _mm512_setzero_pd();
    v8df vpzi = _mm512_setzero_pd();
    auto vk_idx = _mm512_set_epi64(7LL, 6LL, 5LL, 4LL, 3LL, 2LL, 1LL, 0LL);
    const int kp = pointer[i];
    int k = 0;
    auto mask_loop = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
    auto vindex = _mm256_lddqu_si256(( __m256i*)(&sorted_list[kp + k]));
    v8df vqxj = _mm512_i32gather_pd(vindex, &(q[X][0]), 8);
    v8df vqyj = _mm512_i32gather_pd(vindex, &(q[Y][0]), 8);
    v8df vqzj = _mm512_i32gather_pd(vindex, &(q[Z][0]), 8);

    v8df vpxj = _mm512_i32gather_pd(vindex, &(p[X][0]), 8);
    v8df vpyj = _mm512_i32gather_pd(vindex, &(p[Y][0]), 8);
    v8df vpzj = _mm512_i32gather_pd(vindex, &(p[Z][0]), 8);

    v8df vdx = vqxj - vqxi;
    v8df vdy = vqyj - vqyi;
    v8df vdz = vqzj - vqzi;
    v8df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
    v8df vr6 =  vr2 * vr2 * vr2;
    v8df vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
    auto mask_cutoff = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
    auto mask = _mm512_kand(mask_cutoff, mask_loop);
    vdf = _mm512_mask_blend_pd(mask, vzero, vdf);
    for (k = 8; k < np; k += 8) {
      vk_idx = _mm512_add_epi64(vk_idx, vpitch);
      auto mask_loop_b = _mm512_cmp_epi64_mask(vk_idx, vnp, _MM_CMPINT_LT);
      auto vindex_b = _mm256_lddqu_si256(( __m256i*)(&sorted_list[kp + k]));
      vqxj = _mm512_i32gather_pd(vindex_b, &(q[X][0]), 8);
      vqyj = _mm512_i32gather_pd(vindex_b, &(q[Y][0]), 8);
      vqzj = _mm512_i32gather_pd(vindex_b, &(q[Z][0]), 8);
      auto vdx_b = vqxj - vqxi;
      auto vdy_b = vqyj - vqyi;
      auto vdz_b = vqzj - vqzi;
      vr2 = vdx_b * vdx_b + vdy_b * vdy_b + vdz_b * vdz_b;
      vr6 =  vr2 * vr2 * vr2;
      mask_cutoff = _mm512_cmp_pd_mask(vr2, vcl2, _CMP_LE_OS);
      mask = _mm512_kand(mask_cutoff, mask_loop_b);
      // ---- 8<  ------
      vpxj -= vdf * vdx;
      vpyj -= vdf * vdy;
      vpzj -= vdf * vdz;

      vpxi += vdf * vdx;
      vpyi += vdf * vdy;
      vpzi += vdf * vdz;

      _mm512_mask_i32scatter_pd(&(p[X][0]), mask_loop, vindex, vpxj, 8);
      _mm512_mask_i32scatter_pd(&(p[Y][0]), mask_loop, vindex, vpyj, 8);
      _mm512_mask_i32scatter_pd(&(p[Z][0]), mask_loop, vindex, vpzj, 8);
      // ---- 8<  ------

      vpxj = _mm512_i32gather_pd(vindex_b, &(p[X][0]), 8);
      vpyj = _mm512_i32gather_pd(vindex_b, &(p[Y][0]), 8);
      vpzj = _mm512_i32gather_pd(vindex_b, &(p[Z][0]), 8);

      vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2);
      vdf = _mm512_mask_blend_pd(mask, vzero, vdf);
      // To next loop
      mask_loop = mask_loop_b;
      vindex = vindex_b;
      vdx = vdx_b;
      vdy = vdy_b;
      vdz = vdz_b;
    }
    vpxj -= vdf * vdx;
    vpyj -= vdf * vdy;
    vpzj -= vdf * vdz;

    vpxi += vdf * vdx;
    vpyi += vdf * vdy;
    vpzi += vdf * vdz;

    _mm512_mask_i32scatter_pd(&(p[X][0]), mask_loop, vindex, vpxj, 8);
    _mm512_mask_i32scatter_pd(&(p[Y][0]), mask_loop, vindex, vpyj, 8);
    _mm512_mask_i32scatter_pd(&(p[Z][0]), mask_loop, vindex, vpzj, 8);

    p[X][i] += _mm512_reduce_add_pd(vpxi);
    p[Y][i] += _mm512_reduce_add_pd(vpyi);
    p[Z][i] += _mm512_reduce_add_pd(vpzi);
  }
}
//----------------------------------------------------------------------
#endif //AVX512
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
#elif AVX512_SIMPLE
  measure(&force_avx512, "avx512", particle_number);
  soadm.print_results(particle_number);
#elif AVX512_LOOPOPT
  measure(&force_avx512_loopopt, "avx512_loopopt", particle_number);
  soadm.print_results(particle_number);
#elif AVX512_LOOPOPT_SWP
  measure(&force_avx512_loopopt_swp, "avx512_loopopt_swp", particle_number);
  soadm.print_results(particle_number);
#else
  measure(&force_pair, "pair", particle_number);
  measure(&force_sorted, "sorted", particle_number);
  measure(&force_swp, "sorted_swp", particle_number);
  measure(&force_avx2, "avx2", particle_number);
#ifdef AVX512
  measure(&force_avx512, "avx512", particle_number);
  measure(&force_avx512_loopopt, "avx512_loopopt", particle_number);
  measure(&force_avx512_loopopt_swp, "avx512_loopopt_swp", particle_number);
#endif //AVX512
#endif
}
//----------------------------------------------------------------------
