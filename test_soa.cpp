#include <stdio.h>
#include <immintrin.h>
#include <random>
#include <math.h>
#include <sys/time.h>
//----------------------------------------------------------------------
const double density = 1.0;
const int N = 400000;
const int MAX_PAIRS = 30 * N;
double L = 50.0;
const double dt = 0.001;
const int D = 4;
enum {X, Y, Z};
//double q[N][D];
//double p[N][D];
double qx[N],qy[N],qz[N];
double px[N],py[N],pz[N];

int particle_number = 0;
int number_of_pairs = 0;
int number_of_partners[N];
int i_particles[MAX_PAIRS];
int j_particles[MAX_PAIRS];
int pointer[N], pointer2[N];
int sorted_list[MAX_PAIRS];

const double CUTOFF_LENGTH = 3.0;
const double SEARCH_LENGTH = 3.3;
const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
//----------------------------------------------------------------------
typedef double v4df __attribute__((vector_size(32)));
//----------------------------------------------------------------------
void
print256(v4df r) {
  double *a = (double*)(&r);
  printf("%.10f %.10f %.10f %.10f\n", a[0], a[1], a[2], a[3]);
}
//----------------------------------------------------------------------
void
add_particle(double x, double y, double z) {
  static std::mt19937 mt(2);
  std::uniform_real_distribution<double> ud(0.0, 0.1);
  qx[particle_number] = x + ud(mt);
  qy[particle_number] = y + ud(mt);
  qz[particle_number] = z + ud(mt);
  particle_number++;
}
//----------------------------------------------------------------------
double
myclock(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}
//----------------------------------------------------------------------
void
register_pair(int index1, int index2) {
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
void
makepair(void) {
  const double SL2 = SEARCH_LENGTH * SEARCH_LENGTH;
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    number_of_partners[i] = 0;
  }
  for (int i = 0; i < particle_number - 1; i++) {
    for (int j = i + 1; j < particle_number; j++) {
      const double dx = qx[i] - qx[j];
      const double dy = qy[i] - qy[j];
      const double dz = qz[i] - qz[j];
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < SL2) {
        register_pair(i, j);
      }
    }
  }
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
init(void) {
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L / s);
  int sy = static_cast<int>(L / s);
  int sz = static_cast<int>(L / s);
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        double x = ix*s;
        double y = iy*s;
        double z = iz*s;
        add_particle(x     ,y   ,z);
        add_particle(x     ,y+hs,z+hs);
        add_particle(x+hs  ,y   ,z+hs);
        add_particle(x+hs  ,y+hs,z);
      }
    }
  }
  for (int i = 0; i < particle_number; i++) {
    px[i] = 0.0;
    py[i] = 0.0;
    pz[i] = 0.0;
  }
}
//----------------------------------------------------------------------
void
force_pair(void){
  for(int k=0;k<number_of_pairs;k++){
    const int i = i_particles[k];
    const int j = j_particles[k];
    double dx = qx[j] - qx[i];
    double dy = qy[j] - qy[i];
    double dz = qz[j] - qz[i];
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
    px[i] += df * dx;
    py[i] += df * dy;
    pz[i] += df * dz;
    px[j] -= df * dx;
    py[j] -= df * dy;
    pz[j] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
force_sorted(void){
  const int pn =particle_number;
  for (int i=0; i<pn; i++) {
    const double qx_key = qx[i];
    const double qy_key = qy[i];
    const double qz_key = qz[i];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    for (int k=0; k<np; k++) {
      const int j = sorted_list[kp + k];
      double dx = qx[j] - qx_key;
      double dy = qy[j] - qy_key;
      double dz = qz[j] - qz_key;
      double r2 = (dx*dx + dy*dy + dz*dz);
      if (r2 > CL2) continue;
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2))*dt;
      pfx += df*dx;
      pfy += df*dy;
      pfz += df*dz;
      px[j] -= df*dx;
      py[j] -= df*dy;
      pz[j] -= df*dz;
    }
    px[i] += pfx;
    py[i] += pfy;
    pz[i] += pfz;
  }
}
//----------------------------------------------------------------------
void
force_intrin(void){
  const int pn =particle_number;
  const v4df vzero = _mm256_set_pd(0, 0, 0, 0);
  const v4df vcl2 = _mm256_set_pd(CL2, CL2, CL2, CL2);
  const v4df vc24 = _mm256_set_pd(24 * dt, 24 * dt, 24 * dt, 24 * dt);
  const v4df vc48 = _mm256_set_pd(48 * dt, 48 * dt, 48 * dt, 48 * dt);
  for (int i=0; i<pn; i++) {
    const double qx_key = qx[i];
    const double qy_key = qy[i];
    const double qz_key = qz[i];
    const int np = number_of_partners[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    const v4df vqix = _mm256_set_pd(qx[i],qx[i],qx[i],qx[i]);
    const v4df vqiy = _mm256_set_pd(qy[i],qy[i],qy[i],qy[i]);
    const v4df vqiz = _mm256_set_pd(qz[i],qz[i],qz[i],qz[i]);

    const v4df vpix = _mm256_set_pd(px[i],px[i],px[i],px[i]);
    const v4df vpiy = _mm256_set_pd(py[i],py[i],py[i],py[i]);
    const v4df vpiz = _mm256_set_pd(pz[i],pz[i],pz[i],pz[i]);

    for (int k=0; k<(np/4)*4; k+=4) {
      const int j_a = sorted_list[kp + k];
      const int j_b = sorted_list[kp + k + 1];
      const int j_c = sorted_list[kp + k + 2];
      const int j_d = sorted_list[kp + k + 3];

      const v4df vqjx = _mm256_set_pd(qx[j_a],qx[j_b],qx[j_c],qx[j_d]);
      const v4df vqjy = _mm256_set_pd(qy[j_a],qy[j_b],qy[j_c],qy[j_d]);
      const v4df vqjz = _mm256_set_pd(qz[j_a],qz[j_b],qz[j_c],qz[j_d]);

      const v4df vdx = vqjx - vqix;
      const v4df vdy = vqjy - vqiy;
      const v4df vdz = vqjz - vqiz;
      const v4df vr2 = vdx*vdx+vdy*vdy+vdz*vdz;
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

      px[j_a] -= dfx[3];
      py[j_a] -= dfy[3];
      pz[j_a] -= dfz[3];

      pfx += dfx[2];
      pfy += dfy[2];
      pfz += dfz[2];
      px[j_b] -= dfx[2];
      py[j_b] -= dfy[2];
      pz[j_b] -= dfz[2];

      pfx += dfx[1];
      pfy += dfy[1];
      pfz += dfz[1];
      px[j_c] -= dfx[1];
      py[j_c] -= dfy[1];
      pz[j_c] -= dfz[1];

      pfx += dfx[0];
      pfy += dfy[0];
      pfz += dfz[0];
      px[j_d] -= dfx[0];
      py[j_d] -= dfy[0];
      pz[j_d] -= dfz[0];
    }
    px[i] += pfx;
    py[i] += pfy;
    pz[i] += pfz;
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = qx[j] - qx_key;
      double dy = qy[j] - qy_key;
      double dz = qz[j] - qz_key;
      double r2 = (dx*dx + dy*dy + dz*dz);
      if (r2 > CL2) continue;
      double r6 = r2*r2*r2;
      double df = ((24.0*r6-48.0)/(r6*r6*r2))*dt;
      px[i] += df*dx;
      py[i] += df*dy;
      pz[i] += df*dz;
      px[j] -= df*dx;
      py[j] -= df*dy;
      pz[j] -= df*dz;
    }
  }
}
//----------------------------------------------------------------------
void
force_next(void) {
  const int pn = particle_number;
  for (int i = 0; i < pn; i++) {
    const double qx_key = qx[i];
    const double qy_key = qy[i];
    const double qz_key = qz[i];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = pointer[i];
    int ja = sorted_list[kp];
    double dxa = qx[ja] - qx_key;
    double dya = qy[ja] - qy_key;
    double dza = qz[ja] - qz_key;
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
      dxa = qx[ja] - qx_key;
      dya = qy[ja] - qy_key;
      dza = qz[ja] - qz_key;
      if (r2 > CL2)continue;
      pfx += df * dxb;
      pfy += df * dyb;
      pfz += df * dzb;
      px[jb] -= df * dxb;
      py[jb] -= df * dyb;
      pz[jb] -= df * dzb;
      const double r6 = r2 * r2 * r2;
      df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt;
      jb = j;
      dxb = dx;
      dyb = dy;
      dzb = dz;
    }
    px[jb] -= df * dxb;
    py[jb] -= df * dyb;
    pz[jb] -= df * dzb;
    px[i] += pfx + df * dxb;
    py[i] += pfy + df * dyb;
    pz[i] += pfz + df * dzb;
  }
}
//----------------------------------------------------------------------
void
measure(void(*pfunc)(), const char *name) {
  double st = myclock();
  const int LOOP = 100;
  for (int i = 0; i < LOOP; i++) {
    pfunc();
  }
  double t = myclock() - st;
  fprintf(stderr, "N=%d, %s %f [sec]\n", particle_number, name, t);
}
//----------------------------------------------------------------------
int
main(void) {
  init();
  makepair();
#ifdef PAIR
  measure(&force_pair, "pair");
  for (int i = 0; i < 10; i++) {
    printf("%.10f %.10f %.10f\n", px[i], py[i], pz[i]);
  }
#elif INTRIN
  measure(&force_intrin, "intrin");
  for (int i = 0; i < 10; i++) {
    printf("%.10f %.10f %.10f\n", px[i], py[i], pz[i]);
  }
#else
  measure(&force_pair, "pair");
  measure(&force_sorted, "sorted");
  measure(&force_next, "sorted_next");
  measure(&force_intrin, "intrin");
#endif
}
//----------------------------------------------------------------------
