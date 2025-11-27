#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pn_eom.h"
#include "pn_eom_hamiltonians.h"
#include "utils.h"

#define PI 3.14159265359
#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif

// Integral parameters
#define NDIM 3
#define NCOMP 1
#define NVEC 1
#define EPSREL 1e-16
#define EPSABS 1e-20
#define VERBOSE 0
#define SEED 42
#define MINEVAL 50000
#define MAXEVAL 50000
#define KEY 11

typedef struct {
    double **pos;
} IntegralParams;


static int ln_integral(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
    // Userdata
    IntegralParams *params = (IntegralParams*)userdata;

    // Positions passed via userdata
    double pos_a[] = {params->pos[0][0], params->pos[0][1], params->pos[0][2]};
    double pos_b[] = {params->pos[1][0], params->pos[1][1], params->pos[1][2]};
    double pos_c[] = {params->pos[2][0], params->pos[2][1], params->pos[2][2]};
    double pos_d[] = {params->pos[3][0], params->pos[3][1], params->pos[3][2]};

    // New variables for the transformation to an integral over the unit hypercube
    double x_trans[3];
    cubareal rho = xx[0];
    cubareal u = xx[1];
    cubareal v = xx[2];
    if (rho > 1 - 1e-10) rho = 1 - 1e-10;
    x_trans[0] = rho/(1-rho*rho) * sin(PI*u) * cos(2*PI*v);
    x_trans[1] = rho/(1-rho*rho) * sin(PI*u) * sin(2*PI*v);
    x_trans[2] = rho/(1-rho*rho) * cos(PI*u);

    // Jacobian determinant for the transformation to an integral over the unit hypercube
    double jacobian = (1+rho*rho)/pow(1-rho*rho, 4) * 2*PI*PI * rho*rho * sin(PI*u);

    // Quantities occuring in the integrand
    double r_a = sqrt(pow(x_trans[0] - pos_a[0], 2) + pow(x_trans[1] - pos_a[1], 2) + pow(x_trans[2] - pos_a[2], 2));
    double r_b = sqrt(pow(x_trans[0] - pos_b[0], 2) + pow(x_trans[1] - pos_b[1], 2) + pow(x_trans[2] - pos_b[2], 2));
    double r_c = sqrt(pow(x_trans[0] - pos_c[0], 2) + pow(x_trans[1] - pos_c[1], 2) + pow(x_trans[2] - pos_c[2], 2));
    double r_d = sqrt(pow(x_trans[0] - pos_d[0], 2) + pow(x_trans[1] - pos_d[1], 2) + pow(x_trans[2] - pos_d[2], 2));
    double r_ab = sqrt(pow(pos_a[0] - pos_b[0], 2) + pow(pos_a[1] - pos_b[1], 2) + pow(pos_a[2] - pos_b[2], 2));
    double s_ab = r_a + r_b + r_ab;
    double n_a[] = {(x_trans[0] - pos_a[0])/r_a, (x_trans[1] - pos_a[1])/r_a, (x_trans[2] - pos_a[2])/r_a};
    double n_b[] = {(x_trans[0] - pos_b[0])/r_b, (x_trans[1] - pos_b[1])/r_b, (x_trans[2] - pos_b[2])/r_b};
    double n_c[] = {(x_trans[0] - pos_c[0])/r_c, (x_trans[1] - pos_c[1])/r_c, (x_trans[2] - pos_c[2])/r_c};
    double n_d[] = {(x_trans[0] - pos_d[0])/r_d, (x_trans[1] - pos_d[1])/r_d, (x_trans[2] - pos_d[2])/r_d};
    double n_ab[] = {(pos_a[0] - pos_b[0])/r_ab, (pos_a[1] - pos_b[1])/r_ab, (pos_a[2] - pos_b[2])/r_ab};
    double n_a_dot_n_c = 0.0;
    double n_b_dot_n_d = 0.0;
    double n_c_dot_n_d = 0.0;
    double n_c_dot_n_ab = 0.0;
    double n_d_dot_n_ab = 0.0;
    for (int i = 0; i < 3; i++) {
        n_a_dot_n_c += n_a[i] * n_c[i];
        n_b_dot_n_d += n_b[i] * n_d[i];
        n_c_dot_n_d += n_c[i] * n_d[i];
        n_c_dot_n_ab += n_c[i] * n_ab[i];
        n_d_dot_n_ab += n_d[i] * n_ab[i];
    }

    ff[0] = 1/(r_c * r_c * r_d * r_d) * ( (n_c_dot_n_ab - n_a_dot_n_c) * (n_d_dot_n_ab + n_b_dot_n_d) / (s_ab * s_ab) - (n_c_dot_n_d - n_c_dot_n_ab * n_d_dot_n_ab) / (r_ab * s_ab) ) ;
    ff[0] *= jacobian;
    return 0;
}


double H0PN(double* w, struct ode_params* params) {
     int num_bodies = params->num_bodies;
     int num_dim = params->num_dim;
     int array_half = num_dim * num_bodies;
     double rel_dist2, p2;
     double H = 0;

     for (int a = 0; a < num_bodies; a++) {
          // Kinetic energy
          p2 = 0.0;
          for (int i = 0; i < num_dim; i++)
               p2 += pow(w[array_half + num_dim * a + i], 2);
          H += p2/(2 * params->masses[a]);

          // Potential energy
          for (int b = a+1; b < num_bodies; b++) {
               rel_dist2 = 0.0;
               for (int i = 0; i < num_dim; i++)
                    rel_dist2 += pow(w[num_dim * a + i] - w[num_dim * b + i], 2);
               H -= params->masses[a] * params->masses[b] / sqrt(rel_dist2);
          }
     }
     return H;
}


double H1PN(double* w, struct ode_params* params) {
     int num_bodies = params->num_bodies;
     int num_dim = params->num_dim;
     int a, b, c, i;
     double m_a, m_b, m_c;
     double pa_dot_pa, pb_dot_pb, pa_dot_pb;
     double dx, r_ab, r_ac, ni, na_dot_pa, na_dot_pb;
     double p_ai, p_bi, dx_ac;
     int array_half = num_dim * num_bodies;
     double H = 0.0;

     // Compute kinetic energy and potential energy
     for (a = 0; a < num_bodies; a++) {
          m_a = params->masses[a];
          pa_dot_pa = 0.0;

          for (i = 0; i < num_dim; i++) {
               p_ai = w[array_half + a * num_dim + i];
               pa_dot_pa += p_ai * p_ai;
          }

          H += -0.125 * m_a * pow(pa_dot_pa / (m_a * m_a), 2);

          for (b = 0; b < num_bodies; b++) {
               if (b == a) continue;

               m_b = params->masses[b];
               r_ab = 0.0;
               pb_dot_pb = 0.0;
               pa_dot_pb = 0.0;
               na_dot_pa = 0.0;
               na_dot_pb = 0.0;

               for (i = 0; i < num_dim; i++) {
                    dx = w[a * num_dim + i] - w[b * num_dim + i];
                    r_ab += dx * dx;
               }

               r_ab = sqrt(r_ab);

               for (i = 0; i < num_dim; i++) {
                    p_ai = w[array_half + a * num_dim + i];
                    p_bi = w[array_half + b * num_dim + i];

                    dx = w[a * num_dim + i] - w[b * num_dim + i];
                    ni = dx / r_ab;

                    pb_dot_pb += p_bi * p_bi;
                    pa_dot_pb += p_ai * p_bi;
                    na_dot_pa += ni * p_ai;
                    na_dot_pb += ni * p_bi;
               }

               H += -0.25 * m_a * m_b / r_ab * (6 * pa_dot_pa / (m_a * m_a) - 7 * pa_dot_pb / (m_a * m_b) - 
                    (na_dot_pa * na_dot_pb) / (m_a * m_b));

               if (num_bodies > 0) {
                    for (c = 0; c < num_bodies; c++) {
                         if (c == a) continue;

                         m_c = params->masses[c];
                         r_ac = 0.0;

                         for (i = 0; i < num_dim; i++) {
                         dx_ac = w[a * num_dim + i] - w[c * num_dim + i];
                         r_ac += dx_ac * dx_ac;
                         }
                         r_ac = sqrt(r_ac);

                         H += 0.5 * m_a * m_b * m_c / (r_ab * r_ac);
                    }
               }
          }
     }
     return H;
}


double H2PN_threebody(double* w, struct ode_params* params, int p_flag) {
     int num_bodies = params->num_bodies;
     int num_dim = params->num_dim;
     int array_half = num_bodies * num_dim; 

     // Masses
     double* m;
     allocate_vector(&m, num_bodies);
     for (int a = 0; a < num_bodies; a++)
          m[a] = params->masses[a];
     
     // Momenta
     double** p;
     allocate_2d_array(&p, num_bodies, num_dim);
     for (int a = 0; a < num_bodies; a++)
          for (int i = 0; i < num_dim; i++)
               p[a][i] = w[array_half + a * num_dim + i];
     
     // Relative positions and distances
     double*** x_rel, **r, ***n, *n_ab_ac, *n_ab_cb;
     allocate_3d_array(&x_rel, num_bodies, num_bodies, num_dim);
     allocate_3d_array(&n, num_bodies, num_bodies, num_dim);
     allocate_vector(&n_ab_ac, num_dim);
     allocate_vector(&n_ab_cb, num_dim);
     allocate_2d_array(&r, num_bodies, num_bodies);
     for (int a = 0; a < num_bodies; a++) {
          for (int b = a; b < num_bodies; b++) {
               for (int i = 0; i < num_dim; i++){
                    x_rel[a][b][i] = w[a * num_dim + i] - w[b * num_dim + i];
                    x_rel[b][a][i] = -x_rel[a][b][i];
               } 
               r[a][b] = norm(x_rel[a][b], num_dim);
               r[b][a] = r[a][b];
               for (int i = 0; i < num_dim; i++){
                    if (a == b){
                         n[a][b][i] = 0.0;
                         n[b][a][i] = 0.0;
                    } else {
                         n[a][b][i] = x_rel[a][b][i] / r[a][b];
                         n[b][a][i] = -n[a][b][i];
                    }
               }
          }
     }

     // Compute H
     double H = 0.0;
     for (int a = 0; a < num_bodies; a++)
          H += 0.0625 * m[a] * pow(dot_product(p[a], p[a], num_dim) / (m[a] * m[a]), 3);

     for (int a = 0; a < num_bodies; a++) {
          for (int b = 0; b < num_bodies; b++) {
               if (b != a){
                    H += 0.0625 * 1 / r[a][b] * (10 * m[a] * m[b] * pow(dot_product(p[a], p[a], num_dim) / (m[a] * m[a]), 2) 
                         - 11 * dot_product(p[a], p[a], num_dim) * dot_product(p[b], p[b], num_dim) / (m[a] * m[b])
                         - 2 * dot_product(p[a], p[a], num_dim) * dot_product(p[a], p[a], num_dim) / (m[a] * m[b])
                         + 10 * dot_product(p[a], p[a], num_dim) * pow(dot_product(n[a][b], p[b], num_dim), 2) / (m[a] * m[b]) 
                         - 12 * dot_product(p[a], p[b], num_dim) * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim) / (m[a] * m[b])
                         - 3 * pow(dot_product(n[a][b], p[a], num_dim), 2) * pow(dot_product(n[a][b], p[b], num_dim), 2) / (m[a] * m[b]));
                    H += 0.25 * m[a] * m[a] * m[b] / (r[a][b] * r[a][b]) * (dot_product(p[a], p[a], num_dim) / (m[a] * m[a])
                         + dot_product(p[b], p[b], num_dim) / (m[b] * m[b])
                         - 2 * dot_product(p[a], p[b], num_dim) / (m[a] * m[b]));
                    if (!p_flag) {
                         H += 0.5 * m[a] * m[a] * m[a] * m[b] / (r[a][b] * r[a][b] * r[a][b]);
                         H += 0.375 * m[a] * m[a] * m[b] * m[b] / (r[a][b] * r[a][b] * r[a][b]);
                         H += -0.25 * m[a] * m[a] * m[b] * m[b] / (r[a][b] * r[a][b] * r[a][b]);
                    }
               }
          } 
     }

     if (num_bodies > 0){
          for (int a = 0; a < num_bodies; a++) {
               for (int b = 0; b < num_bodies; b++) {
                    for (int c = 0; c < num_bodies; c++) {
                         if (b != a && c != a){
                              H += 0.125 * m[a] * m[b] * m[c] / (r[a][b] * r[a][c]) * (18 * dot_product(p[a], p[a], num_dim) / (m[a] * m[a]) 
                                   + 14 * dot_product(p[b], p[b], num_dim) / (m[b] * m[b]) 
                                   - 2 * pow(dot_product(n[a][b], p[b], num_dim), 2) / (m[b] * m[b])
                                   - 50 * dot_product(p[a], p[b], num_dim) / (m[a] * m[b]) 
                                   + 17 * dot_product(p[b], p[c], num_dim) / (m[b] * m[c]) 
                                   - 14 * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim) / (m[a] * m[b]) 
                                   + 14 * dot_product(n[a][b], p[b], num_dim) * dot_product(n[a][b], p[c], num_dim) / (m[b] * m[c]) 
                                   + dot_product(n[a][b], n[a][c], num_dim) * dot_product(n[a][b], p[b], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[b] * m[c]));
                              H += 0.125 * m[a] * m[b] * m[c] / (r[a][b] * r[a][b]) * (2 * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[a] * m[c])
                                   + 2 * dot_product(n[a][b], p[b], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[a] * m[c])
                                   + 5 * dot_product(n[a][b], n[a][c], num_dim) * dot_product(p[c], p[c], num_dim) / (m[c] * m[c])
                                   - dot_product(n[a][b], n[a][c], num_dim) * dot_product(n[a][c], p[c], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[c] * m[c])
                                   - 14 * dot_product(n[a][b], p[c], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[c] * m[c]));
                              if (!p_flag){
                                   H += -0.25 * m[a] * m[b] * m[c] * m[c] / (r[a][b] * r[a][c] * r[a][c]);
                                   H += -0.75 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][b] * r[a][c]);
                              }
                         }
                         if (b != a && c != a && c != b) {
                              for (int i = 0; i < num_dim; i++) {
                                   n_ab_ac[i] = n[a][b][i] + n[a][c][i];
                                   n_ab_cb[i] = n[a][b][i] + n[c][b][i];
                              }
                              H += 0.5 * m[a] * m[b] * m[c] / pow(r[a][b] + r[b][c] + r[c][a], 2) * (
                                   8 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[c], num_dim) / (m[a] * m[c])
                                   - 16 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[a], num_dim) / (m[a] * m[c])
                                   + 3 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[b], num_dim) / (m[a] * m[b])
                                   + 4 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[c], num_dim) / (m[c] * m[c])
                                   + dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[a], num_dim) / (m[a] * m[a]));
                              H += 0.5 * m[a] * m[b] * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * r[a][b]) * (
                                   8 * (dot_product(p[a], p[c], num_dim) - dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[c], num_dim)) / (m[a] * m[c])
                                   - 3 * (dot_product(p[a], p[b], num_dim) - dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim)) / (m[a] * m[b])
                                   - 4 * (dot_product(p[c], p[c], num_dim) - dot_product(n[a][b], p[c], num_dim) * dot_product(n[a][b], p[c], num_dim)) / (m[c] * m[c])
                                   - (dot_product(p[a], p[a], num_dim) - dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[a], num_dim)) / (m[a] * m[a]));
                              if (!p_flag) {
                                   H += -0.375 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][c] * r[b][c]);
                                   H += -0.015625 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][b] * r[a][b] * r[a][c] * r[a][c] * r[a][c] * r[b][c]) * (
                                        18 * r[a][b] * r[a][b] * r[a][c] * r[a][c] - 60 * r[a][b] * r[a][b] * r[b][c] * r[b][c]
                                        - 24 * r[a][b] * r[a][b] * r[a][c] * (r[a][b] + r[b][c])
                                        + 60 * r[a][b] * r[a][c] * r[b][c] * r[b][c] + 56 * r[a][b] * r[a][b] * r[a][b] * r[b][c]
                                        - 72 * r[a][b] * r[b][c] * r[b][c] * r[b][c] + 35 * r[b][c] * r[b][c] * r[b][c] * r[b][c] + 6 * r[a][b] * r[a][b] * r[a][b] * r[a][b]);
                              }
                         }
                         if (b != a && c != b)
                              if (!p_flag)
                                   H += -0.5 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][b] * r[b][c]);
                    }
               }
          }
     }
     free_vector(m);
     free_2d_array(p, num_bodies);
     free_3d_array(x_rel, num_bodies, num_bodies);
     free_3d_array(n, num_bodies, num_bodies);
     free_vector(n_ab_ac);
     free_vector(n_ab_cb);
     free_2d_array(r, num_bodies);
     return H;
}


double H2PN_nbody(double* w, struct ode_params* params, int p_flag) {
     int num_bodies = params->num_bodies;
     int num_dim = params->num_dim;
     int array_half = num_bodies * num_dim; 

     // Masses
     double* m;
     allocate_vector(&m, num_bodies);
     for (int a = 0; a < num_bodies; a++)
          m[a] = params->masses[a];
     
     // Momenta
     double** p;
     allocate_2d_array(&p, num_bodies, num_dim);
     for (int a = 0; a < num_bodies; a++)
          for (int i = 0; i < num_dim; i++)
               p[a][i] = w[array_half + a * num_dim + i];

     // Positions
     IntegralParams integral_params;
     integral_params.pos = malloc(num_bodies * sizeof(double*));
     for (int i = 0; i < num_bodies; i++)
          integral_params.pos[i] = malloc(num_dim * sizeof(double));
     
     // Relative positions and distances
     double*** x_rel, **r, ***n, *n_ab_ac, *n_ab_cb;
     allocate_3d_array(&x_rel, num_bodies, num_bodies, num_dim);
     allocate_3d_array(&n, num_bodies, num_bodies, num_dim);
     allocate_vector(&n_ab_ac, num_dim);
     allocate_vector(&n_ab_cb, num_dim);
     allocate_2d_array(&r, num_bodies, num_bodies);
     for (int a = 0; a < num_bodies; a++) {
          for (int b = a; b < num_bodies; b++) {
               for (int i = 0; i < num_dim; i++){
                    x_rel[a][b][i] = w[a * num_dim + i] - w[b * num_dim + i];
                    x_rel[b][a][i] = -x_rel[a][b][i];
               } 
               r[a][b] = norm(x_rel[a][b], num_dim);
               r[b][a] = r[a][b];
               for (int i = 0; i < num_dim; i++){
                    if (a == b){
                         n[a][b][i] = 0.0;
                         n[b][a][i] = 0.0;
                    } else {
                         n[a][b][i] = x_rel[a][b][i] / r[a][b];
                         n[b][a][i] = -n[a][b][i];
                    }
               }
          }
     }

     // Compute H
     double H = 0.0;
     for (int a = 0; a < num_bodies; a++)
          H += 0.0625 * m[a] * pow(dot_product(p[a], p[a], num_dim) / (m[a] * m[a]), 3);

     for (int a = 0; a < num_bodies; a++) {
          for (int b = 0; b < num_bodies; b++) {
               if (b != a){
                    H += 0.0625 * 1 / r[a][b] * (10 * m[a] * m[b] * pow(dot_product(p[a], p[a], num_dim) / (m[a] * m[a]), 2) 
                         - 11 * dot_product(p[a], p[a], num_dim) * dot_product(p[b], p[b], num_dim) / (m[a] * m[b])
                         - 2 * dot_product(p[a], p[a], num_dim) * dot_product(p[a], p[a], num_dim) / (m[a] * m[b])
                         + 10 * dot_product(p[a], p[a], num_dim) * pow(dot_product(n[a][b], p[b], num_dim), 2) / (m[a] * m[b]) 
                         - 12 * dot_product(p[a], p[b], num_dim) * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim) / (m[a] * m[b])
                         - 3 * pow(dot_product(n[a][b], p[a], num_dim), 2) * pow(dot_product(n[a][b], p[b], num_dim), 2) / (m[a] * m[b]));
                    H += 0.25 * m[a] * m[a] * m[b] / (r[a][b] * r[a][b]) * (dot_product(p[a], p[a], num_dim) / (m[a] * m[a])
                         + dot_product(p[b], p[b], num_dim) / (m[b] * m[b])
                         - 2 * dot_product(p[a], p[b], num_dim) / (m[a] * m[b]));
                    if (!p_flag)
                         H += -0.25 * m[a] * m[a] * m[b] * m[b] / (r[a][b] * r[a][b] * r[a][b]);
               }
          } 
     }

     for (int a = 0; a < num_bodies; a++) {
          for (int b = 0; b < num_bodies; b++) {
               for (int c = 0; c < num_bodies; c++) {
                    if (b != a && c != a){
                         H += 0.125 * m[a] * m[b] * m[c] / (r[a][b] * r[a][c]) * (18 * dot_product(p[a], p[a], num_dim) / (m[a] * m[a]) 
                              + 14 * dot_product(p[b], p[b], num_dim) / (m[b] * m[b]) 
                              - 2 * pow(dot_product(n[a][b], p[b], num_dim), 2) / (m[b] * m[b])
                              - 50 * dot_product(p[a], p[b], num_dim) / (m[a] * m[b]) 
                              + 17 * dot_product(p[b], p[c], num_dim) / (m[b] * m[c]) 
                              - 14 * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim) / (m[a] * m[b]) 
                              + 14 * dot_product(n[a][b], p[b], num_dim) * dot_product(n[a][b], p[c], num_dim) / (m[b] * m[c]) 
                              + dot_product(n[a][b], n[a][c], num_dim) * dot_product(n[a][b], p[b], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[b] * m[c]));
                         H += 0.125 * m[a] * m[b] * m[c] / (r[a][b] * r[a][b]) * (2 * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[a] * m[c])
                              + 2 * dot_product(n[a][b], p[b], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[a] * m[c])
                              + 5 * dot_product(n[a][b], n[a][c], num_dim) * dot_product(p[c], p[c], num_dim) / (m[c] * m[c])
                              - dot_product(n[a][b], n[a][c], num_dim) * dot_product(n[a][c], p[c], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[c] * m[c])
                              - 14 * dot_product(n[a][b], p[c], num_dim) * dot_product(n[a][c], p[c], num_dim) / (m[c] * m[c]));
                    }
                    if (b != a && c != a && c != b) {
                         for (int i = 0; i < num_dim; i++) {
                              n_ab_ac[i] = n[a][b][i] + n[a][c][i];
                              n_ab_cb[i] = n[a][b][i] + n[c][b][i];
                         }
                         H += 0.5 * m[a] * m[b] * m[c] / pow(r[a][b] + r[b][c] + r[c][a], 2) * (
                              8 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[c], num_dim) / (m[a] * m[c])
                              - 16 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[a], num_dim) / (m[a] * m[c])
                              + 3 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[b], num_dim) / (m[a] * m[b])
                              + 4 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[c], num_dim) / (m[c] * m[c])
                              + dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[a], num_dim) / (m[a] * m[a]));
                         H += 0.5 * m[a] * m[b] * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * r[a][b]) * (
                              8 * (dot_product(p[a], p[c], num_dim) - dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[c], num_dim)) / (m[a] * m[c])
                              - 3 * (dot_product(p[a], p[b], num_dim) - dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim)) / (m[a] * m[b])
                              - 4 * (dot_product(p[c], p[c], num_dim) - dot_product(n[a][b], p[c], num_dim) * dot_product(n[a][b], p[c], num_dim)) / (m[c] * m[c])
                              - (dot_product(p[a], p[a], num_dim) - dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[a], num_dim)) / (m[a] * m[a]));
                         if (!p_flag)
                              H += -0.015625 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][b] * r[a][b] * r[a][c] * r[a][c] * r[a][c] * r[b][c]) * (
                                   18 * r[a][b] * r[a][b] * r[a][c] * r[a][c] - 60 * r[a][b] * r[a][b] * r[b][c] * r[b][c]
                                   - 24 * r[a][b] * r[a][b] * r[a][c] * (r[a][b] + r[b][c])
                                   + 60 * r[a][b] * r[a][c] * r[b][c] * r[b][c] + 56 * r[a][b] * r[a][b] * r[a][b] * r[b][c]
                                   - 72 * r[a][b] * r[b][c] * r[b][c] * r[b][c] + 35 * r[b][c] * r[b][c] * r[b][c] * r[b][c] + 6 * r[a][b] * r[a][b] * r[a][b] * r[a][b]);
                    }
               }
          }
     }

     // G^3 terms
     if (!p_flag) {
          cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];
          int neval, fail, nregions;
          for (int a = 0; a < num_bodies; a++) {
               for (int b = 0; b < num_bodies; b++) {
                    for (int c = 0; c < num_bodies; c++) {
                         for (int d = 0; d < num_bodies; d++) {
                              if (b != a && c != b && d != c) {
                                   H += - 0.375 * m[a] * m[b] * m[c] * m[d] / (r[a][b] * r[b][c] * r[c][d]);
                              }
                              if (b != a && c != a && d != a) {
                                   H += - 0.25 * m[a] * m[b] * m[c] * m[d] / (r[a][b] * r[a][c] * r[a][d]);
                              }
                              if (b != a && c != a && c != b && d != a && d != b && d != c) {
                                   H += - 0.0625 * m[a] * m[b] * m[c] * m[d] / (pow(r[a][b], 3) * pow(r[c][d], 3) * pow(r[a][d], 3) * pow(r[b][c], 3)) * (
                                        4 * pow(r[a][b], 3) * pow(r[b][c], 3) * pow(r[c][d], 2) * pow(r[a][d], 2) / r[b][d]
                                        -6* pow(r[b][c], 3) * pow(r[a][d], 2) * (pow(r[a][b], 2) * pow(r[c][d], 2) + pow(r[a][d], 2) * (pow(r[a][d], 2) + pow(r[b][c], 2) - pow(r[a][c], 2) - pow(r[b][d], 2)))
                                        - pow(r[a][b], 2) * (pow(r[b][d], 2) - pow(r[b][c], 2) - pow(r[c][d], 2)) * (2 * r[a][b] * pow(r[a][d], 3) * pow(r[b][c], 2) / (r[a][c] + r[b][c] + r[a][b]) - pow(r[a][d], 2) * pow(r[b][c], 2) - r[a][b] * pow(r[c][d], 2) * (pow(r[a][c], 2) - pow(r[a][d], 2) - pow(r[c][d], 2))) );

                                   // ln integral
                                   for (int i = 0; i < num_dim; i++){
                                        integral_params.pos[0][i] = w[a * num_dim + i];
                                        integral_params.pos[1][i] = w[b * num_dim + i];
                                        integral_params.pos[2][i] = w[c * num_dim + i];
                                        integral_params.pos[3][i] = w[d * num_dim + i];
                                   }
                                   Cuhre(NDIM, NCOMP, ln_integral, &integral_params, NVEC,
                                        EPSREL, EPSABS, 0,
                                        MINEVAL, MAXEVAL, KEY,
                                        NULL, NULL,
                                        &nregions, &neval, &fail, integral, error, prob);
                                   H += 0.25/PI * m[a] * m[b] * m[c] * m[d] * integral[0];
                              }
                         }
                    }
               }
          }
     }
     free_vector(m);
     free_2d_array(p, num_bodies);
     free_3d_array(x_rel, num_bodies, num_bodies);
     free_3d_array(n, num_bodies, num_bodies);
     free_vector(n_ab_ac);
     free_vector(n_ab_cb);
     free_2d_array(r, num_bodies);
     return H;
}


void update_eom_hamiltonian(double *w, double *dwdt, double (*hamiltonian)(double*, struct ode_params*, int p_flag), double h, struct ode_params* params) {
     int array_half = params->num_dim * params->num_bodies;
     double w_copy[2*array_half];
     double dHdw[2*array_half];
     double forward_value, backward_value;
     int p_flag;

     // Copy original array to w_copy
     for (int i = 0; i < 2 * array_half; i++) {
          w_copy[i] = w[i];
     }

     for (int i = 0; i < 2 * array_half; i++) {
          p_flag = (i < array_half) ? 0 : 1;
               
          // Perturb the i-th variable forward
          w_copy[i] = w[i] + h;
          forward_value = hamiltonian(w_copy, params, p_flag);

          // Perturb the i-th variable backward
          w_copy[i] = w[i] - h;
          backward_value = hamiltonian(w_copy, params, p_flag);

          // Restore original value
          w_copy[i] = w[i];

          // Compute the partial derivative
          dHdw[i] = (forward_value - backward_value) / (2 * h);
     }

     // Compute dwdt
     for (int i = 0; i < 2 * array_half; i++) {
          if (i < array_half)
               dwdt[i] += dHdw[i + array_half];
          else
               dwdt[i] += -dHdw[i - array_half];
     }
}