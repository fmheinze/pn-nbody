#include <math.h>
#include "pn_eom.h"
#include "pn_eom_hamiltonians.h"
#include "utils.h"


double H0PN(double* w, struct ode_params* params) {
     int array_half = params->dim * params->num_bodies;
     double rel_dist2, p2;
     double H = 0;

     for (int a = 0; a < params->num_bodies; a++) {
          // Kinetic energy
          p2 = 0.0;
          for (int i = 0; i < params->dim; i++)
               p2 += pow(w[array_half + params->dim * a + i], 2);
          H += p2/(2 * params->masses[a]);

          // Potential energy
          for (int b = a+1; b < params->num_bodies; b++) {
               rel_dist2 = 0.0;
               for (int i = 0; i < params->dim; i++)
                    rel_dist2 += pow(w[params->dim * a + i] - w[params->dim * b + i], 2);
               H -= params->masses[a] * params->masses[b] / sqrt(rel_dist2);
          }
     }
     return H;
}


double H1PN(double* w, struct ode_params* params) {
    int a, b, c, i;
    double m_a, m_b, m_c;
    double pa_dot_pa, pb_dot_pb, pa_dot_pb;
    double dx, r_ab, r_ac, ni, na_dot_pa, na_dot_pb;
    double p_ai, p_bi, dx_ac;
    int array_half = params->dim * params->num_bodies;
    double H = 0.0;

    // Compute kinetic energy and potential energy
    for (a = 0; a < params->num_bodies; a++) {
        m_a = params->masses[a];
        pa_dot_pa = 0.0;

        for (i = 0; i < params->dim; i++) {
            p_ai = w[array_half + a * params->dim + i];
            pa_dot_pa += p_ai * p_ai;
        }

        H += -0.125 * m_a * pow(pa_dot_pa / (m_a * m_a), 2);

        for (b = 0; b < params->num_bodies; b++) {
            if (b == a) continue;

            m_b = params->masses[b];
            r_ab = 0.0;
            pb_dot_pb = 0.0;
            pa_dot_pb = 0.0;
            na_dot_pa = 0.0;
            na_dot_pb = 0.0;

            for (i = 0; i < params->dim; i++) {
                dx = w[a * params->dim + i] - w[b * params->dim + i];
                r_ab += dx * dx;
            }

            r_ab = sqrt(r_ab);

            for (i = 0; i < params->dim; i++) {
                p_ai = w[array_half + a * params->dim + i];
                p_bi = w[array_half + b * params->dim + i];

                dx = w[a * params->dim + i] - w[b * params->dim + i];
                ni = dx / r_ab;

                pb_dot_pb += p_bi * p_bi;
                pa_dot_pb += p_ai * p_bi;
                na_dot_pa += ni * p_ai;
                na_dot_pb += ni * p_bi;
            }

            H += -0.25 * m_a * m_b / r_ab * (6 * pa_dot_pa / (m_a * m_a) - 7 * pa_dot_pb / (m_a * m_b) - 
                  (na_dot_pa * na_dot_pb) / (m_a * m_b));

            if (params->num_bodies > 0) {
                for (c = 0; c < params->num_bodies; c++) {
                    if (c == a) continue;

                    m_c = params->masses[c];
                    r_ac = 0.0;

                    for (i = 0; i < params->dim; i++) {
                        dx_ac = w[a * params->dim + i] - w[c * params->dim + i];
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


double H2PN(double* w, struct ode_params* params) {
     int array_half = params->num_bodies * params->dim; 

     // Masses
     double* m;
     allocate_vector(&m, params->num_bodies);
     for (int a = 0; a < params->num_bodies; a++)
          m[a] = params->masses[a];
     
     // Momenta
     double** p;
     allocate_array(&p, params->num_bodies, params->dim);
     for (int a = 0; a < params->num_bodies; a++)
          for (int i = 0; i < params->dim; i++)
               p[a][i] = w[array_half + a * params->dim + i];
     
     // Relative positions and distances
     double*** x_rel, **r, ***n, *n_ab_ac, *n_ab_cb;
     allocate_3d_array(&x_rel, params->num_bodies, params->num_bodies, params->dim);
     allocate_3d_array(&n, params->num_bodies, params->num_bodies, params->dim);
     allocate_vector(&n_ab_ac, params->dim);
     allocate_vector(&n_ab_cb, params->dim);
     allocate_array(&r, params->num_bodies, params->num_bodies);
     for (int a = 0; a < params->num_bodies; a++) {
          for (int b = a; b < params->num_bodies; b++) {
               for (int i = 0; i < params->dim; i++){
                    x_rel[a][b][i] = w[a * params->dim + i] - w[b * params->dim + i];
                    x_rel[b][a][i] = -x_rel[a][b][i];
               } 
               r[a][b] = norm(x_rel[a][b], params->dim);
               r[b][a] = r[a][b];
               for (int i = 0; i < params->dim; i++){
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
     for (int a = 0; a < params->num_bodies; a++)
          H += 0.0625 * m[a] * pow(dot_product(p[a], p[a], params->dim) / (m[a] * m[a]), 3);

     for (int a = 0; a < params->num_bodies; a++) {
          for (int b = 0; b < params->num_bodies; b++) {
               if (b != a){
                    H += 0.0625 * 1 / r[a][b] * (10 * m[a] * m[b] * pow(dot_product(p[a], p[a], params->dim) / (m[a] * m[a]), 2) 
                         - 11 * dot_product(p[a], p[a], params->dim) * dot_product(p[b], p[b], params->dim) / (m[a] * m[b])
                         - 2 * dot_product(p[a], p[a], params->dim) * dot_product(p[a], p[a], params->dim) / (m[a] * m[b])
                         + 10 * dot_product(p[a], p[a], params->dim) * pow(dot_product(n[a][b], p[b], params->dim), 2) / (m[a] * m[b]) 
                         - 12 * dot_product(p[a], p[b], params->dim) * dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][b], p[b], params->dim) / (m[a] * m[b])
                         - 3 * pow(dot_product(n[a][b], p[a], params->dim), 2) * pow(dot_product(n[a][b], p[b], params->dim), 2) / (m[a] * m[b]));
                    H += 0.25 * m[a] * m[a] * m[b] / (r[a][b] * r[a][b]) * (dot_product(p[a], p[a], params->dim) / (m[a] * m[a])
                         + dot_product(p[b], p[b], params->dim) / (m[b] * m[b])
                         - 2 * dot_product(p[a], p[b], params->dim) / (m[a] * m[b]));
                    H += 0.5 * m[a] * m[a] * m[a] * m[b] / (r[a][b] * r[a][b] * r[a][b]);
                    H += 0.375 * m[a] * m[a] * m[b] * m[b] / (r[a][b] * r[a][b] * r[a][b]);
                    H += -0.25 * m[a] * m[a] * m[b] * m[b] / (r[a][b] * r[a][b] * r[a][b]);
               }
          } 
     }

     if (params->num_bodies > 0){
          for (int a = 0; a < params->num_bodies; a++) {
               for (int b = 0; b < params->num_bodies; b++) {
                    for (int c = 0; c < params->num_bodies; c++) {
                         if (b != a && c != a){
                              H += 0.125 * m[a] * m[b] * m[c] / (r[a][b] * r[a][c]) * (18 * dot_product(p[a], p[a], params->dim) / (m[a] * m[a]) 
                                   + 14 * dot_product(p[b], p[b], params->dim) / (m[b] * m[b]) 
                                   - 2 * pow(dot_product(n[a][b], p[b], params->dim), 2) / (m[b] * m[b])
                                   - 50 * dot_product(p[a], p[b], params->dim) / (m[a] * m[b]) 
                                   + 17 * dot_product(p[b], p[c], params->dim) / (m[b] * m[c]) 
                                   - 14 * dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][b], p[b], params->dim) / (m[a] * m[b]) 
                                   + 14 * dot_product(n[a][b], p[b], params->dim) * dot_product(n[a][b], p[c], params->dim) / (m[b] * m[c]) 
                                   + dot_product(n[a][b], n[a][c], params->dim) * dot_product(n[a][b], p[b], params->dim) * dot_product(n[a][c], p[c], params->dim) / (m[b] * m[c]));
                              H += 0.125 * m[a] * m[b] * m[c] / (r[a][b] * r[a][b]) * (2 * dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][c], p[c], params->dim) / (m[a] * m[c])
                                   + 2 * dot_product(n[a][b], p[b], params->dim) * dot_product(n[a][c], p[c], params->dim) / (m[a] * m[c])
                                   + 5 * dot_product(n[a][b], n[a][c], params->dim) * dot_product(p[c], p[c], params->dim) / (m[c] * m[c])
                                   - dot_product(n[a][b], n[a][c], params->dim) * dot_product(n[a][c], p[c], params->dim) * dot_product(n[a][c], p[c], params->dim) / (m[c] * m[c])
                                   - 14 * dot_product(n[a][b], p[c], params->dim) * dot_product(n[a][c], p[c], params->dim) / (m[c] * m[c]));
                              H += -0.25 * m[a] * m[b] * m[c] * m[c] / (r[a][b] * r[a][c] * r[a][c]);
                              H += -0.75 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][b] * r[a][c]);
                         }
                         if (b != a && c != a && c != b) {
                              for (int i = 0; i < params->dim; i++) {
                                   n_ab_ac[i] = n[a][b][i] + n[a][c][i];
                                   n_ab_cb[i] = n[a][b][i] + n[c][b][i];
                              }
                              H += 0.5 * m[a] * m[b] * m[c] / pow(r[a][b] + r[b][c] + r[c][a], 2) * (
                                   8 * dot_product(n_ab_ac, p[a], params->dim) * dot_product(n_ab_cb, p[c], params->dim) / (m[a] * m[c])
                                   - 16 * dot_product(n_ab_ac, p[c], params->dim) * dot_product(n_ab_cb, p[a], params->dim) / (m[a] * m[c])
                                   + 3 * dot_product(n_ab_ac, p[a], params->dim) * dot_product(n_ab_cb, p[b], params->dim) / (m[a] * m[b])
                                   + 4 * dot_product(n_ab_ac, p[c], params->dim) * dot_product(n_ab_cb, p[c], params->dim) / (m[c] * m[c])
                                   + dot_product(n_ab_ac, p[a], params->dim) * dot_product(n_ab_cb, p[a], params->dim) / (m[a] * m[a]));
                              H += 0.5 * m[a] * m[b] * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * r[a][b]) * (
                                   8 * (dot_product(p[a], p[c], params->dim) - dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][b], p[c], params->dim)) / (m[a] * m[c])
                                   - 3 * (dot_product(p[a], p[b], params->dim) - dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][b], p[b], params->dim)) / (m[a] * m[b])
                                   - 4 * (dot_product(p[c], p[c], params->dim) - dot_product(n[a][b], p[c], params->dim) * dot_product(n[a][b], p[c], params->dim)) / (m[c] * m[c])
                                   - (dot_product(p[a], p[a], params->dim) - dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][b], p[a], params->dim)) / (m[a] * m[a]));
                              H += -0.375 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][c] * r[b][c]);
                              H += 0.015625 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][b] * r[a][b] * r[a][c] * r[a][c] * r[a][c] * r[b][c]) * (
                                   18 * r[a][b] * r[a][b] * r[a][c] * r[a][c] - 60 * r[a][b] * r[a][b] * r[b][c] * r[b][c]
                                   - 24 * r[a][b] * r[a][b] * r[a][c] * (r[a][b] + r[b][c])
                                   + 60 * r[a][b] * r[a][c] * r[b][c] * r[b][c] + 56 * r[a][b] * r[a][b] * r[a][b] * r[b][c]
                                   - 72 * r[a][b] * r[b][c] * r[b][c] * r[b][c] + 35 * r[b][c] * r[b][c] * r[b][c] * r[b][c] + 6 * r[a][b] * r[a][b] * r[a][b] * r[a][b]);
                         }
                         if (b != a && c != b) {
                              H += -0.5 * m[a] * m[a] * m[b] * m[c] / (r[a][b] * r[a][b] * r[b][c]);
                         }
                    }
               }
          }
     }

     return H;
}


void update_eom_hamiltonian(double *w, double *dwdt, double (*hamiltonian)(double*, struct ode_params*), double h, struct ode_params* params) {
     int array_half = params->dim * params->num_bodies;
     double w_copy[2*array_half];
     double dHdw[2*array_half];
     double forward_value, backward_value;

     // Copy original array to w_copy
     for (int i = 0; i < 2 * array_half; i++) {
          w_copy[i] = w[i];
     }

     for (int i = 0; i < 2 * array_half; i++) {
          // Perturb the i-th variable forward
          w_copy[i] = w[i] + h;
          forward_value = hamiltonian(w_copy, params);

          // Perturb the i-th variable backward
          w_copy[i] = w[i] - h;
          backward_value = hamiltonian(w_copy, params);

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