#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "pn_eom.h"
#include "pn_eom_hamiltonians.h"
#include "utils.h"

#define PI 3.1415926535897932384626433832795
#define INVPI 0.31830988618379067153776752674503
#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif

// Userdata for the 2PN ln-integral
#define NUM_LOCAL_BODIES 4
#define NCOMP_H_DERIV 72
typedef struct {
    double pos[NUM_LOCAL_BODIES][3];
} IntegralParams;

// Unique permutation representatives for the 2PN ln-integral
#define NUM_PERMS 6
static const int perms[NUM_PERMS][NUM_LOCAL_BODIES] = {
    {0,1,2,3},
    {0,1,3,2},
    {0,2,1,3},
    {0,2,3,1},
    {0,3,1,2},
    {0,3,2,1}
};

static inline int role_of_body(int p, int body) {
    for (int r = 0; r < 4; ++r)
        if (perms[p][r] == body) return r;
    return -1;
}

// Integral parameters
#define NVEC 1
#define EPSREL 1e-16
#define EPSABS 1e-20
#define VERBOSE 0
#define SEED 42
#define MINEVAL 500000
#define MAXEVAL 500000
#define KEY 11


static int ln_integral(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
    if (*ndim != 3)
        errorexit("The 2PN ln-integral can only be computed in 3D! Please use num_dim = 3");

    // Userdata
    IntegralParams *params = (IntegralParams*)userdata;

    // New variables and jacobian determinant for the transformation to an integral over the unit hypercube
    double x_trans[3];
    cubareal rho = xx[0];
    cubareal u = xx[1];
    cubareal v = xx[2];
    if (rho > 1 - 1e-10) rho = 1 - 1e-10;

    cubareal PIu = PI * u;
    cubareal twoPIv = 2.0 * PI * v;
    cubareal rho2 = rho * rho;
    cubareal one_minus_rho2 = 1.0 - rho2;
    cubareal factor = rho / one_minus_rho2;
    x_trans[0] = factor * sin(PIu) * cos(twoPIv);
    x_trans[1] = factor * sin(PIu) * sin(twoPIv);
    x_trans[2] = factor * cos(PIu);

    cubareal jacobian = (1.0 + rho2) / (one_minus_rho2 * one_minus_rho2 * one_minus_rho2 * one_minus_rho2) * 2.0 * PI * PI * rho2 * sin(PIu);

    for (int p = 0; p < 6; p++) {
        int A_idx, B_idx, C_idx, D_idx;
        A_idx = perms[p][0];
        B_idx = perms[p][1];
        C_idx = perms[p][2];
        D_idx = perms[p][3];

        double *A = params->pos[A_idx];
        double *B = params->pos[B_idx];
        double *C = params->pos[C_idx];
        double *D = params->pos[D_idx];

        // Quantities occuring in the integrand
        double dx_a[3], dx_b[3], dx_c[3], dx_d[3], ab[3];
        double n_a[3], n_b[3], n_c[3], n_d[3], n_ab[3];
        double r_a, r_b, r_c, r_d, r_ab, s_ab;
        double n_c_dot_n_ab = 0.0;
        double n_c_dot_n_d  = 0.0;
        double n_a_dot_n_c  = 0.0;
        double n_d_dot_n_ab = 0.0;
        double n_b_dot_n_d  = 0.0;

        for (int i = 0; i < 3; ++i) {
            dx_a[i] = x_trans[i] - A[i];
            dx_b[i] = x_trans[i] - B[i];
            dx_c[i] = x_trans[i] - C[i];
            dx_d[i] = x_trans[i] - D[i];
            ab[i] = A[i] - B[i];
        }
        r_a = norm(dx_a, 3);
        r_b = norm(dx_b, 3);
        r_c = norm(dx_c, 3);
        r_d = norm(dx_d, 3);
        r_ab = norm(ab, 3);
        s_ab = r_a + r_b + r_ab;
        for (int i = 0; i < 3; ++i) {
            n_a[i]  = dx_a[i] / r_a;
            n_b[i]  = dx_b[i] / r_b;
            n_c[i]  = dx_c[i] / r_c;
            n_d[i]  = dx_d[i] / r_d;
            n_ab[i] = ab[i]   / r_ab;

            n_c_dot_n_ab += n_c[i] * n_ab[i];
            n_c_dot_n_d  += n_c[i] * n_d[i];
            n_a_dot_n_c  += n_a[i] * n_c[i];
            n_d_dot_n_ab += n_d[i] * n_ab[i];
            n_b_dot_n_d  += n_b[i] * n_d[i];
        }

        ff[p] = (1.0 / (r_c*r_c * r_d*r_d) *
            ( (n_c_dot_n_ab - n_a_dot_n_c) * (n_d_dot_n_ab + n_b_dot_n_d) / (s_ab*s_ab)
            - (n_c_dot_n_d - n_c_dot_n_ab * n_d_dot_n_ab) / (r_ab * s_ab)
            ));
        ff[p] *= jacobian;
    }
    return 0;
}


static int ln_integral_complex(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    if (*ndim != 3)
        errorexit("The 2PN ln-integral can only be computed in 3D! Please use num_dim = 3");

    IntegralParams *params = (IntegralParams*)userdata;
    const double h = 1e-20;  // step size for complex-step differentiation

    // New variables and jacobian determinant for the transformation to an integral over the unit hypercube (can remain real-valued)
    cubareal rho = xx[0];
    cubareal u   = xx[1];
    cubareal v   = xx[2];
    if (rho > 1 - 1e-10) rho = 1 - 1e-10;

    cubareal PIu = PI * u;
    cubareal twoPIv = 2.0 * PI * v;
    cubareal rho2 = rho * rho;
    cubareal one_minus_rho2 = 1.0 - rho2;
    cubareal factor = rho / one_minus_rho2;
    cubareal x_trans[3];
    x_trans[0] = factor * sin(PIu) * cos(twoPIv);
    x_trans[1] = factor * sin(PIu) * sin(twoPIv);
    x_trans[2] = factor * cos(PIu);

    cubareal jacobian = (1.0 + rho2) / (one_minus_rho2 * one_minus_rho2 * one_minus_rho2 * one_minus_rho2) * 2.0 * PI * PI * rho2 * sin(PIu);

    // Make a complex copy of all positions (will re-init per derivative)
    double complex pos_cmplx[3 * NUM_LOCAL_BODIES];
    #define BODY(i) (&pos_cmplx[3*(i)])

    for (int a = 0; a < NUM_LOCAL_BODIES; a++)
        for (int i = 0; i < 3; ++i)
            pos_cmplx[3 * a + i] = (double complex)params->pos[a][i];

    for (int p = 0; p < NUM_PERMS; ++p) {   // Loop over the 6 unique permutations
        for (int body = 0; body < NUM_LOCAL_BODIES; ++body) {
            int r = role_of_body(p, body);  // r = position of body in permutation p, tells us whether to use I_{ab;cd} or I_{cd;ab} later
            for (int axis = 0; axis < 3; ++axis) {

                // Complex-step perturbation for this derivative
                int coord_index = 3 * body + axis;
                pos_cmplx[coord_index] += I * h;

                // Decide which bodies are A,B,C,D in this eval
                int A_idx, B_idx, C_idx, D_idx;
                if (r == 0 || r == 1) {
                    // body is in first pair: use I_{ab;cd}
                    A_idx = perms[p][0];
                    B_idx = perms[p][1];
                    C_idx = perms[p][2];
                    D_idx = perms[p][3];
                } else {
                    // body is in second pair: use I_{cd;ab}
                    A_idx = perms[p][2];
                    B_idx = perms[p][3];
                    C_idx = perms[p][0];
                    D_idx = perms[p][1];
                }
                double complex *A = BODY(A_idx);
                double complex *B = BODY(B_idx);
                double complex *C = BODY(C_idx);
                double complex *D = BODY(D_idx);

                // Geometry
                double complex dx_a[3], dx_b[3], dx_c[3], dx_d[3], ab[3];
                double complex n_a[3], n_b[3], n_c[3], n_d[3], n_ab[3];
                double complex r_a, r_b, r_c, r_d, r_ab, s_ab;
                double complex n_c_dot_n_ab = 0.0;
                double complex n_c_dot_n_d  = 0.0;
                double complex n_a_dot_n_c  = 0.0;
                double complex n_d_dot_n_ab = 0.0;
                double complex n_b_dot_n_d  = 0.0;

                for (int i = 0; i < 3; ++i) {
                    dx_a[i] = x_trans[i] - A[i];
                    dx_b[i] = x_trans[i] - B[i];
                    dx_c[i] = x_trans[i] - C[i];
                    dx_d[i] = x_trans[i] - D[i];
                    ab[i] = A[i] - B[i];
                }

                r_a  = csqrt(dx_a[0]*dx_a[0] + dx_a[1]*dx_a[1] + dx_a[2]*dx_a[2]);
                r_b  = csqrt(dx_b[0]*dx_b[0] + dx_b[1]*dx_b[1] + dx_b[2]*dx_b[2]);
                r_c  = csqrt(dx_c[0]*dx_c[0] + dx_c[1]*dx_c[1] + dx_c[2]*dx_c[2]);
                r_d  = csqrt(dx_d[0]*dx_d[0] + dx_d[1]*dx_d[1] + dx_d[2]*dx_d[2]);
                r_ab = csqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
                s_ab = r_a + r_b + r_ab;

                for (int i = 0; i < 3; ++i) {
                    n_a[i]  = dx_a[i] / r_a;
                    n_b[i]  = dx_b[i] / r_b;
                    n_c[i]  = dx_c[i] / r_c;
                    n_d[i]  = dx_d[i] / r_d;
                    n_ab[i] = ab[i]   / r_ab;

                    n_c_dot_n_ab += n_c[i] * n_ab[i];
                    n_c_dot_n_d  += n_c[i] * n_d[i];
                    n_a_dot_n_c  += n_a[i] * n_c[i];
                    n_d_dot_n_ab += n_d[i] * n_ab[i];
                    n_b_dot_n_d  += n_b[i] * n_d[i];
                }

                // Compute the integrand and multiply by jaconian determinant
                double complex F = (1.0 / (r_c*r_c * r_d*r_d) *
                    ( (n_c_dot_n_ab - n_a_dot_n_c) * (n_d_dot_n_ab + n_b_dot_n_d) / (s_ab*s_ab)
                    - (n_c_dot_n_d - n_c_dot_n_ab * n_d_dot_n_ab) / (r_ab * s_ab)
                    ));
                F *= jacobian;

                // Undo complex step in the given direction
                pos_cmplx[coord_index] -= I * h;

                // Store derivative
                int outIndex = NUM_LOCAL_BODIES * 3 * p + 3 * body + axis;
                ff[outIndex] = (cubareal)(cimag(F) / h);
            }
        }
    }
    return 0;
}


double ln_integral_sum(double* w, struct ode_params* params) {
    IntegralParams integral_params;
    int num_bodies = params->num_bodies;
    int num_dim = params->num_dim;
    double integral_sum_value = 0.0;
    cubareal integral[6], error[6], prob[6];
    int neval, fail, nregions;

    if (num_dim != 3)
        errorexit("The 2PN ln-integral can only be computed in 3D! Please use num_dim = 3");

    // Loop over unordered quadruples a<b<c<d
    for (int a = 0; a < num_bodies; ++a) {
        for (int b = a+1; b < num_bodies; ++b) {
            for (int c = b+1; c < num_bodies; ++c) {
                for (int d = c+1; d < num_bodies; ++d) {

                    // Prepare local positions for this quadruple
                    for (int j = 0; j < 3; ++j) {
                        integral_params.pos[0][j] = w[3 * a + j];
                        integral_params.pos[1][j] = w[3 * b + j];
                        integral_params.pos[2][j] = w[3 * c + j];
                        integral_params.pos[3][j] = w[3 * d + j];
                    }

                    // Integrate 6 symmetry-inequivalent terms
                    Cuhre(num_dim, 6, ln_integral, &integral_params, NVEC,
                          EPSREL, EPSABS, 0,
                          MINEVAL, MAXEVAL, KEY,
                          NULL, NULL,
                          &nregions, &neval, &fail,
                          integral, error, prob);

                    // Accumulate contributions to full quadruple sum
                    // Each value should be multiplied by the symmetry factor 4, but in the Hamiltonian there is a factor 1/(4*PI), so they cancel
                    for (int i = 0; i < 6; ++i)
                        integral_sum_value += params->masses[a] * params->masses[b] * params->masses[c] * params->masses[d] * integral[i];
                }
            }
        }
    }
    return integral_sum_value;
}


complex double UTT4_without_ln_integral_complex(complex double* w, struct ode_params* params) {
    int num_bodies = params->num_bodies;
    int num_dim = params->num_dim;
    int array_half = num_bodies * num_dim; 
    complex double temp0, temp1, temp2, temp3;
    double ma_inv, mb_inv, mc_inv;

    // Masses
    double m[num_bodies];
    for (int a = 0; a < num_bodies; a++)
        m[a] = params->masses[a];
    
    // Momenta
    complex double p[num_bodies][num_dim];
    for (int a = 0; a < num_bodies; a++)
        for (int i = 0; i < num_dim; i++)
            p[a][i] = w[array_half + a * num_dim + i];
    
    // Relative positions and distances
    complex double x_rel[num_bodies][num_bodies][num_dim]; 
    complex double n[num_bodies][num_bodies][num_dim];
    complex double r[num_bodies][num_bodies];
    complex double n_ab_ac[num_dim];
    complex double n_ab_cb[num_dim];
    for (int a = 0; a < num_bodies; a++) {
        for (int b = a; b < num_bodies; b++) {
            for (int i = 0; i < num_dim; i++){
                x_rel[a][b][i] = w[a * num_dim + i] - w[b * num_dim + i];
                x_rel[b][a][i] = -x_rel[a][b][i];
            } 
            r[a][b] = norm_c(x_rel[a][b], num_dim);
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

    //Compute UTT4 without the ln-integral
    complex double UTT4 = 0.0;

    for (int a = 0; a < num_bodies; ++a) {
        for (int b = 0; b < num_bodies; ++b) {
            for (int c = 0; c < num_bodies; ++c) {
                for (int d = 0; d < num_bodies; ++d) {
                    if (b != a && c != a && c != b && d != a && d != b && d != c) {
                        temp0 = r[a][b] * r[a][b];
                        temp1 = r[b][c] * r[b][c];
                        temp2 = r[c][d] * r[c][d];
                        temp3 = r[a][d] * r[a][d];
                        UTT4 += - 0.0625 * m[a] * m[b] * m[c] * m[d] / (temp0*r[a][b] * temp2*r[c][d] * temp3*r[a][d] * temp1*r[b][c]) * (
                                4 * temp0*r[a][b] * temp1*r[b][c] * temp2 * temp3 / r[b][d]
                                -6* temp1*r[b][c] * temp3 * (temp0 * temp2 + temp3 * (temp3 + temp1 - r[a][c]*r[a][c] - r[b][d]*r[b][d]))
                                - temp0 * (r[b][d]*r[b][d] - temp1 - temp2) * (2 * r[a][b] * temp3*r[a][d] * temp1 / (r[a][c] + r[b][c] + r[a][b]) 
                                - temp3 * temp1 - r[a][b] * temp2 * (r[a][c]*r[a][c] - temp3 - temp2)) );
                    }
                }
            }
        }
    }
    return UTT4;
}


void compute_dUTT4_dx(double*w, struct ode_params* params, double *dUdx) {
    IntegralParams integral_params;
    int num_bodies = params->num_bodies;
    int num_dim = params->num_dim;
    int array_half = num_dim * num_bodies;
    complex double w_c[array_half];
    complex double UTT4_cs_val;

    const double h = 1e-20;

    if (num_dim != 3)
        errorexit("The 2PN ln-integral can only be computed in 3D! Please use num_dim = 3");

    for (int i = 0; i < array_half; ++i)
        dUdx[i] = 0.0;

    // UTT4 part without the ln-integral
    // Copy position part of the original array to w_c
    for (int i = 0; i < array_half; ++i)
        w_c[i] = (complex double)w[i];

    for (int i = 0; i < array_half; ++i) {
        // Add tiny imaginary step in coordinate i
        w_c[i] += I * h; 

        // Compute derivative
        UTT4_cs_val = UTT4_without_ln_integral_complex(w_c, params);
        dUdx[i] = cimag(UTT4_cs_val) / h;

        // Restore original value
        w_c[i] = (complex double)w[i];
    }


    // Add the ln-integral part
    cubareal integral[NCOMP_H_DERIV], error[NCOMP_H_DERIV], prob[NCOMP_H_DERIV];
    int neval, fail, nregions;

    // Loop over unordered quadruples a<b<c<d
    for (int a = 0; a < num_bodies; ++a) {
        for (int b = a+1; b < num_bodies; ++b) {
            for (int c = b+1; c < num_bodies; ++c) {
                for (int d = c+1; d < num_bodies; ++d) {

                    double mass_factor = params->masses[a] * params->masses[b] * params->masses[c] * params->masses[d];

                    // Local -> global index mapping
                    int slot2global[4] = { a, b, c, d };

                    // Prepare local positions for this quadruple
                    for (int j = 0; j < 3; ++j) {
                        integral_params.pos[0][j] = w[3 * a + j];
                        integral_params.pos[1][j] = w[3 * b + j];
                        integral_params.pos[2][j] = w[3 * c + j];
                        integral_params.pos[3][j] = w[3 * d + j];
                    }

                    // Integrate: one call, 72 outputs
                    Cuhre(num_dim, NCOMP_H_DERIV, ln_integral_complex, &integral_params, NVEC,
                          EPSREL, EPSABS, 0,
                          MINEVAL, MAXEVAL, KEY,
                          NULL, NULL,
                          &nregions, &neval, &fail,
                          integral, error, prob);

                    // Accumulate derivatives into global gradient
                    for (int p = 0; p < 6; ++p) {
                        for (int body = 0; body < 4; ++body) {
                            int global_idx = slot2global[body];
                            for (int axis = 0; axis < 3; ++axis) {
                                int local_comp = 12*p + 3*body + axis;
                                // Each value should be multiplied by the symmetry factor 4, but in the Hamiltonian there is a factor 1/(4*PI), so they cancel
                                dUdx[3*global_idx + axis] += INVPI * mass_factor * integral[local_comp];
                            }
                        }
                    }
                }
            }
        }
    }
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
    double m[num_bodies];
    for (int a = 0; a < num_bodies; a++)
        m[a] = params->masses[a];
    
    // Momenta
    double p[num_bodies][num_dim];
    for (int a = 0; a < num_bodies; a++)
        for (int i = 0; i < num_dim; i++)
            p[a][i] = w[array_half + a * num_dim + i];
    
    // Relative positions and distances
    double x_rel[num_bodies][num_bodies][num_dim]; 
    double n[num_bodies][num_bodies][num_dim];
    double r[num_bodies][num_bodies];
    double n_ab_ac[num_dim];
    double n_ab_cb[num_dim];
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
    return H;
}


double H2PN_nbody(double* w, struct ode_params* params, int p_flag) {
    int num_bodies = params->num_bodies;
    int num_dim = params->num_dim;
    int array_half = num_bodies * num_dim; 
    double temp, temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13;
    double ma_inv, mb_inv, mc_inv;

    // Masses
    double m[num_bodies];
    for (int a = 0; a < num_bodies; a++)
        m[a] = params->masses[a];
    
    // Momenta
    double p[num_bodies][num_dim];
    for (int a = 0; a < num_bodies; a++)
        for (int i = 0; i < num_dim; i++)
            p[a][i] = w[array_half + a * num_dim + i];
    
    // Relative positions and distances
    double x_rel[num_bodies][num_bodies][num_dim]; 
    double n[num_bodies][num_bodies][num_dim];
    double r[num_bodies][num_bodies];
    double n_ab_ac[num_dim];
    double n_ab_cb[num_dim];
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
    for (int a = 0; a < num_bodies; a++) {
        ma_inv = 1/m[a];
        temp = dot_product(p[a], p[a], num_dim);
        temp2 = temp * ma_inv * ma_inv;

        H += 0.0625 * m[a] * temp2 * temp2 * temp2;

        for (int b = 0; b < num_bodies; b++) {
            mb_inv = 1/m[b];
            temp0 = r[a][b] * r[a][b];
            temp1 = m[a] * m[b];
            temp3 = temp * ma_inv * mb_inv;
            temp4 = dot_product(n[a][b], p[b], num_dim);
            temp5 = dot_product(n[a][b], p[a], num_dim);
            temp6 = dot_product(p[b], p[b], num_dim);
            temp7 = dot_product(p[a], p[b], num_dim) * ma_inv * mb_inv;

            if (b != a){
                H += 0.0625 * 1 / r[a][b] * (10 * temp1 * temp2 * temp2
                    - 11 * temp3 * temp6
                    - 2 * temp3 * temp
                    + 10 * temp3 * temp4 * temp4
                    - 12 * temp7 * temp5 * temp4
                    - 3 * temp5 * temp5 * temp4 * temp4 * ma_inv * mb_inv);
                H += 0.25 * m[a] * temp1 / temp0 * (temp2
                    + temp6 * mb_inv * mb_inv
                    - 2 * temp7);
                if (!p_flag)
                    H += -0.25 * temp1 * temp1 / (temp0 * r[a][b]);
            }
            for (int c = 0; c < num_bodies; c++) {
                mc_inv = 1/m[c];
                temp8 = dot_product(n[a][c], p[c], num_dim);
                temp9 = dot_product(n[a][b], n[a][c], num_dim);
                temp10 = dot_product(n[a][b], p[c], num_dim);
                temp11 = r[b][c] * r[b][c];

                if (b != a && c != a) {
                    H += 0.125 * temp1 * m[c] / (r[a][b] * r[a][c]) * (18 * temp2 
                        + 14 * temp6 * mb_inv * mb_inv
                        - 2 * temp4 * temp4 * mb_inv * mb_inv
                        - 50 * temp7
                        + 17 * dot_product(p[b], p[c], num_dim) * mb_inv * mc_inv
                        - 14 * temp5 * temp4 * ma_inv * mb_inv
                        + 14 * temp4 * temp10 * mb_inv * mc_inv
                        + temp9 * temp4 * temp8 * mb_inv * mc_inv);
                    H += 0.125 * temp1 * m[c] / (r[a][b] * r[a][b]) * (2 * temp5 * temp8 * ma_inv * mc_inv
                        + 2 * temp4 * temp8 * ma_inv * mc_inv
                        + 5 * temp9 * dot_product(p[c], p[c], num_dim) * mc_inv * mc_inv
                        - temp9 * temp8 * temp8 * mc_inv * mc_inv
                        - 14 * temp10 * temp8 * mc_inv * mc_inv);
                }
                if (b != a && c != a && c != b) {
                    for (int i = 0; i < num_dim; i++) {
                        n_ab_ac[i] = n[a][b][i] + n[a][c][i];
                        n_ab_cb[i] = n[a][b][i] + n[c][b][i];
                    }
                    H += 0.5 * temp1 * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * (r[a][b] + r[b][c] + r[c][a])) * (
                        8 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[c], num_dim) * ma_inv * mc_inv
                        - 16 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[a], num_dim) * ma_inv * mc_inv
                        + 3 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[b], num_dim) * ma_inv * mb_inv
                        + 4 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[c], num_dim) * mc_inv * mc_inv
                        + dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[a], num_dim) * ma_inv * ma_inv);
                    H += 0.5 * temp1 * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * r[a][b]) * (
                        8 * (dot_product(p[a], p[c], num_dim) - temp5 * temp10) * ma_inv * mc_inv
                        - 3 * (temp7 - temp5 * temp4 * ma_inv * mb_inv) 
                        - 4 * (dot_product(p[c], p[c], num_dim) - temp10 * temp10) * mc_inv * mc_inv
                        - (temp - temp5 * temp5) * ma_inv * ma_inv);
                    if (!p_flag)
                        H += -0.015625 * m[a] * temp1 * m[c] / (temp0 * r[a][b] * r[a][c] * r[a][c] * r[a][c] * r[b][c]) * (
                            18 * temp0 * r[a][c] * r[a][c] - 60 * temp0 * temp11
                            - 24 * temp0 * r[a][c] * (r[a][b] + r[b][c])
                            + 60 * r[a][b] * r[a][c] * temp11 + 56 * temp0 * r[a][b] * r[b][c]
                            - 72 * r[a][b] * r[b][c] * temp11 + 35 * temp11 * temp11 + 6 * temp0 * temp0);
                }
                // G^3 quadruple sum terms
                if (!p_flag) {
                    for (int d = 0; d < num_bodies; d++) {
                        temp12 = r[c][d]*r[c][d];
                        temp13 = r[a][d]*r[a][d];

                        if (b != a && c != b && d != c) {
                            H += - 0.375 * temp1 * m[c] * m[d] / (r[a][b] * r[b][c] * r[c][d]);
                        }
                        if (b != a && c != a && d != a) {
                            H += - 0.25 * temp1 * m[c] * m[d] / (r[a][b] * r[a][c] * r[a][d]);
                        }
                        if (b != a && c != a && c != b && d != a && d != b && d != c) {
                            H += - 0.0625 * temp1 * m[c] * m[d] / (temp0*r[a][b] * temp12*r[c][d] * temp13*r[a][d] * temp11*r[b][c]) * (
                                4 * temp0*r[a][b] * temp11*r[b][c] * temp12 * temp13 / r[b][d]
                                -6* temp11*r[b][c] * temp13 * (temp0 * temp12 + temp13 * (temp13 + temp11 - r[a][c]*r[a][c] - r[b][d]*r[b][d]))
                                - temp0 * (r[b][d]*r[b][d] - temp11 - temp12) * (2 * r[a][b] * temp13*r[a][d] * temp11 / (r[a][c] + r[b][c] + r[a][b]) 
                                - temp13 * temp11 - r[a][b] * temp12 * (r[a][c]*r[a][c] - temp13 - temp12)) );
                        }
                    }
                }
            }
        } 
    }
    // ln-integral sum
    if (!p_flag)
        H += INVPI * ln_integral_sum(w, params);
    return H;
}


complex double H2PN_nbody_base_complex(complex double* w, struct ode_params* params, int p_flag) {
    int num_bodies = params->num_bodies;
    int num_dim = params->num_dim;
    int array_half = num_bodies * num_dim; 
    complex double temp, temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13;
    double ma_inv, mb_inv, mc_inv;

    // Masses
    double m[num_bodies];
    for (int a = 0; a < num_bodies; a++)
        m[a] = params->masses[a];
    
    // Momenta
    complex double p[num_bodies][num_dim];
    for (int a = 0; a < num_bodies; a++)
        for (int i = 0; i < num_dim; i++)
            p[a][i] = w[array_half + a * num_dim + i];
    
    // Relative positions and distances
    complex double x_rel[num_bodies][num_bodies][num_dim]; 
    complex double n[num_bodies][num_bodies][num_dim];
    complex double r[num_bodies][num_bodies];
    complex double n_ab_ac[num_dim];
    complex double n_ab_cb[num_dim];
    for (int a = 0; a < num_bodies; a++) {
        for (int b = a; b < num_bodies; b++) {
            for (int i = 0; i < num_dim; i++){
                x_rel[a][b][i] = w[a * num_dim + i] - w[b * num_dim + i];
                x_rel[b][a][i] = -x_rel[a][b][i];
            } 
            r[a][b] = norm_c(x_rel[a][b], num_dim);
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
    complex double H = 0.0;
    for (int a = 0; a < num_bodies; a++) {
        ma_inv = 1/m[a];
        temp = dot_product_c(p[a], p[a], num_dim);
        temp2 = temp * ma_inv * ma_inv;

        H += 0.0625 * m[a] * temp2 * temp2 * temp2;

        for (int b = 0; b < num_bodies; b++) {
            mb_inv = 1/m[b];
            temp0 = r[a][b] * r[a][b];
            temp1 = m[a] * m[b];
            temp3 = temp * ma_inv * mb_inv;
            temp4 = dot_product_c(n[a][b], p[b], num_dim);
            temp5 = dot_product_c(n[a][b], p[a], num_dim);
            temp6 = dot_product_c(p[b], p[b], num_dim);
            temp7 = dot_product_c(p[a], p[b], num_dim) * ma_inv * mb_inv;

            if (b != a){
                H += 0.0625 * 1 / r[a][b] * (10 * temp1 * temp2 * temp2
                    - 11 * temp3 * temp6
                    - 2 * temp3 * temp
                    + 10 * temp3 * temp4 * temp4
                    - 12 * temp7 * temp5 * temp4
                    - 3 * temp5 * temp5 * temp4 * temp4 * ma_inv * mb_inv);
                H += 0.25 * m[a] * temp1 / temp0 * (temp2
                    + temp6 * mb_inv * mb_inv
                    - 2 * temp7);
                if (!p_flag)
                    H += -0.25 * temp1 * temp1 / (temp0 * r[a][b]);
            }
            for (int c = 0; c < num_bodies; c++) {
                mc_inv = 1/m[c];
                temp8 = dot_product_c(n[a][c], p[c], num_dim);
                temp9 = dot_product_c(n[a][b], n[a][c], num_dim);
                temp10 = dot_product_c(n[a][b], p[c], num_dim);
                temp11 = r[b][c] * r[b][c];

                if (b != a && c != a) {
                    H += 0.125 * temp1 * m[c] / (r[a][b] * r[a][c]) * (18 * temp2 
                        + 14 * temp6 * mb_inv * mb_inv
                        - 2 * temp4 * temp4 * mb_inv * mb_inv
                        - 50 * temp7
                        + 17 * dot_product_c(p[b], p[c], num_dim) * mb_inv * mc_inv
                        - 14 * temp5 * temp4 * ma_inv * mb_inv
                        + 14 * temp4 * temp10 * mb_inv * mc_inv
                        + temp9 * temp4 * temp8 * mb_inv * mc_inv);
                    H += 0.125 * temp1 * m[c] / (r[a][b] * r[a][b]) * (2 * temp5 * temp8 * ma_inv * mc_inv
                        + 2 * temp4 * temp8 * ma_inv * mc_inv
                        + 5 * temp9 * dot_product_c(p[c], p[c], num_dim) * mc_inv * mc_inv
                        - temp9 * temp8 * temp8 * mc_inv * mc_inv
                        - 14 * temp10 * temp8 * mc_inv * mc_inv);
                }
                if (b != a && c != a && c != b) {
                    for (int i = 0; i < num_dim; i++) {
                        n_ab_ac[i] = n[a][b][i] + n[a][c][i];
                        n_ab_cb[i] = n[a][b][i] + n[c][b][i];
                    }
                    H += 0.5 * temp1 * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * (r[a][b] + r[b][c] + r[c][a])) * (
                        8 * dot_product_c(n_ab_ac, p[a], num_dim) * dot_product_c(n_ab_cb, p[c], num_dim) * ma_inv * mc_inv
                        - 16 * dot_product_c(n_ab_ac, p[c], num_dim) * dot_product_c(n_ab_cb, p[a], num_dim) * ma_inv * mc_inv
                        + 3 * dot_product_c(n_ab_ac, p[a], num_dim) * dot_product_c(n_ab_cb, p[b], num_dim) * ma_inv * mb_inv
                        + 4 * dot_product_c(n_ab_ac, p[c], num_dim) * dot_product_c(n_ab_cb, p[c], num_dim) * mc_inv * mc_inv
                        + dot_product_c(n_ab_ac, p[a], num_dim) * dot_product_c(n_ab_cb, p[a], num_dim) * ma_inv * ma_inv);
                    H += 0.5 * temp1 * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * r[a][b]) * (
                        8 * (dot_product_c(p[a], p[c], num_dim) - temp5 * temp10) * ma_inv * mc_inv
                        - 3 * (temp7 - temp5 * temp4 * ma_inv * mb_inv) 
                        - 4 * (dot_product_c(p[c], p[c], num_dim) - temp10 * temp10) * mc_inv * mc_inv
                        - (temp - temp5 * temp5) * ma_inv * ma_inv);
                    if (!p_flag)
                        H += -0.015625 * m[a] * temp1 * m[c] / (temp0 * r[a][b] * r[a][c] * r[a][c] * r[a][c] * r[b][c]) * (
                            18 * temp0 * r[a][c] * r[a][c] - 60 * temp0 * temp11
                            - 24 * temp0 * r[a][c] * (r[a][b] + r[b][c])
                            + 60 * r[a][b] * r[a][c] * temp11 + 56 * temp0 * r[a][b] * r[b][c]
                            - 72 * r[a][b] * r[b][c] * temp11 + 35 * temp11 * temp11 + 6 * temp0 * temp0);
                }
                // G^3 quadruple sum terms
                if (!p_flag) {
                    for (int d = 0; d < num_bodies; d++) {
                        temp12 = r[c][d]*r[c][d];
                        temp13 = r[a][d]*r[a][d];

                        if (b != a && c != b && d != c) {
                            H += - 0.375 * temp1 * m[c] * m[d] / (r[a][b] * r[b][c] * r[c][d]);
                        }
                        if (b != a && c != a && d != a) {
                            H += - 0.25 * temp1 * m[c] * m[d] / (r[a][b] * r[a][c] * r[a][d]);
                        }
                    }
                }
            }
        } 
    }
    return H;
}


void update_eom_hamiltonian_fd(double *w, double *dwdt, double (*hamiltonian)(double*, struct ode_params*, int), double h, struct ode_params* params) {
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


void update_eom_hamiltonian_cs(double *w, double *dwdt, complex double (*hamiltonian)(complex double*, struct ode_params*, int), double h, struct ode_params* params) {
    int array_half = params->num_dim * params->num_bodies;
    int total_dim = 2 * array_half;
    complex double w_c[total_dim];
    complex double H_cs_val;
    double dHdw[total_dim];
    int p_flag;

    // Copy original array to w_copy
    for (int i = 0; i < total_dim; ++i)
        w_c[i] = (complex double)w[i];

    for (int i = 0; i < total_dim; ++i) {
        p_flag = (i < array_half) ? 0 : 1;

        // Add tiny imaginary step in coordinate i
        w_c[i] += I * h; 

        // Compute derivative
        H_cs_val = hamiltonian(w_c, params, p_flag);
        dHdw[i] = cimag(H_cs_val) / h;

        // Restore original value
        w_c[i] = (complex double)w[i];
    }

    // Compute dwdt
    for (int i = 0; i < 2 * array_half; i++) {
        if (i < array_half)
            dwdt[i] += dHdw[i + array_half];
        else
            dwdt[i] += -dHdw[i - array_half];
    }
}