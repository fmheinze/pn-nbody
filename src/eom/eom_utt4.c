/**
 * @file eom_utt4.c
 * @brief Analytic Cartesian gradient of the four-body UTT4 potential.
 */

#include <math.h>

#include "cache.h"
#include "eom.h"
#include "hamiltonian.h"
#include "utils.h"

#define UTT4_NUM_LOCAL_BODIES 4
#define UTT4_NUM_LOCAL_PAIRS 6

// ------------------------------------------------------------------------------------------------
// UTT4 building blocks
// ------------------------------------------------------------------------------------------------

// Local pair ordering: 01, 02, 03, 12, 13, 23
static const int utt4_pair_i[UTT4_NUM_LOCAL_PAIRS] = {0, 0, 0, 1, 1, 2};
static const int utt4_pair_j[UTT4_NUM_LOCAL_PAIRS] = {1, 2, 3, 2, 3, 3};


static int utt4_local_pair_index(int a, int b)
{
    if (a > b) {
        const int tmp = a;
        a = b;
        b = tmp;
    }

    if (a == 0 && b == 1) return 0;
    if (a == 0 && b == 2) return 1;
    if (a == 0 && b == 3) return 2;
    if (a == 1 && b == 2) return 3;
    if (a == 1 && b == 3) return 4;
    if (a == 2 && b == 3) return 5;
    return -1;
}


/**
 * @brief Add the analytic distance derivatives of one ordered non-logarithmic UTT4 term.
 *
 * @param r             Six local pair distances in ordering 01,02,03,12,13,23.
 * @param a,b,c,d       Ordered local body labels (all distinct, each in [0,3]).
 * @param dF_dr         Accumulator for derivatives with respect to the six local pair distances.
 */
static void utt4_nonlog_ordered_distance_gradient(const double r[UTT4_NUM_LOCAL_PAIRS],
    int a, int b, int c, int d, double dF_dr[UTT4_NUM_LOCAL_PAIRS])
{
    const int pab = utt4_local_pair_index(a, b);
    const int pac = utt4_local_pair_index(a, c);
    const int pad = utt4_local_pair_index(a, d);
    const int pbc = utt4_local_pair_index(b, c);
    const int pbd = utt4_local_pair_index(b, d);
    const int pcd = utt4_local_pair_index(c, d);

    const double x = r[pab];
    const double y = r[pac];
    const double z = r[pad];
    const double u = r[pbc];
    const double v = r[pbd];
    const double w = r[pcd];

    const double x2 = x*x, x3 = x2*x;
    const double y2 = y*y;
    const double z2 = z*z, z3 = z2*z, z4 = z2*z2;
    const double u2 = u*u, u3 = u2*u;
    const double v2 = v*v;
    const double w2 = w*w, w3 = w2*w;

    const double S = x + y + u;
    const double invS = 1.0 / S;
    const double invS2 = invS * invS;

    const double R = z2 + u2 - y2 - v2;
    const double P = v2 - u2 - w2;
    const double Q = y2 - z2 - w2;

    const double C = -8.0 + 16.0*x*invS;
    const double C_x =  16.0*(y + u)*invS2;
    const double C_y = -16.0*x*invS2;
    const double C_u = C_y;

    const double K = z3*u2*C + x*w2*Q;
    const double K_x = z3*u2*C_x + w2*Q;
    const double K_y = z3*u2*C_y + 2.0*x*w2*y;
    const double K_z = 3.0*z2*u2*C - 2.0*x*w2*z;
    const double K_u = 2.0*z3*u*C + z3*u2*C_u;
    const double K_v = 0.0;
    const double K_w = 2.0*x*w*(Q - w2);

    const double T1 =  16.0*x3*u3*w2*z2/v;
    const double T2 = -24.0*u3*x2*w2*z2;
    const double T3 = -30.0*z4*u3*R;
    const double T4 = x2*P*K;
    const double B = T1 + T2 + T3 + T4;

    // Derivatives of B with respect to x, y, z, u, v, w
    const double B_x = 3.0*T1/x + 2.0*T2/x + 2.0*x*P*K + x2*P*K_x;
    const double B_y = 60.0*y*z4*u3 + x2*P*K_y;
    const double B_z = 2.0*T1/z + 2.0*T2/z - 30.0*z3*u3*(4.0*R + 2.0*z2) + x2*P*K_z;
    const double B_u = 3.0*T1/u + 3.0*T2/u - 30.0*z4*u2*(3.0*R + 2.0*u2) + x2*(-2.0*u*K + P*K_u);
    const double B_v = -T1/v + 60.0*v*z4*u3 + 2.0*x2*v*K + x2*P*K_v;
    const double B_w = 2.0*T1/w + 2.0*T2/w + x2*(-2.0*w*K + P*K_w);

    const double D = x3*z3*u3*w3;
    const double prefactor = -1.0 / (64.0*D);

    // Derivatives with respect to x, y, z, u, v, w
    const double F_x = prefactor * (B_x - 3.0*B/x);
    const double F_y = prefactor * B_y;
    const double F_z = prefactor * (B_z - 3.0*B/z);
    const double F_u = prefactor * (B_u - 3.0*B/u);
    const double F_v = prefactor * B_v;
    const double F_w = prefactor * (B_w - 3.0*B/w);

    dF_dr[pab] += F_x;
    dF_dr[pac] += F_y;
    dF_dr[pad] += F_z;
    dF_dr[pbc] += F_u;
    dF_dr[pbd] += F_v;
    dF_dr[pcd] += F_w;
}


/**
 * @brief Gradient of the complete 24-permutation non-logarithmic UTT4 contribution.
 *
 * The common four-mass product is not included. The result is reconstructed from the exact
 * distance derivatives through grad_{x_i} F = sum_{j != i} (dF/dr_ij) (x_i-x_j)/r_ij.
 */
static void utt4_nonlog_quadruple_gradient(const double pos[UTT4_NUM_LOCAL_BODIES][3],
    double grad[UTT4_NUM_LOCAL_BODIES][3])
{
    double r[UTT4_NUM_LOCAL_PAIRS];
    double n[UTT4_NUM_LOCAL_PAIRS][3];

    for (int p = 0; p < UTT4_NUM_LOCAL_PAIRS; ++p) {
        const int i = utt4_pair_i[p];
        const int j = utt4_pair_j[p];
        double r2 = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
            const double dx = pos[i][axis] - pos[j][axis];
            n[p][axis] = dx;
            r2 += dx*dx;
        }
        r[p] = sqrt(r2);
        for (int axis = 0; axis < 3; ++axis)
            n[p][axis] /= r[p];
    }

    double dF_dr[UTT4_NUM_LOCAL_PAIRS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int a = 0; a < UTT4_NUM_LOCAL_BODIES; ++a) {
        for (int b = 0; b < UTT4_NUM_LOCAL_BODIES; ++b) {
            if (b == a) continue;
            for (int c = 0; c < UTT4_NUM_LOCAL_BODIES; ++c) {
                if (c == a || c == b) continue;
                for (int d = 0; d < UTT4_NUM_LOCAL_BODIES; ++d) {
                    if (d == a || d == b || d == c) continue;
                    utt4_nonlog_ordered_distance_gradient(r, a, b, c, d, dF_dr);
                }
            }
        }
    }

    for (int local = 0; local < UTT4_NUM_LOCAL_BODIES; ++local)
        for (int axis = 0; axis < 3; ++axis)
            grad[local][axis] = 0.0;

    for (int p = 0; p < UTT4_NUM_LOCAL_PAIRS; ++p) {
        const int i = utt4_pair_i[p];
        const int j = utt4_pair_j[p];
        for (int axis = 0; axis < 3; ++axis) {
            const double contribution = dF_dr[p] * n[p][axis];
            grad[i][axis] += contribution;
            grad[j][axis] -= contribution;
        }
    }
}


/**
 * @brief Compute the Cartesian gradient of the complete UTT4 potential.
 *
 * The closed-form/non-logarithmic part is differentiated with explicit analytic derivatives
 * with respect to the six mutual distances. The logarithmic-integral part uses the system-level
 * cached gradient shared with the Hamiltonian energy evaluation.
 */
void compute_dUTT4_dx(double *w, struct ode_params *ode_params, double *dUdx)
{
    const int num_bodies = ode_params->num_bodies_initial;
    const int num_dim = ode_params->num_dim;
    const int array_half = num_dim * num_bodies;

    if (num_dim != 3)
        errorexit("UTT4 can only be computed in 3D! Please use num_dim = 3");

    PairCache *cache = pair_cache_get_workspace(ode_params);
    active_list_refresh(&cache->active, ode_params);
    const ActiveList *active = &cache->active;

    for (int i = 0; i < array_half; ++i)
        dUdx[i] = 0.0;

    // Closed-form/non-logarithmic contribution: exact analytic distance derivatives.
    for (int ia = 0; ia < active->num_active; ++ia) {
        const int a = active->ids[ia];
        for (int ib = ia + 1; ib < active->num_active; ++ib) {
            const int b = active->ids[ib];
            for (int ic = ib + 1; ic < active->num_active; ++ic) {
                const int c = active->ids[ic];
                for (int id = ic + 1; id < active->num_active; ++id) {
                    const int d = active->ids[id];
                    const int body[UTT4_NUM_LOCAL_BODIES] = {a, b, c, d};
                    double pos[UTT4_NUM_LOCAL_BODIES][3];

                    for (int local = 0; local < UTT4_NUM_LOCAL_BODIES; ++local)
                        for (int axis = 0; axis < 3; ++axis)
                            pos[local][axis] = w[num_dim * body[local] + axis];

                    const double mass_fac = ode_params->masses[a] * ode_params->masses[b]
                        * ode_params->masses[c] * ode_params->masses[d];

                    double nonlog_grad[UTT4_NUM_LOCAL_BODIES][3];
                    utt4_nonlog_quadruple_gradient(pos, nonlog_grad);
                    for (int local = 0; local < UTT4_NUM_LOCAL_BODIES; ++local) {
                        const int global = body[local];
                        for (int axis = 0; axis < 3; ++axis)
                            dUdx[num_dim * global + axis] += mass_fac * nonlog_grad[local][axis];
                    }
                }
            }
        }
    }

    // Logarithmic-integral contribution. The same cached value/gradient is shared with H2PN.
    double ln_grad[array_half];
    utt4_ln_integral_cached(w, ode_params, NULL, ln_grad);
    for (int i = 0; i < array_half; ++i)
        dUdx[i] += ln_grad[i];
}
