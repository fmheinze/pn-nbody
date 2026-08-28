/**
 * @file eom_dissipative.c
 * @brief Dissipative 2.5PN equations-of-motion building blocks.
 */

#include "cache.h"
#include "eom_internal.h"
#include "utils.h"


// ------------------------------------------------------------------------------------------------
// 2.5PN analytic building blocks
// ------------------------------------------------------------------------------------------------

// Compute the tensor chi_dot entering the exact 2.5PN contribution.
void pn_compute_25pn_chi_dot(const PairCache *cache, int num_bodies, int num_dim,
    double x_dot[num_bodies][num_dim], double p_dot[num_bodies][num_dim],
    double C[num_dim][num_dim])
{
    const ActiveList *active = &cache->active;

    for (int i = 0; i < num_dim; i++)
        for (int j = 0; j < num_dim; j++)
            C[i][j] = 0.0;

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];

        double pdot_dot_p = 0.0;
        for (int k = 0; k < num_dim; k++)
            pdot_dot_p += p_dot[a][k] * pair_cache_p(cache, a, k);

        const double pref = 2.0 * cache->inv_m[a];

        for (int i = 0; i < num_dim; i++) {
            const double pdot_ai = p_dot[a][i];
            const double p_ai = pair_cache_p(cache, a, i);

            for (int j = 0; j < num_dim; j++) {
                C[i][j] += pref * (
                    2.0 * pdot_dot_p * delta(i, j)
                  - 3.0 * (pdot_ai * pair_cache_p(cache, a, j) + p_ai * p_dot[a][j])
                );
            }
        }
    }

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];

        for (int ib = ia + 1; ib < active->num_active; ib++) {
            const int b = active->ids[ib];

            double v_ab[num_dim];
            for (int k = 0; k < num_dim; k++)
                v_ab[k] = x_dot[a][k] - x_dot[b][k];

            double n_dot_v = 0.0;
            for (int k = 0; k < num_dim; k++)
                n_dot_v += pair_cache_n(cache, a, b, k) * v_ab[k];

            const double inv_r = pair_cache_inv_r(cache, a, b);
            const double pref = 2.0 * cache->m[a] * cache->m[b] * inv_r * inv_r;

            for (int i = 0; i < num_dim; i++) {
                const double ni = pair_cache_n(cache, a, b, i);

                for (int j = 0; j < num_dim; j++) {
                    const double nj = pair_cache_n(cache, a, b, j);

                    C[i][j] += pref * (
                        3.0 * (v_ab[i] * nj + ni * v_ab[j])
                      + n_dot_v * (delta(i, j) - 9.0 * ni * nj)
                    );
                }
            }
        }
    }
}


// Compute the 2.5PN velocity contribution to dw/dt.
void pn_add_25pn_velocity_contribution(const PairCache *cache, int num_dim,
    double C[num_dim][num_dim], double *velocity)
{
    const ActiveList *active = &cache->active;

    double trC = 0.0;
    for (int i = 0; i < num_dim; i++)
        trC += C[i][i];

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];

        for (int k = 0; k < num_dim; k++) {
            double Cp_k = 0.0;

            for (int j = 0; j < num_dim; j++)
                Cp_k += C[k][j] * pair_cache_p(cache, a, j);

            // Direct contraction: (1/45) C_ij d chi_ij / d p_a^k
            velocity[a * num_dim + k] += (4.0/45.0) * cache->inv_m[a]
              * pair_cache_p(cache, a, k) * trC
              - (12.0 / 45.0) * cache->inv_m[a] * Cp_k;
        }
    }
}


// Compute the 2.5PN momentum contribution to dw/dt.
void pn_add_25pn_momentum_contribution(const PairCache *cache, int num_dim,
    double C[num_dim][num_dim], double *momentum)
{
    const ActiveList *active = &cache->active;

    double trC = 0.0;
    for (int i = 0; i < num_dim; i++)
        trC += C[i][i];

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];

        for (int ib = ia + 1; ib < active->num_active; ib++) {
            const int b = active->ids[ib];

            const double inv_r = pair_cache_inv_r(cache, a, b);
            const double pref = 2.0 * cache->m[a] * cache->m[b] * inv_r * inv_r / 45.0;

            double Cn[num_dim];

            for (int k = 0; k < num_dim; k++) {
                Cn[k] = 0.0;

                for (int j = 0; j < num_dim; j++)
                    Cn[k] += C[k][j] * pair_cache_n(cache, a, b, j);
            }

            double nCn = 0.0;
            for (int k = 0; k < num_dim; k++)
                nCn += pair_cache_n(cache, a, b, k) * Cn[k];

            for (int k = 0; k < num_dim; k++) {
                const double nk = pair_cache_n(cache, a, b, k);

                // For symmetric C: G_k = 6(Cn)_k + (trC - 9 n.C.n) n_k
                const double Gk = 6.0 * Cn[k] + (trC - 9.0 * nCn) * nk;

                momentum[a * num_dim + k] -= pref * Gk;
                momentum[b * num_dim + k] += pref * Gk;
            }
        }
    }
}
