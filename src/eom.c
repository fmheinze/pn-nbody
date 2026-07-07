/**
 * @file eom.c
 * @brief Functions for the right-hand side of the post-Newtonian equations of motion
 */

#include <complex.h>
#include "utils.h"
#include "eom.h"
#include "hamiltonian.h"
#include "pair_cache.h"


/**
 * @brief Right-hand side of the N-body equations of motion up to 2.5PN order
 * 
 * @param[in]       t           Time (currently not used, but kept for completeness)
 * @param[in]       w           State of the system, w = [positions, momenta]
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[out]      dwdt        Right-hand side of the equations of motion
 */
void rhs_pn_nbody(double t, double* w, struct ode_params* ode_params, double* dwdt)
{
    (void)t;    // Unused

    // --------------------------------------------------------------------------------------------
    // Initialize the arrays
    // --------------------------------------------------------------------------------------------
    int num_bodies = ode_params->num_bodies_initial;
    int num_dim = ode_params->num_dim;
    int array_half = num_bodies * num_dim; 
    double result, temp;

    // Cached active-body, momentum and pair-geometry data.
    //
    // PairCache owns an ActiveList in cache.active.  All expensive RHS loops below now iterate over
    // cache.active.ids instead of scanning num_bodies_initial and skipping inactive bodies.
    PairCache cache;
    pair_cache_build(&cache, w, ode_params);
    const ActiveList *active = &cache.active;

    // Masses
    double m[num_bodies];
    double inv_m[num_bodies];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        m[a] = cache.m[a];
        inv_m[a] = cache.inv_m[a];
    }

    // Momenta
    double p[num_bodies][num_dim];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        for (int i = 0; i < num_dim; i++)
            p[a][i] = pair_cache_p(&cache, a, i);
    }

    // Relative positions and distances
    double x_rel[num_bodies][num_bodies][num_dim]; 
    double n[num_bodies][num_bodies][num_dim];
    double r[num_bodies][num_bodies];
    double inv_r[num_bodies][num_bodies];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        for (int ib = ia; ib < active->num_active; ib++) {
            int b = active->ids[ib];

            for (int i = 0; i < num_dim; i++) {
                x_rel[a][b][i] = pair_cache_xrel(&cache, a, b, i);
                x_rel[b][a][i] = pair_cache_xrel(&cache, b, a, i);
                n[a][b][i] = pair_cache_n(&cache, a, b, i);
                n[b][a][i] = pair_cache_n(&cache, b, a, i);
            }

            r[a][b] = pair_cache_r(&cache, a, b);
            r[b][a] = pair_cache_r(&cache, b, a);
            inv_r[a][b] = pair_cache_inv_r(&cache, a, b);
            inv_r[b][a] = pair_cache_inv_r(&cache, b, a);
        }
    }

    // Time derivatives
    double p_dot[num_bodies][num_dim];
    double x_dot[num_bodies][num_dim];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        for (int i = 0; i < num_dim; i++) {
            p_dot[a][i] = 0.0;
            x_dot[a][i] = 0.0;
        }
    }

    // Set ODE right-hand side initially to zero.
    // This also guarantees that inactive bodies have zero RHS, because all contribution loops below
    // iterate only over active->ids.
    for (int i = 0; i < 2 * array_half; i++)   
        dwdt[i] = 0.0;

    // --------------------------------------------------------------------------------------------
    // Add 0PN (Newtonian) terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[0]) {
        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];

            // Velocities
            for (int i = 0; i < num_dim; i++) {
                result = w[array_half + a * num_dim + i] * inv_m[a];
                dwdt[a * num_dim + i] += result;
                x_dot[a][i] += result;
            }

            // Accelerations.  Use ib = ia + 1 to keep the original pair-once symmetric update.
            for (int ib = ia + 1; ib < active->num_active; ib++) {
                int b = active->ids[ib];

                for (int i = 0; i < num_dim; i++) {
                    result = -m[a] * m[b] * inv_r[a][b] * inv_r[a][b] * n[a][b][i];
                    dwdt[array_half + a * num_dim + i] += result;
                    dwdt[array_half + b * num_dim + i] += -result;
                    p_dot[a][i] += result;
                    p_dot[b][i] += -result;
                }
            }
        }
    }

    // --------------------------------------------------------------------------------------------
    // Add 1PN terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[1]) {
        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];

            double pa_dot_pa = cache.p2[a];

            for (int i = 0; i < num_dim; i++) {

                // Velocities
                result = -0.5 * pa_dot_pa * inv_m[a] * inv_m[a] * inv_m[a] * p[a][i];
                x_dot[a][i] += result;
                dwdt[a * num_dim + i] += result;

                for (int ib = 0; ib < active->num_active; ib++) {
                    int b = active->ids[ib];

                    double pa_dot_pb = pair_cache_p_dot(&cache, a, b);
                    double pb_dot_pb = cache.p2[b];
                    double nab_dot_pa = pair_cache_n_dot_p(&cache, a, b, a);
                    double nab_dot_pb = pair_cache_n_dot_p(&cache, a, b, b);

                    if (b != a) {
                        temp = -0.5 * inv_r[a][b];
                        result = temp * (6 * m[b] * inv_m[a] * p[a][i] 
                            - 7 * p[b][i] - nab_dot_pb * n[a][b][i]);
                        x_dot[a][i] += result;
                        dwdt[a * num_dim + i] += result;

                        // Accelerations 
                        temp *= inv_r[a][b];
                        result = temp * (3 * m[b] * inv_m[a] * pa_dot_pa 
                            + 3 * m[a] * inv_m[b] * pb_dot_pb 
                            - 7 * pa_dot_pb - 3 * nab_dot_pa * nab_dot_pb) * n[a][b][i];
                        p_dot[a][i] += result;                       
                        dwdt[array_half + a * num_dim + i] += result;

                        result = temp * (nab_dot_pb * p[a][i] + nab_dot_pa * p[b][i]);
                        p_dot[a][i] += result;
                        dwdt[array_half + a * num_dim + i] += result;
                    }

                    for (int ic = 0; ic < active->num_active; ic++) {
                        int c = active->ids[ic];

                        temp = m[a] * m[b] * m[c] * inv_r[a][b] * inv_r[a][b] * n[a][b][i];

                        if (b != a && c != a) {
                            result = temp * inv_r[a][c];
                            p_dot[a][i] += result;
                            dwdt[array_half + a * num_dim + i] += result;
                        }

                        if (b != a && c != b) {
                            result = temp * inv_r[b][c];
                            p_dot[a][i] += result;
                            dwdt[array_half + a * num_dim + i] += result;
                        }
                    }
                }
            }
        }
    }

    // --------------------------------------------------------------------------------------------
    // Add 2PN terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[2]) {
        // Add the contributions from H2PN without UTT4
        update_eom_hamiltonian_cs(w, H2PN_base_complex, 1e-30, ode_params, dwdt);

        // If not using impulse splitting, add UTT4 contributions directly to dp/dt
        #if HAVE_CUBA
        if (ode_params->include_utt4 && !ode_params->use_impulse_method)
        {
            double dUdx[array_half];
            compute_dUTT4_dx(w, ode_params, dUdx);
            for (int ia = 0; ia < active->num_active; ia++) {
                int a = active->ids[ia];
                for (int i = 0; i < num_dim; i++) {
                    int idx = a * num_dim + i;
                    dwdt[array_half + idx] -= dUdx[idx];
                }
            }
        }
        #endif
    }

    // --------------------------------------------------------------------------------------------
    // Add 2.5PN terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[3]) {

        // Initialize arrays
        double x_rel_dot[num_bodies][num_bodies][num_dim];

        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];
            for (int ib = 0; ib < active->num_active; ib++) {
                int b = active->ids[ib];

                for (int i = 0; i < num_dim; i++)
                    x_rel_dot[a][b][i] = x_dot[a][i] - x_dot[b][i];
            }
        }

        double chi_dot[num_dim][num_dim];
        double dp_chi[num_bodies][num_dim][num_dim][num_dim];
        double dx_chi[num_bodies][num_dim][num_dim][num_dim];

        for (int i = 0; i < num_dim; i++)
            for (int j = 0; j < num_dim; j++)
                chi_dot[i][j] = 0.0;

        for (int ia = 0; ia < active->num_active; ia++) { 
            int a = active->ids[ia];
            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        dp_chi[a][i][j][k] = 0.0;
                        dx_chi[a][i][j][k] = 0.0;
                    }
                }
            }
        }

        // Compute Chi values
        for (int i = 0; i < num_dim; i++) {
            for (int j = 0; j < num_dim; j++) {
                for (int ia = 0; ia < active->num_active; ia++) {
                    int a = active->ids[ia];

                    chi_dot[i][j] += 2.0 * inv_m[a] * (
                        2 * dot_product(p_dot[a], p[a], num_dim) * delta(i, j) 
                        - 3 * (p_dot[a][i] * p[a][j] + p[a][i] * p_dot[a][j]) );

                    for (int ib = 0; ib < active->num_active; ib++) {
                        int b = active->ids[ib];

                        if (b != a) {
                            chi_dot[i][j] += m[a] * m[b] * inv_r[a][b] * inv_r[a][b] * (3 * 
                                (x_rel_dot[a][b][i] * n[a][b][j] + n[a][b][i] * x_rel_dot[a][b][j])
                                + dot_product(n[a][b], x_rel_dot[a][b], num_dim) * (delta(i, j) 
                                - 9 * n[a][b][i] * n[a][b][j]) );
                        }
                    }
                }
            }
        }

        for (int ic = 0; ic < active->num_active; ic++) {
            int c = active->ids[ic];

            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        dp_chi[c][i][j][k] += 2.0 / m[c] * (2 * p[c][k] * delta(i, j) 
                            - 3 * (p[c][j] * delta(i, k) + p[c][i] * delta(j, k)));

                        for (int ia = 0; ia < active->num_active; ia++) {
                            int a = active->ids[ia];

                            for (int ib = 0; ib < active->num_active; ib++) {
                                int b = active->ids[ib];

                                if (b != a) {
                                    dx_chi[c][i][j][k] += m[a] * m[b] * inv_r[a][b] * inv_r[a][b] *
                                        (delta(a, c) - delta(b, c)) * 
                                        (3 * (delta(i, k) * n[a][b][j] 
                                        + delta(j, k) * n[a][b][i]) 
                                        - 9 * n[a][b][k] * n[a][b][i] * n[a][b][j] 
                                        + delta(i, j) * n[a][b][k]);
                                }
                            }
                        }
                    }    
                }
            }
        }

        // Add contribution to the ODE right-hand side
        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];

            for (int k = 0; k < num_dim; k++) {
                for (int i = 0; i < num_dim; i++) {
                    for (int j = 0; j < num_dim; j++) {
                        dwdt[a * num_dim + k] += 1/45.0 * chi_dot[i][j] * dp_chi[a][i][j][k];
                        dwdt[array_half + a * num_dim + k] += -1/45.0 * chi_dot[i][j] 
                            * dx_chi[a][i][j][k];
                    }
                }
            }
        }
    }

    pair_cache_free(&cache);
}


/**
 * @brief Adds contribution from a Hamiltonian to the right-hand side of the equations of motion.
 * 
 * Adds contribution from a Hamiltonian to the right-hand side of the equations of motion,
 * according to dx/dt = dH/dp, dp/dt = -dH/dx. The derivatives of the Hamiltonian are computed
 * numerically using a complex-step derivative. The Hamiltonian must be of type complex double
 * with arguments (w, ode_params, p_flag), where p_flag just ignores all the terms that do not
 * have a momentum dependence for the computation of dH/dp.
 * 
 * @param[in]       w           State of the system, w = [positions, momenta]
 * @param[in]       H           Complex-valued Hamiltonian
 * @param[in]       h           Complex step size
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[out]      dwdt        Right-hand side of the ODE
 */
void update_eom_hamiltonian_cs(double *w, c_hamiltonian H, double h, struct ode_params* ode_params,
    double *dwdt)
{
    int num_dim = ode_params->num_dim;
    int num_bodies = ode_params->num_bodies_initial;
    int array_half = num_dim * num_bodies;
    int total_dim = 2 * array_half;
    complex double w_c[total_dim];
    complex double H_cs_val;
    double dHdw[total_dim];

    ActiveList active;
    active_list_init(&active, ode_params);

    // Copy original array to w_c and initialize dHdw
    for (int i = 0; i < total_dim; ++i) {
        w_c[i] = (complex double)w[i];
        dHdw[i] = 0.0;
    }

    // Position derivatives: dH/dx.  These are needed for dp/dt = -dH/dx.
    for (int ia = 0; ia < active.num_active; ia++) {
        int a = active.ids[ia];

        for (int k = 0; k < num_dim; k++) {
            int idx = a * num_dim + k;

            // Add tiny imaginary step in coordinate idx
            w_c[idx] += I * h;

            // p_flag = 0: keep position-only Hamiltonian terms
            H_cs_val = H(w_c, ode_params, 0);
            dHdw[idx] = cimag(H_cs_val) / h;

            // Restore original value
            w_c[idx] = (complex double)w[idx];
        }
    }

    // Momentum derivatives: dH/dp.  These are needed for dx/dt = dH/dp.
    for (int ia = 0; ia < active.num_active; ia++) {
        int a = active.ids[ia];

        for (int k = 0; k < num_dim; k++) {
            int idx = array_half + a * num_dim + k;

            // Add tiny imaginary step in momentum component idx
            w_c[idx] += I * h;

            // p_flag = 1: skip Hamiltonian terms without momentum dependence
            H_cs_val = H(w_c, ode_params, 1);
            dHdw[idx] = cimag(H_cs_val) / h;

            // Restore original value
            w_c[idx] = (complex double)w[idx];
        }
    }

    // Compute dwdt for active bodies only
    for (int ia = 0; ia < active.num_active; ia++) {
        int a = active.ids[ia];

        for (int k = 0; k < num_dim; k++) {
            int x_idx = a * num_dim + k;
            int p_idx = array_half + a * num_dim + k;

            dwdt[x_idx] += dHdw[p_idx];
            dwdt[p_idx] += -dHdw[x_idx];
        }
    }

    active_list_free(&active);
}
