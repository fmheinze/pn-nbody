/**
 * @file eom.c
 * @brief Functions for the right-hand side of the post-Newtonian equations of motion
 */

#include <complex.h>
#include "utils.h"
#include "eom.h"
#include "hamiltonian.h"


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
    int num_bodies = ode_params->num_bodies;
    int num_dim = ode_params->num_dim;
    int array_half = num_bodies * num_dim; 
    double result, temp;

    // Masses
    double m[num_bodies];
    double inv_m[num_bodies];
    for (int a = 0; a < num_bodies; a++) {
        m[a] = ode_params->masses[a];
        inv_m[a] = 1.0 / m[a];
    }
    
    // Momenta
    double p[num_bodies][num_dim];
    for (int a = 0; a < num_bodies; a++)
        for (int i = 0; i < num_dim; i++)
            p[a][i] = w[array_half + a * num_dim + i];
    
    // Relative positions and distances
    double x_rel[num_bodies][num_bodies][num_dim]; 
    double n[num_bodies][num_bodies][num_dim];
    double r[num_bodies][num_bodies];
    double inv_r[num_bodies][num_bodies];
    for (int a = 0; a < num_bodies; a++) {
        for (int b = a; b < num_bodies; b++) {
            for (int i = 0; i < num_dim; i++){
                x_rel[a][b][i] = w[a * num_dim + i] - w[b * num_dim + i];
                x_rel[b][a][i] = -x_rel[a][b][i];
            } 
            if (a == b) {
                r[a][b] = 0.0;
                inv_r[a][b] = 0.0;
            } else {
                r[a][b] = norm(x_rel[a][b], num_dim);
                inv_r[a][b] = 1.0 / r[a][b];
            }
            r[b][a] = r[a][b];
            inv_r[b][a] = inv_r[a][b];

            for (int i = 0; i < num_dim; i++){
                if (a == b){
                        n[a][b][i] = n[b][a][i] = 0.0;
                } else {
                        n[a][b][i] = x_rel[a][b][i] * inv_r[a][b];
                        n[b][a][i] = -n[a][b][i];
                }
            }
        }
    }

    // Time derivatives
    double p_dot[num_bodies][num_dim];
    double x_dot[num_bodies][num_dim];
    for (int a = 0; a < num_bodies; a++) {
        for (int i = 0; i < num_dim; i++) {
            p_dot[a][i] = 0.0;
            x_dot[a][i] = 0.0;
        }
    }
            
    // Set ODE right-hand side initially to zero
    for (int i = 0; i < 2 * array_half; i++)   
        dwdt[i] = 0.0;

    // --------------------------------------------------------------------------------------------
    // Add 0PN (Newtonian) terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[0]) {
        for (int a = 0; a < num_bodies; a++) {
            // Velocities
            for (int i = 0; i < num_dim; i++) {
                result = w[array_half + a * num_dim + i] * inv_m[a];;
                dwdt[a * num_dim + i] += result;
                x_dot[a][i] += result;
            }

            // Accelerations
            for (int b = a + 1; b < num_bodies; b++) {
                for (int i = 0; i < num_dim; i++) {
                    result = - m[a] * m[b] * inv_r[a][b] * inv_r[a][b] * n[a][b][i];
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
        for (int a = 0; a < num_bodies; a++) {
            double pa_dot_pa = dot_product(p[a], p[a], num_dim);
            for (int i = 0; i < num_dim; i++) {

                // Velocities
                result = -0.5 * pa_dot_pa * inv_m[a] * inv_m[a] * inv_m[a] * p[a][i];
                x_dot[a][i] += result;
                dwdt[a * num_dim + i] += result;
                for (int b = 0; b < num_bodies; b++) {
                    double pa_dot_pb = dot_product(p[a], p[b], num_dim);
                    double pb_dot_pb = dot_product(p[b], p[b], num_dim);
                    double nab_dot_pa = dot_product(n[a][b], p[a], num_dim);
                    double nab_dot_pb = dot_product(n[a][b], p[b], num_dim);

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
                    for (int c = 0; c < num_bodies; c++) {
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
        if (ode_params->include_utt4 && !ode_params->use_impulse_method)
        {
            double dUdx[array_half];
            compute_dUTT4_dx(w, ode_params, dUdx);
            for (int i = 0; i < array_half; i++)
                dwdt[array_half + i] -= dUdx[i];
        }
    }

    // --------------------------------------------------------------------------------------------
    // Add 2.5PN terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[3]) {

        // Initialize arrays
        double x_rel_dot[num_bodies][num_bodies][num_dim];

        for (int a = 0; a < num_bodies; a++)
            for (int b = 0; b < num_bodies; b++)
                for (int i = 0; i < num_dim; i++)
                    x_rel_dot[a][b][i] = x_dot[a][i] - x_dot[b][i];

        double chi_dot[num_dim][num_dim];
        double dp_chi[num_bodies][num_dim][num_dim][num_dim];
        double dx_chi[num_bodies][num_dim][num_dim][num_dim];

        for (int i = 0; i < num_dim; i++)
            for (int j = 0; j < num_dim; j++)
                chi_dot[i][j] = 0.0;

        for (int a = 0; a < num_bodies; a++) { 
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
                for (int a = 0; a < num_bodies; a++) {
                    chi_dot[i][j] += 2.0 * inv_m[a] * (
                        2 * dot_product(p_dot[a], p[a], num_dim) * delta(i, j) 
                        - 3 * (p_dot[a][i] * p[a][j] + p[a][i] * p_dot[a][j]) );
                    
                    for (int b = 0; b < num_bodies; b++) {
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

        for (int c = 0; c < num_bodies; c++) {
            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        dp_chi[c][i][j][k] += 2.0 / m[c] * (2 * p[c][k] * delta(i, j) 
                            - 3 * (p[c][j] * delta(i, k) + p[c][i] * delta(j, k)));
                        
                        for (int a = 0; a < num_bodies; a++) {
                            for (int b = 0; b < num_bodies; b++) {
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
        for (int a = 0; a < num_bodies; a++) {
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
    int array_half = ode_params->num_dim * ode_params->num_bodies;
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
        H_cs_val = H(w_c, ode_params, p_flag);
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
