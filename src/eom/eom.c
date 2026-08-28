/**
 * @file eom.c
 * @brief Public API and orchestration for the post-Newtonian equations of motion.
 */

#include <complex.h>
#include <string.h>

#include "cache.h"
#include "eom.h"
#include "eom_internal.h"
#include "utils.h"


static unsigned int pair_cache_levels_for_dynamics(const struct ode_params *params,
    int need_momentum, int need_chi_dot)
{
    unsigned int levels = PAIR_CACHE_LEVEL_NONE;

    if (params->pn_terms[0] || params->pn_terms[1] || params->pn_terms[2] || params->pn_terms[3])
        levels |= PAIR_CACHE_LEVEL_MASS_MOMENTUM;

    if ((params->pn_terms[0] && (need_momentum || need_chi_dot))
        || params->pn_terms[1] || params->pn_terms[2] || params->pn_terms[3])
        levels |= PAIR_CACHE_LEVEL_GEOMETRY;

    if (params->pn_terms[1] || params->pn_terms[2])
        levels |= PAIR_CACHE_LEVEL_P2 | PAIR_CACHE_LEVEL_PAIR_DOTS;

    return levels;
}


// ------------------------------------------------------------------------------------------------
// Main EOM right-hand side
// ------------------------------------------------------------------------------------------------

// Evaluate only the requested halves of the analytic PN dynamics. Lower-order derivatives are
// additionally evaluated when they are physically required to construct the 2.5PN chi-dot tensor.
static void evaluate_pn_dynamics(double* w, struct ode_params* ode_params, double* velocity,
    double* momentum, DynamicsCache *dynamics_cache)
{
    // --------------------------------------------------------------------------------------------
    // Initialize the arrays
    // --------------------------------------------------------------------------------------------

    const int num_bodies = ode_params->num_bodies_initial;
    const int num_dim = ode_params->num_dim;
    const int array_half = num_bodies * num_dim;
    const int need_velocity = velocity != NULL;
    const int need_momentum = momentum != NULL;

    const int need_chi_dot = ode_params->pn_terms[3] && !dynamics_cache->chi_dot_valid;
    PairCache *cache = pair_cache_get_workspace(ode_params);
    pair_cache_refresh(cache, w, ode_params,
        pair_cache_levels_for_dynamics(ode_params, need_momentum, need_chi_dot));
    const ActiveList *active = &cache->active;

    double *lower_xdot = NULL;
    double *lower_pdot = NULL;
    if (need_chi_dot) {
        double *const lower_rhs = dynamics_cache_get_scratch(dynamics_cache);
        memset(lower_rhs, 0, (size_t)2 * (size_t)array_half * sizeof(*lower_rhs));
        lower_xdot = lower_rhs;
        lower_pdot = lower_rhs + array_half;
    }

    if (need_velocity)
        memset(velocity, 0, (size_t)array_half * sizeof(*velocity));
    if (need_momentum)
        memset(momentum, 0, (size_t)array_half * sizeof(*momentum));

    // --------------------------------------------------------------------------------------------
    // Add 0PN (Newtonian) terms
    // --------------------------------------------------------------------------------------------

    if (ode_params->pn_terms[0]) {
        if (need_velocity || need_chi_dot)
            add_velocity_0pn(cache, velocity, lower_xdot);

        if (need_momentum || need_chi_dot)
            add_momentum_0pn(cache, momentum, lower_pdot);
    }

    // --------------------------------------------------------------------------------------------
    // Add 1PN terms
    // --------------------------------------------------------------------------------------------

    if (ode_params->pn_terms[1]) {
        if (need_velocity || need_chi_dot)
            add_velocity_1pn(cache, velocity, lower_xdot);

        if (need_momentum || need_chi_dot)
            add_momentum_1pn(cache, momentum, lower_pdot);
    }

    // chi_dot depends only on the 0PN+1PN dynamics. Build it before 2PN so the same persistent
    // scratch can subsequently be reused for a UTT4 gradient.
    if (need_chi_dot) {
        if (num_dim != 3)
            errorexit("2.5PN dynamics currently assumes num_dim = 3");

        pn_compute_25pn_chi_dot(cache, num_bodies, num_dim,
            (double (*)[num_dim])lower_xdot, (double (*)[num_dim])lower_pdot,
            dynamics_cache->chi_dot);
        dynamics_cache->chi_dot_valid = 1;
    }

    // --------------------------------------------------------------------------------------------
    // Add 2PN terms
    // --------------------------------------------------------------------------------------------

    if (ode_params->pn_terms[2]) {
        if (need_velocity) {
            add_velocity_2pn_onebody_analytic(cache, velocity);
            add_velocity_2pn_pair_analytic(cache, velocity);
            add_velocity_2pn_triple_analytic(cache, velocity);
        }

        if (need_momentum) {
            add_momentum_2pn_pair_analytic(cache, momentum);
            add_momentum_2pn_triple_analytic(cache, momentum);
            add_momentum_2pn_fourbody_analytic(cache, momentum);
        }

        // If not using impulse splitting, add UTT4 contributions directly to dp/dt.
        if (need_momentum && ode_params->include_utt4
            && !ode_params->use_impulse_method)
        {
            double *const dUdx = dynamics_cache_get_scratch(dynamics_cache);

            compute_dUTT4_dx(w, ode_params, dUdx);

            for (int ia = 0; ia < active->num_active; ia++) {
                const int a = active->ids[ia];

                for (int i = 0; i < num_dim; i++) {
                    const int idx = a * num_dim + i;
                    momentum[idx] -= dUdx[idx];
                }
            }
        }
    }

    // --------------------------------------------------------------------------------------------
    // Add 2.5PN terms
    // --------------------------------------------------------------------------------------------

    if (ode_params->pn_terms[3]) {
        if (num_dim != 3)
            errorexit("2.5PN dynamics currently assumes num_dim = 3");

        if (need_velocity)
            pn_add_25pn_velocity_contribution(cache, num_dim, dynamics_cache->chi_dot, velocity);

        if (need_momentum)
            pn_add_25pn_momentum_contribution(cache, num_dim, dynamics_cache->chi_dot, momentum);
    }
}


/**
 * @brief Compute the coordinate velocities dx/dt at the enabled PN orders.
 *
 * This evaluates dx/dt = dH/dp for the conservative 0PN, 1PN, and 2PN terms and includes the
 * explicit dissipative 2.5PN velocity contribution when enabled. Position-only 2PN terms,
 * including UTT4, are deliberately skipped. Results are cached for the exact phase-space state so
 * output and other derived quantities share the evaluation, and a following RHS call computes
 * only a missing momentum half.
 *
 * @param[in]       w           State of the system, w = [positions, momenta]
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[out]      velocities  Coordinate velocities, one component per position component
 */
void compute_coordinate_velocities(double* w, struct ode_params* ode_params, double* velocities)
{
    const int array_half = ode_params->num_bodies_initial * ode_params->num_dim;
    DynamicsCache *cache = dynamics_cache_prepare(ode_params, w);

    if (!cache->velocity_valid) {
        evaluate_pn_dynamics(w, ode_params, cache->velocity, NULL, cache);
        cache->velocity_valid = 1;
    }

    memcpy(velocities, cache->velocity, (size_t)array_half*sizeof(*velocities));
}


/**
 * @brief Right-hand side of the N-body post-Newtonian equations of motion
 *
 * @param[in]       t           Time (currently not used, but kept for completeness)
 * @param[in]       w           State of the system, w = [positions, momenta]
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[out]      dwdt        Right-hand side of the equations of motion
 */
void rhs_pn_nbody(double t, double* w, struct ode_params* ode_params, double* dwdt)
{
    (void)t;
    const int array_half = ode_params->num_bodies_initial * ode_params->num_dim;
    DynamicsCache *cache = dynamics_cache_prepare(ode_params, w);

    if (cache->rhs_valid) {
        memcpy(dwdt, cache->rhs, (size_t)(2*array_half)*sizeof(*dwdt));
        return;
    }

    if (cache->velocity_valid) {
        evaluate_pn_dynamics(w, ode_params, NULL, dwdt + array_half, cache);
        memcpy(dwdt, cache->velocity, (size_t)array_half*sizeof(*dwdt));
    } else {
        evaluate_pn_dynamics(w, ode_params, dwdt, dwdt + array_half, cache);
        memcpy(cache->velocity, dwdt, (size_t)array_half*sizeof(*cache->velocity));
        cache->velocity_valid = 1;
    }

    memcpy(cache->rhs, dwdt, (size_t)(2*array_half)*sizeof(*cache->rhs));
    cache->rhs_valid = 1;
}


// ------------------------------------------------------------------------------------------------
// Complex-step differentiation of Hamiltonian (mainly for cross-checks / validation)
// ------------------------------------------------------------------------------------------------

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

    // Copy original array to w_c and initialize dHdw
    for (int i = 0; i < total_dim; ++i) {
        w_c[i] = (complex double)w[i];
        dHdw[i] = 0.0;
    }

    // Position derivatives: dH/dx. These are needed for dp/dt = -dH/dx.
    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a])
            continue;

        for (int k = 0; k < num_dim; k++) {
            int idx = a * num_dim + k;

            // Add tiny imaginary step in coordinate idx
            w_c[idx] += (complex double)I * h;

            // p_flag = 0: keep position-only Hamiltonian terms
            H_cs_val = H(w_c, ode_params, 0);
            dHdw[idx] = cimag(H_cs_val) / h;

            // Restore original value
            w_c[idx] = (complex double)w[idx];
        }
    }

    // Momentum derivatives: dH/dp. These are needed for dx/dt = dH/dp.
    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a])
            continue;

        for (int k = 0; k < num_dim; k++) {
            int idx = array_half + a * num_dim + k;

            // Add tiny imaginary step in momentum component idx
            w_c[idx] += (complex double)I * h;

            // p_flag = 1: skip Hamiltonian terms without momentum dependence
            H_cs_val = H(w_c, ode_params, 1);
            dHdw[idx] = cimag(H_cs_val) / h;

            // Restore original value
            w_c[idx] = (complex double)w[idx];
        }
    }

    // Compute dwdt for active bodies only
    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a])
            continue;

        for (int k = 0; k < num_dim; k++) {
            int x_idx = a * num_dim + k;
            int p_idx = array_half + a * num_dim + k;

            dwdt[x_idx] += dHdw[p_idx];
            dwdt[p_idx] += -dHdw[x_idx];
        }
    }
}
