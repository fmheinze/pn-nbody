/**
 * @file physical_quantities.c
 * @brief Functions for the computation of physical quantities.
 */

#include "eom.h"
#include "hamiltonian.h"
#include "cache.h"
#include "physical_quantities.h"
#include "pn_binary.h"

#include <float.h>
#include <math.h>


/**
 * @brief Computes the total energy of a system.
 *
 * Computes the total energy of a system based on the specified post-Newtonian approximation.
 *
 * @param[in]   w              Full state vector, w = [positions, momenta]
 * @param[in]   ode_params     Parameter struct containing general information about the system
 * @returns Total energy of the system.
 */
double total_energy_conservative(double* w, struct ode_params* ode_params)
{
    unsigned int levels = PAIR_CACHE_LEVEL_NONE;

    if (ode_params->pn_terms[0])
        levels |= PAIR_CACHE_LEVEL_GEOMETRY | PAIR_CACHE_LEVEL_P2;

    if (ode_params->pn_terms[1] || ode_params->pn_terms[2])
        levels |= PAIR_CACHE_LEVEL_GEOMETRY | PAIR_CACHE_LEVEL_P2 | PAIR_CACHE_LEVEL_PAIR_DOTS;

    PairCache *cache = pair_cache_get_workspace(ode_params);
    pair_cache_refresh(cache, w, ode_params, levels);

    double H = 0.0;
    if (ode_params->pn_terms[0] == 1)
        H += H0PN_cached(cache);

    if (ode_params->pn_terms[1] == 1)
        H += H1PN_cached(cache);

    if (ode_params->pn_terms[2] == 1) {
        if (ode_params->include_utt4)
            H += H2PN_cached(w, ode_params, cache, 1);
        else
            H += H2PN_cached(w, ode_params, cache, 0);
    }
    return H;
}


/**
 * @brief Prepare the pairwise kinematics shared by Newtonian osculating-element definitions.
 *
 * The supplied coordinate velocities may contain every enabled post-Newtonian order. Preparing
 * this state once lets independently selected orbital definitions reuse the same relative
 * position, relative velocity, separation and mass calculations.
 */
int compute_newtonian_relative_orbit_state(const double *w, const double *velocities,
    const struct ode_params *ode_params, int body_a, int body_b,
    NewtonianRelativeOrbitState *state)
{
    if (state == NULL)
        return 0;

    *state = (NewtonianRelativeOrbitState){0};

    if (w == NULL || velocities == NULL || ode_params == NULL)
        return 0;
    if (body_a < 0 || body_b < 0 || body_a >= ode_params->num_bodies_initial
        || body_b >= ode_params->num_bodies_initial || body_a == body_b)
        return 0;
    if (!ode_params->active[body_a] || !ode_params->active[body_b])
        return 0;

    const int D = ode_params->num_dim;
    if (D < 1 || D > 3)
        return 0;

    state->num_dim = D;
    state->total_mass = ode_params->masses[body_a] + ode_params->masses[body_b];
    if (!(state->total_mass > 0.0) || !isfinite(state->total_mass))
        return 0;

    double r2 = 0.0;
    for (int axis = 0; axis < D; ++axis) {
        const int idx_a = body_a*D + axis;
        const int idx_b = body_b*D + axis;
        state->relative_position[axis] = w[idx_a] - w[idx_b];
        state->relative_velocity[axis] = velocities[idx_a] - velocities[idx_b];
        r2 += state->relative_position[axis]*state->relative_position[axis];
        state->relative_speed_squared +=
            state->relative_velocity[axis]*state->relative_velocity[axis];
        state->radial_velocity_product +=
            state->relative_position[axis]*state->relative_velocity[axis];
    }

    if (!(r2 > 0.0) || !isfinite(r2) || !isfinite(state->relative_speed_squared)
        || !isfinite(state->radial_velocity_product))
        return 0;

    state->separation = sqrt(r2);
    state->valid = 1;
    return 1;
}


/**
 * @brief Return the Newtonian osculating semimajor axis from the specific orbital energy.
 *
 * Hyperbolic orbits have a negative value and an exactly parabolic state returns infinity.
 */
double semimajor_axis_newtonian(const NewtonianRelativeOrbitState *state)
{
    if (state == NULL || !state->valid)
        return NAN;

    const double specific_energy =
        0.5*state->relative_speed_squared - state->total_mass/state->separation;
    if (specific_energy == 0.0)
        return INFINITY;

    const double semimajor_axis = -state->total_mass/(2.0*specific_energy);
    return isnan(semimajor_axis) ? NAN : semimajor_axis;
}


/**
 * @brief Return the magnitude of the Newtonian osculating eccentricity vector.
 */
double eccentricity_newtonian(const NewtonianRelativeOrbitState *state)
{
    if (state == NULL || !state->valid)
        return NAN;

    const double radial_coefficient =
        state->relative_speed_squared - state->total_mass/state->separation;
    double eccentricity_squared = 0.0;
    for (int axis = 0; axis < state->num_dim; ++axis) {
        const double component =
            (radial_coefficient*state->relative_position[axis]
             - state->radial_velocity_product*state->relative_velocity[axis])
            / state->total_mass;
        eccentricity_squared += component*component;
    }

    return sqrt(eccentricity_squared);
}


/**
 * @brief Prepare radial quasi-Keplerian elements of an instantaneous isolated pair.
 *
 * The pair COM momentum is removed with the linear canonical Jacobi decomposition
 *
 *     P_rel = (m_b p_a - m_a p_b) / (m_a + m_b).
 *
 * From r and P_rel, the function evaluates the reduced conservative ADM binary Hamiltonian
 * through the enabled 1PN/2PN orders. The dissipative 2.5PN force and all other bodies are
 * deliberately excluded. The two roots of H(r, p_r=0, J) = E that enclose the current radius
 * define the instantaneous radial pericenter and apocenter.
 */
int compute_pn_radial_orbit_state(const double *w, const struct ode_params *ode_params,
    int body_a, int body_b, PNRadialOrbitState *state)
{
    if (state == NULL)
        return 0;

    *state = (PNRadialOrbitState){0};

    if (w == NULL || ode_params == NULL)
        return 0;
    if (body_a < 0 || body_b < 0 || body_a >= ode_params->num_bodies_initial
        || body_b >= ode_params->num_bodies_initial || body_a == body_b)
        return 0;
    if (!ode_params->active[body_a] || !ode_params->active[body_b])
        return 0;

    const int D = ode_params->num_dim;
    const int N = ode_params->num_bodies_initial;
    if (D < 2 || D > 3)
        return 0;

    const double mass_a = ode_params->masses[body_a];
    const double mass_b = ode_params->masses[body_b];
    const double total_mass = mass_a + mass_b;
    if (!(mass_a > 0.0) || !(mass_b > 0.0) || !isfinite(total_mass))
        return 0;

    const double reduced_mass = mass_a * mass_b / total_mass;
    const double nu = reduced_mass / total_mass;
    if (!(reduced_mass > 0.0) || !isfinite(reduced_mass)
        || !(nu > 0.0) || nu > 0.25)
        return 0;

    const int momentum_offset = N * D;
    double separation_squared = 0.0;
    double relative_momentum_squared = 0.0;
    double radial_momentum_product = 0.0;

    for (int axis = 0; axis < D; ++axis) {
        const int idx_a = body_a * D + axis;
        const int idx_b = body_b * D + axis;
        const double relative_position = w[idx_a] - w[idx_b];
        const double relative_momentum =
            (mass_b * w[momentum_offset + idx_a]
             - mass_a * w[momentum_offset + idx_b]) / total_mass;

        separation_squared += relative_position * relative_position;
        relative_momentum_squared += relative_momentum * relative_momentum;
        radial_momentum_product += relative_position * relative_momentum;
    }

    if (!(separation_squared > 0.0) || !isfinite(separation_squared)
        || !isfinite(relative_momentum_squared) || !isfinite(radial_momentum_product))
        return 0;

    double angular_momentum_squared =
        separation_squared * relative_momentum_squared
        - radial_momentum_product * radial_momentum_product;
    const double angular_scale =
        separation_squared * relative_momentum_squared
        + radial_momentum_product * radial_momentum_product;
    if (angular_momentum_squared < 0.0
        && angular_momentum_squared >= -64.0 * DBL_EPSILON * angular_scale)
        angular_momentum_squared = 0.0;
    if (!(angular_momentum_squared > 0.0) || !isfinite(angular_momentum_squared))
        return 0;

    const double separation = sqrt(separation_squared);
    const double x = separation / total_mass;
    const double pr_hat = radial_momentum_product / (separation * reduced_mass);
    const double j = sqrt(angular_momentum_squared) / (reduced_mass * total_mass);
    const int use_1pn = ode_params->pn_terms[1] == 1;
    const int use_2pn = ode_params->pn_terms[2] == 1;
    const double reduced_energy = pn_binary_reduced_hamiltonian(
        x, pr_hat, j, nu, use_1pn, use_2pn);

    double x_pericenter;
    double x_apocenter;
    if (!pn_binary_solve_turning_points(reduced_energy, j, nu, use_1pn, use_2pn, x,
            &x_pericenter, &x_apocenter))
        return 0;

    state->total_mass = total_mass;
    state->symmetric_mass_ratio = nu;
    state->reduced_energy = reduced_energy;
    state->reduced_angular_momentum = j;
    state->pericenter = total_mass * x_pericenter;
    state->apocenter = total_mass * x_apocenter;
    if (!isfinite(state->pericenter) || !isfinite(state->apocenter)
        || !(state->pericenter > 0.0) || state->apocenter < state->pericenter)
        return 0;

    state->valid = 1;
    return 1;
}


double pericenter_radial_pn(const PNRadialOrbitState *state)
{
    return state != NULL && state->valid ? state->pericenter : NAN;
}


double apocenter_radial_pn(const PNRadialOrbitState *state)
{
    return state != NULL && state->valid ? state->apocenter : NAN;
}


double semimajor_axis_radial_pn(const PNRadialOrbitState *state)
{
    if (state == NULL || !state->valid)
        return NAN;
    return 0.5 * (state->pericenter + state->apocenter);
}


double eccentricity_radial_pn(const PNRadialOrbitState *state)
{
    if (state == NULL || !state->valid)
        return NAN;

    const double radius_sum = state->pericenter + state->apocenter;
    if (!(radius_sum > 0.0) || !isfinite(radius_sum))
        return NAN;

    const double eccentricity = (state->apocenter - state->pericenter) / radius_sum;
    if (eccentricity < 0.0 && eccentricity > -64.0 * DBL_EPSILON)
        return 0.0;
    return eccentricity;
}
