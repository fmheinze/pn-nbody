/**
 * @file physical_quantities.c
 * @brief Functions for preparing and evaluating physical quantities.
 *
 * System-wide quantities are evaluated directly from the active canonical state. Pair quantities
 * are evaluated through small prepared-state structs built for one selected pair at one
 * phase-space point. The base helpers collect shared relative kinematics or canonical PN invariants
 * once; optional helpers then add turning points, frequencies, vectors, or time derivatives only
 * when a requested output needs them. The individual quantity functions consume these prepared
 * states, which avoids repeating expensive work such as a turning-point solve or a full
 * right-hand-side evaluation and centralizes validity and NAN handling.
 */

#include "eom.h"
#include "hamiltonian.h"
#include "cache.h"
#include "physical_quantities.h"
#include "pn_binary.h"

#include <float.h>
#include <math.h>


// ------------------------------------------------------------------------------------------------
// Shared helpers
// ------------------------------------------------------------------------------------------------

// Solve for the radial turning points enclosing the state's current separation.
int compute_pn_radial_turning_points(PNOrbitState *state)
{
    if (state == NULL || !state->valid)
        return 0;

    state->radial_turning_points_valid = 0;
    state->pericenter = 0.0;
    state->apocenter = 0.0;

    double x_pericenter;
    double x_apocenter;
    if (!pn_binary_solve_turning_points(state->reduced_energy,
            state->reduced_angular_momentum, state->symmetric_mass_ratio,
            state->use_1pn, state->use_2pn, state->reduced_separation,
            &x_pericenter, &x_apocenter))
        return 0;

    state->pericenter = state->total_mass * x_pericenter;
    state->apocenter = state->total_mass * x_apocenter;
    if (!isfinite(state->pericenter) || !isfinite(state->apocenter)
        || !(state->pericenter > 0.0) || state->apocenter < state->pericenter)
        return 0;

    state->radial_turning_points_valid = 1;
    return 1;
}


static int pn_orbit_bound_invariants(const PNOrbitState *state, long double *energy_scale,
    long double *energy_angular_momentum_squared, long double *nu)
{
    if (state == NULL || !state->valid || energy_scale == NULL
        || energy_angular_momentum_squared == NULL || nu == NULL)
        return 0;

    *energy_scale = -2.0L * state->reduced_energy;
    const long double angular_momentum = state->reduced_angular_momentum;
    *energy_angular_momentum_squared =
        *energy_scale * angular_momentum * angular_momentum;
    *nu = state->symmetric_mass_ratio;

    return *energy_scale > 0.0L && *energy_angular_momentum_squared > 0.0L
        && isfinite(*energy_scale) && isfinite(*energy_angular_momentum_squared)
        && *nu > 0.0L && *nu <= 0.25L && isfinite(*nu);
}


// Return the dimensionless mean motion M*n from the ADM quasi-Keplerian expansion through 2PN.
static long double pn_orbit_dimensionless_mean_motion(const PNOrbitState *state,
    long double energy_scale, long double energy_angular_momentum_squared, long double nu)
{
    long double expansion = 1.0L;
    if (state->use_1pn)
        expansion += energy_scale / 8.0L * (-15.0L + nu);
    if (state->use_2pn)
        expansion += energy_scale*energy_scale / 128.0L
            * (555.0L + 30.0L*nu + 11.0L*nu*nu
               + 192.0L/sqrtl(energy_angular_momentum_squared)*(-5.0L + 2.0L*nu));

    const long double mean_motion = energy_scale*sqrtl(energy_scale)*expansion;
    return isfinite(mean_motion) && mean_motion > 0.0L ? mean_motion : NAN;
}


// Return K = Phi/(2*pi), the azimuth accumulated per radial period, through 2PN.
static long double pn_orbit_periastron_advance_factor(const PNOrbitState *state,
    long double energy_scale, long double energy_angular_momentum_squared, long double nu)
{
    long double factor = 1.0L;
    if (state->use_1pn)
        factor += 3.0L*energy_scale/energy_angular_momentum_squared;
    if (state->use_2pn)
        factor += energy_scale*energy_scale / 4.0L
            * (3.0L/energy_angular_momentum_squared*(-5.0L + 2.0L*nu)
               + 15.0L/(energy_angular_momentum_squared*energy_angular_momentum_squared)
                 *(7.0L - 2.0L*nu));

    return isfinite(factor) && factor > 0.0L ? factor : NAN;
}


// ------------------------------------------------------------------------------------------------
// State preparation
// ------------------------------------------------------------------------------------------------

/**
 * @brief Prepare the origin-based and center-of-mass-frame angular momenta of the active system.
 *
 * The orbital part is sum_i r_i x p_i about the coordinate origin. Stored spins are dimensionless
 * Kerr vectors chi_i, so their physical contribution in G = c = 1 units is m_i^2 chi_i. The
 * intrinsic orbital part uses the mass-weighted COM convention R = sum_i m_i r_i / sum_i m_i and
 * L_COM = L - R x P. Every scalar/vector output reuses this single active-body traversal.
 */
int compute_system_angular_momentum_state(const double *w,
    const struct ode_params *ode_params, SystemAngularMomentumState *state)
{
    if (state == NULL)
        return 0;
    *state = (SystemAngularMomentumState){0};

    if (w == NULL || ode_params == NULL || ode_params->active == NULL
        || ode_params->masses == NULL)
        return 0;

    const int D = ode_params->num_dim;
    const int N = ode_params->num_bodies_initial;
    if (D < 1 || D > 3 || N < 1)
        return 0;

    const int momentum_offset = D*N;
    const int spin_offset = 2*momentum_offset;
    long double orbital[3] = {0.0L, 0.0L, 0.0L};
    long double spin[3] = {0.0L, 0.0L, 0.0L};
    long double mass_position[3] = {0.0L, 0.0L, 0.0L};
    long double total_momentum[3] = {0.0L, 0.0L, 0.0L};
    long double total_mass = 0.0L;

    for (int body = 0; body < N; ++body) {
        if (!ode_params->active[body])
            continue;

        const double mass = ode_params->masses[body];
        if (!(mass > 0.0) || !isfinite(mass))
            return 0;

        long double position[3] = {0.0L, 0.0L, 0.0L};
        long double momentum[3] = {0.0L, 0.0L, 0.0L};
        for (int axis = 0; axis < D; ++axis) {
            const int component = body*D + axis;
            position[axis] = w[component];
            momentum[axis] = w[momentum_offset + component];
            if (!isfinite(position[axis]) || !isfinite(momentum[axis]))
                return 0;
        }

        orbital[0] += position[1]*momentum[2] - position[2]*momentum[1];
        orbital[1] += position[2]*momentum[0] - position[0]*momentum[2];
        orbital[2] += position[0]*momentum[1] - position[1]*momentum[0];

        total_mass += mass;
        for (int axis = 0; axis < 3; ++axis) {
            mass_position[axis] += (long double)mass*position[axis];
            total_momentum[axis] += momentum[axis];
        }

        const long double mass_squared = (long double)mass*mass;
        for (int axis = 0; axis < 3; ++axis) {
            const double dimensionless_spin = w[spin_offset + 3*body + axis];
            if (!isfinite(dimensionless_spin))
                return 0;
            spin[axis] += mass_squared*dimensionless_spin;
        }
    }

    if (!(total_mass > 0.0L) || !isfinite(total_mass))
        return 0;

    long double center_of_mass[3];
    for (int axis = 0; axis < 3; ++axis)
        center_of_mass[axis] = mass_position[axis]/total_mass;

    const long double com_cross_total_momentum[3] = {
        center_of_mass[1]*total_momentum[2] - center_of_mass[2]*total_momentum[1],
        center_of_mass[2]*total_momentum[0] - center_of_mass[0]*total_momentum[2],
        center_of_mass[0]*total_momentum[1] - center_of_mass[1]*total_momentum[0]
    };

    long double orbital_com[3];
    long double total[3];
    long double total_com[3];
    for (int axis = 0; axis < 3; ++axis) {
        orbital_com[axis] = orbital[axis] - com_cross_total_momentum[axis];
        total[axis] = orbital[axis] + spin[axis];
        total_com[axis] = orbital_com[axis] + spin[axis];
    }

    const long double orbital_norm = sqrtl(
        orbital[0]*orbital[0] + orbital[1]*orbital[1] + orbital[2]*orbital[2]);
    const long double orbital_com_norm = sqrtl(
        orbital_com[0]*orbital_com[0] + orbital_com[1]*orbital_com[1]
        + orbital_com[2]*orbital_com[2]);
    const long double spin_norm = sqrtl(
        spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
    const long double total_norm = sqrtl(
        total[0]*total[0] + total[1]*total[1] + total[2]*total[2]);
    const long double total_com_norm = sqrtl(
        total_com[0]*total_com[0] + total_com[1]*total_com[1]
        + total_com[2]*total_com[2]);

    for (int axis = 0; axis < 3; ++axis) {
        if (!isfinite(orbital[axis]) || !isfinite(orbital_com[axis])
            || !isfinite(spin[axis]) || !isfinite(total[axis]) || !isfinite(total_com[axis])
            || fabsl(orbital[axis]) > DBL_MAX || fabsl(spin[axis]) > DBL_MAX
            || fabsl(orbital_com[axis]) > DBL_MAX || fabsl(total[axis]) > DBL_MAX
            || fabsl(total_com[axis]) > DBL_MAX)
            return 0;
    }
    if (!isfinite(orbital_norm) || !isfinite(orbital_com_norm) || !isfinite(spin_norm)
        || !isfinite(total_norm) || !isfinite(total_com_norm)
        || orbital_norm > DBL_MAX || orbital_com_norm > DBL_MAX || spin_norm > DBL_MAX
        || total_norm > DBL_MAX || total_com_norm > DBL_MAX)
        return 0;

    for (int axis = 0; axis < 3; ++axis) {
        state->orbital_vector[axis] = (double)orbital[axis];
        state->orbital_com_vector[axis] = (double)orbital_com[axis];
        state->spin_vector[axis] = (double)spin[axis];
        state->total_vector[axis] = (double)total[axis];
        state->total_com_vector[axis] = (double)total_com[axis];
    }
    state->orbital_magnitude = (double)orbital_norm;
    state->orbital_com_magnitude = (double)orbital_com_norm;
    state->spin_magnitude = (double)spin_norm;
    state->total_magnitude = (double)total_norm;
    state->total_com_magnitude = (double)total_com_norm;
    state->valid = 1;
    return 1;
}


/**
 * @brief Prepare the pairwise kinematics shared by Newtonian osculating-element definitions.
 *
 * The optional coordinate velocities may contain every enabled post-Newtonian order. Preparing
 * this state once lets independently selected orbital definitions reuse the same relative
 * position, separation and, when supplied, relative velocity calculations. A position-only state
 * is sufficient for separation output and avoids an otherwise unnecessary dynamics evaluation.
 */
int compute_newtonian_relative_orbit_state(const double *w, const double *velocities,
    const struct ode_params *ode_params, int body_a, int body_b,
    NewtonianRelativeOrbitState *state)
{
    if (state == NULL)
        return 0;

    *state = (NewtonianRelativeOrbitState){0};

    if (w == NULL || ode_params == NULL)
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
        r2 += state->relative_position[axis]*state->relative_position[axis];

        if (velocities != NULL) {
            state->relative_velocity[axis] = velocities[idx_a] - velocities[idx_b];
            state->relative_speed_squared +=
                state->relative_velocity[axis]*state->relative_velocity[axis];
            state->radial_velocity_product +=
                state->relative_position[axis]*state->relative_velocity[axis];
        }
    }

    if (!(r2 > 0.0) || !isfinite(r2))
        return 0;
    if (velocities != NULL
        && (!isfinite(state->relative_speed_squared)
            || !isfinite(state->radial_velocity_product)))
        return 0;

    state->separation = sqrt(r2);
    state->velocity_valid = velocities != NULL;
    state->valid = 1;
    return 1;
}


/**
 * @brief Prepare the invariants used by quasi-Keplerian elements of an isolated pair.
 *
 * The pair COM momentum is removed with the linear canonical Jacobi decomposition
 *
 *     P_rel = (m_b p_a - m_a p_b) / (m_a + m_b).
 *
 * From r and P_rel, the function evaluates the reduced conservative ADM binary Hamiltonian
 * through the enabled 1PN/2PN orders. The dissipative 2.5PN force and all other bodies are
 * deliberately excluded. Radial turning points are prepared separately so invariant preparation
 * and the numerical root solve remain distinct and can be shared selectively.
 */
int compute_pn_orbit_state(const double *w, const struct ode_params *ode_params,
    int body_a, int body_b, PNOrbitState *state)
{
    if (state == NULL)
        return 0;

    *state = (PNOrbitState){0};

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
    double relative_position[3] = {0.0, 0.0, 0.0};
    double relative_momentum[3] = {0.0, 0.0, 0.0};
    double separation_squared = 0.0;
    double relative_momentum_squared = 0.0;
    double radial_momentum_product = 0.0;

    for (int axis = 0; axis < D; ++axis) {
        const int idx_a = body_a * D + axis;
        const int idx_b = body_b * D + axis;
        relative_position[axis] = w[idx_a] - w[idx_b];
        relative_momentum[axis] =
            (mass_b * w[momentum_offset + idx_a]
             - mass_a * w[momentum_offset + idx_b]) / total_mass;

        separation_squared += relative_position[axis] * relative_position[axis];
        relative_momentum_squared += relative_momentum[axis] * relative_momentum[axis];
        radial_momentum_product += relative_position[axis] * relative_momentum[axis];
    }

    if (!(separation_squared > 0.0) || !isfinite(separation_squared)
        || !isfinite(relative_momentum_squared) || !isfinite(radial_momentum_product))
        return 0;

    double angular_momentum_vector[3];
    angular_momentum_vector[0] =
        relative_position[1]*relative_momentum[2]
        - relative_position[2]*relative_momentum[1];
    angular_momentum_vector[1] =
        relative_position[2]*relative_momentum[0]
        - relative_position[0]*relative_momentum[2];
    angular_momentum_vector[2] =
        relative_position[0]*relative_momentum[1]
        - relative_position[1]*relative_momentum[0];

    double angular_momentum_squared = 0.0;
    for (int axis = 0; axis < 3; ++axis)
        angular_momentum_squared += angular_momentum_vector[axis]*angular_momentum_vector[axis];
    if (!isfinite(angular_momentum_squared))
        return 0;

    const double separation = sqrt(separation_squared);
    const double x = separation / total_mass;
    const double pr_hat = radial_momentum_product / (separation * reduced_mass);
    const double j = sqrt(angular_momentum_squared) / (reduced_mass * total_mass);
    const int use_1pn = ode_params->pn_terms[1] == 1;
    const int use_2pn = ode_params->pn_terms[2] == 1;
    const double reduced_energy = pn_binary_reduced_hamiltonian(
        x, pr_hat, j, nu, use_1pn, use_2pn);
    if (!isfinite(reduced_energy))
        return 0;

    state->use_1pn = use_1pn;
    state->use_2pn = use_2pn;
    state->total_mass = total_mass;
    state->symmetric_mass_ratio = nu;
    state->reduced_energy = reduced_energy;
    state->reduced_angular_momentum = j;
    state->reduced_separation = x;
    for (int axis = 0; axis < 3; ++axis) {
        state->relative_position[axis] = relative_position[axis];
        state->relative_momentum[axis] = relative_momentum[axis];
        state->angular_momentum_vector[axis] = angular_momentum_vector[axis];
    }
    state->valid = 1;
    return 1;
}


/**
 * @brief Prepare the radial and orbit-averaged azimuthal frequencies of an isolated PN pair.
 *
 * The quasi-Keplerian mean motion gives the physical radial frequency Omega_r = n. During one
 * radial period the azimuth advances by Phi = 2*pi*K, so Omega_phi = K*n. Preparing this state once
 * lets every frequency-, period-, and precession-based output reuse the same 2PN evaluation.
 */
int compute_pn_orbit_frequency_state(const PNOrbitState *orbit_state,
    PNOrbitFrequencyState *frequency_state)
{
    if (frequency_state == NULL)
        return 0;
    *frequency_state = (PNOrbitFrequencyState){0};

    long double energy_scale;
    long double energy_angular_momentum_squared;
    long double nu;
    if (!pn_orbit_bound_invariants(
            orbit_state, &energy_scale, &energy_angular_momentum_squared, &nu)
        || !(orbit_state->total_mass > 0.0) || !isfinite(orbit_state->total_mass))
        return 0;

    const long double dimensionless_mean_motion = pn_orbit_dimensionless_mean_motion(
        orbit_state, energy_scale, energy_angular_momentum_squared, nu);
    const long double periastron_factor = pn_orbit_periastron_advance_factor(
        orbit_state, energy_scale, energy_angular_momentum_squared, nu);
    const long double radial_frequency =
        dimensionless_mean_motion / orbit_state->total_mass;
    const long double azimuthal_frequency = radial_frequency * periastron_factor;
    if (!(radial_frequency > 0.0L) || !(azimuthal_frequency > 0.0L)
        || !isfinite(radial_frequency) || !isfinite(azimuthal_frequency)
        || !isfinite(periastron_factor))
        return 0;

    frequency_state->total_mass = orbit_state->total_mass;
    frequency_state->radial_frequency = (double)radial_frequency;
    frequency_state->azimuthal_frequency = (double)azimuthal_frequency;
    frequency_state->periastron_advance_factor = (double)periastron_factor;
    frequency_state->valid = 1;
    return 1;
}


/**
 * @brief Differentiate the isolated-pair energy and angular-momentum magnitude along the full RHS.
 *
 * The diagnostic quantities are functions of the pair's relative canonical state, while rhs is
 * the derivative generated by the complete active N-body system, including external bodies and
 * 2.5PN dissipation. Consequently these rates vanish for an isolated conservative pair (up to
 * roundoff) but retain energy and angular-momentum exchange with the environment.
 */
int compute_pn_orbit_rate_state(const double *rhs, const struct ode_params *ode_params,
    int body_a, int body_b, const PNOrbitState *orbit_state, PNOrbitRateState *rate_state)
{
    if (rate_state == NULL)
        return 0;
    *rate_state = (PNOrbitRateState){0};

    if (rhs == NULL || ode_params == NULL || orbit_state == NULL || !orbit_state->valid)
        return 0;
    if (body_a < 0 || body_b < 0 || body_a >= ode_params->num_bodies_initial
        || body_b >= ode_params->num_bodies_initial || body_a == body_b)
        return 0;

    const int D = ode_params->num_dim;
    const int momentum_offset = ode_params->num_bodies_initial*D;
    if (D < 2 || D > 3)
        return 0;

    const double mass_a = ode_params->masses[body_a];
    const double mass_b = ode_params->masses[body_b];
    const double total_mass = mass_a + mass_b;
    if (!(mass_a > 0.0) || !(mass_b > 0.0) || !isfinite(total_mass))
        return 0;
    const double reduced_mass = mass_a*mass_b/total_mass;
    const double separation = orbit_state->reduced_separation*total_mass;
    if (!(total_mass > 0.0) || !(reduced_mass > 0.0) || !(separation > 0.0)
        || !isfinite(reduced_mass) || !isfinite(separation))
        return 0;

    double relative_velocity[3] = {0.0, 0.0, 0.0};
    double relative_momentum_rate[3] = {0.0, 0.0, 0.0};
    double radial_velocity = 0.0;
    double radial_momentum = 0.0;
    for (int axis = 0; axis < D; ++axis) {
        const int idx_a = body_a*D + axis;
        const int idx_b = body_b*D + axis;
        const double radial_unit = orbit_state->relative_position[axis]/separation;

        relative_velocity[axis] = rhs[idx_a] - rhs[idx_b];
        relative_momentum_rate[axis] =
            (mass_b*rhs[momentum_offset + idx_a]
             - mass_a*rhs[momentum_offset + idx_b])/total_mass;
        radial_velocity += radial_unit*relative_velocity[axis];
        radial_momentum += radial_unit*orbit_state->relative_momentum[axis];
    }

    const double pr_hat = radial_momentum/reduced_mass;
    double pr_hat_rate = 0.0;
    for (int axis = 0; axis < D; ++axis) {
        const double radial_unit = orbit_state->relative_position[axis]/separation;
        const double radial_unit_rate =
            (relative_velocity[axis] - radial_unit*radial_velocity)/separation;
        pr_hat_rate +=
            (radial_unit_rate*orbit_state->relative_momentum[axis]
             + radial_unit*relative_momentum_rate[axis])/reduced_mass;
    }

    const double *r = orbit_state->relative_position;
    const double *p = orbit_state->relative_momentum;
    const double angular_momentum_rate_vector[3] = {
        relative_velocity[1]*p[2] - relative_velocity[2]*p[1]
            + r[1]*relative_momentum_rate[2] - r[2]*relative_momentum_rate[1],
        relative_velocity[2]*p[0] - relative_velocity[0]*p[2]
            + r[2]*relative_momentum_rate[0] - r[0]*relative_momentum_rate[2],
        relative_velocity[0]*p[1] - relative_velocity[1]*p[0]
            + r[0]*relative_momentum_rate[1] - r[1]*relative_momentum_rate[0]
    };

    double angular_momentum_squared = 0.0;
    double angular_momentum_rate_squared = 0.0;
    double angular_momentum_dot_rate = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        const double component = orbit_state->angular_momentum_vector[axis];
        angular_momentum_squared += component*component;
        angular_momentum_rate_squared +=
            angular_momentum_rate_vector[axis]*angular_momentum_rate_vector[axis];
        angular_momentum_dot_rate += component*angular_momentum_rate_vector[axis];
    }

    double angular_momentum_rate;
    double reduced_angular_momentum_rate;
    if (angular_momentum_squared > 0.0) {
        angular_momentum_rate =
            angular_momentum_dot_rate/sqrt(angular_momentum_squared);
        reduced_angular_momentum_rate =
            angular_momentum_rate/(reduced_mass*total_mass);
    }
    else {
        /* |L| has only a forward derivative when a torque creates L from an exactly radial state. */
        angular_momentum_rate = sqrt(angular_momentum_rate_squared);
        reduced_angular_momentum_rate = 0.0;
    }

    const double x = orbit_state->reduced_separation;
    const double j = orbit_state->reduced_angular_momentum;
    const double nu = orbit_state->symmetric_mass_ratio;
    const double dh_dx = pn_binary_reduced_hamiltonian_dx(
        x, pr_hat, j, nu, orbit_state->use_1pn, orbit_state->use_2pn);
    const double dh_dpr = pn_binary_reduced_hamiltonian_dpr(
        x, pr_hat, j, nu, orbit_state->use_1pn, orbit_state->use_2pn);
    const double dh_dj = pn_binary_reduced_hamiltonian_dj(
        x, pr_hat, j, nu, orbit_state->use_1pn, orbit_state->use_2pn);
    const double energy_rate = reduced_mass
        * (dh_dx*radial_velocity/total_mass
           + dh_dpr*pr_hat_rate
           + dh_dj*reduced_angular_momentum_rate);

    if (!isfinite(energy_rate) || !isfinite(angular_momentum_rate))
        return 0;

    rate_state->energy_rate = energy_rate;
    rate_state->angular_momentum_rate = angular_momentum_rate;
    rate_state->valid = 1;
    return 1;
}


/**
 * @brief Prepare canonical Newtonian vector diagnostics from an isolated pair state.
 *
 * This is separate from compute_pn_orbit_state so the LRL and eccentricity-vector cross product
 * is evaluated only when a vector output was requested. All requested vectors then share it.
 */
int compute_pn_orbit_vector_state(const PNOrbitState *orbit_state,
    PNOrbitVectorState *vector_state)
{
    if (vector_state == NULL)
        return 0;
    *vector_state = (PNOrbitVectorState){0};
    if (orbit_state == NULL || !orbit_state->valid)
        return 0;

    const double reduced_mass =
        orbit_state->total_mass*orbit_state->symmetric_mass_ratio;
    const double lrl_normalization =
        reduced_mass*reduced_mass*orbit_state->total_mass;
    const double separation = orbit_state->reduced_separation*orbit_state->total_mass;
    if (!(lrl_normalization > 0.0) || !(separation > 0.0)
        || !isfinite(lrl_normalization) || !isfinite(separation))
        return 0;

    const double *p = orbit_state->relative_momentum;
    const double *angular_momentum = orbit_state->angular_momentum_vector;
    double momentum_cross_angular_momentum[3];
    momentum_cross_angular_momentum[0] =
        p[1]*angular_momentum[2] - p[2]*angular_momentum[1];
    momentum_cross_angular_momentum[1] =
        p[2]*angular_momentum[0] - p[0]*angular_momentum[2];
    momentum_cross_angular_momentum[2] =
        p[0]*angular_momentum[1] - p[1]*angular_momentum[0];

    for (int axis = 0; axis < 3; ++axis) {
        vector_state->angular_momentum_vector[axis] = angular_momentum[axis];
        vector_state->laplace_runge_lenz_vector[axis] =
            momentum_cross_angular_momentum[axis]
            - lrl_normalization*orbit_state->relative_position[axis]/separation;
        vector_state->eccentricity_vector[axis] =
            vector_state->laplace_runge_lenz_vector[axis]/lrl_normalization;
        if (!isfinite(vector_state->angular_momentum_vector[axis])
            || !isfinite(vector_state->laplace_runge_lenz_vector[axis])
            || !isfinite(vector_state->eccentricity_vector[axis]))
            return 0;
    }

    vector_state->valid = 1;
    return 1;
}


static int copy_binary_vector(const PNOrbitVectorState *state, const double source[3],
    double vector[3])
{
    if (vector == NULL)
        return 0;
    if (state == NULL || !state->valid || source == NULL) {
        for (int axis = 0; axis < 3; ++axis)
            vector[axis] = NAN;
        return 0;
    }

    for (int axis = 0; axis < 3; ++axis)
        vector[axis] = source[axis];
    return 1;
}


// ------------------------------------------------------------------------------------------------
// Total energy and angular momentum
// ------------------------------------------------------------------------------------------------

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


double orbital_angular_momentum(const SystemAngularMomentumState *state)
{
    return state != NULL && state->valid ? state->orbital_magnitude : NAN;
}


double orbital_angular_momentum_com(const SystemAngularMomentumState *state)
{
    return state != NULL && state->valid ? state->orbital_com_magnitude : NAN;
}


double spin_angular_momentum(const SystemAngularMomentumState *state)
{
    return state != NULL && state->valid ? state->spin_magnitude : NAN;
}


double total_angular_momentum(const SystemAngularMomentumState *state)
{
    return state != NULL && state->valid ? state->total_magnitude : NAN;
}


double total_angular_momentum_com(const SystemAngularMomentumState *state)
{
    return state != NULL && state->valid ? state->total_com_magnitude : NAN;
}


static int copy_system_angular_momentum_vector(
    const SystemAngularMomentumState *state, const double source[3], double vector[3])
{
    if (vector == NULL)
        return 0;
    if (state == NULL || !state->valid || source == NULL) {
        for (int axis = 0; axis < 3; ++axis)
            vector[axis] = NAN;
        return 0;
    }

    for (int axis = 0; axis < 3; ++axis)
        vector[axis] = source[axis];
    return 1;
}


int orbital_angular_momentum_vector(
    const SystemAngularMomentumState *state, double vector[3])
{
    return copy_system_angular_momentum_vector(
        state, state != NULL ? state->orbital_vector : NULL, vector);
}


int orbital_angular_momentum_vector_com(
    const SystemAngularMomentumState *state, double vector[3])
{
    return copy_system_angular_momentum_vector(
        state, state != NULL ? state->orbital_com_vector : NULL, vector);
}


int spin_angular_momentum_vector(
    const SystemAngularMomentumState *state, double vector[3])
{
    return copy_system_angular_momentum_vector(
        state, state != NULL ? state->spin_vector : NULL, vector);
}


int total_angular_momentum_vector(
    const SystemAngularMomentumState *state, double vector[3])
{
    return copy_system_angular_momentum_vector(
        state, state != NULL ? state->total_vector : NULL, vector);
}


int total_angular_momentum_vector_com(
    const SystemAngularMomentumState *state, double vector[3])
{
    return copy_system_angular_momentum_vector(
        state, state != NULL ? state->total_com_vector : NULL, vector);
}


// ------------------------------------------------------------------------------------------------
// Direct pair kinematics
// ------------------------------------------------------------------------------------------------

// Return the current coordinate separation of the selected pair.
double separation_binary(const NewtonianRelativeOrbitState *state)
{
    return state != NULL && state->valid ? state->separation : NAN;
}


// Return dr/dt from the full coordinate velocities of the selected pair.
double radial_velocity_binary(const NewtonianRelativeOrbitState *state)
{
    if (state == NULL || !state->valid || !state->velocity_valid
        || !(state->separation > 0.0))
        return NAN;
    return state->radial_velocity_product/state->separation;
}


// Return the instantaneous orbital frequency |r x v_rel|/r^2 from the full coordinate velocity.
double orbital_frequency_instantaneous_binary(const NewtonianRelativeOrbitState *state)
{
    if (state == NULL || !state->valid || !state->velocity_valid
        || !(state->separation > 0.0))
        return NAN;

    const double *r = state->relative_position;
    const double *v = state->relative_velocity;
    const double cross[3] = {
        r[1]*v[2] - r[2]*v[1],
        r[2]*v[0] - r[0]*v[2],
        r[0]*v[1] - r[1]*v[0]
    };
    const double cross_squared =
        cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2];
    const double frequency = sqrt(cross_squared)/(state->separation*state->separation);
    return isfinite(frequency) ? frequency : NAN;
}


// ------------------------------------------------------------------------------------------------
// Radial turning-point geometry
// ------------------------------------------------------------------------------------------------

double pericenter_radial_pn(const PNOrbitState *state)
{
    return state != NULL && state->radial_turning_points_valid ? state->pericenter : NAN;
}


double apocenter_radial_pn(const PNOrbitState *state)
{
    return state != NULL && state->radial_turning_points_valid ? state->apocenter : NAN;
}


double semilatus_rectum_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->radial_turning_points_valid)
        return NAN;

    if (!(state->pericenter > 0.0) || !(state->apocenter > 0.0)
        || !isfinite(state->pericenter) || !isfinite(state->apocenter))
        return NAN;

    const double semilatus_rectum =
        2.0*state->pericenter/(1.0 + state->pericenter/state->apocenter);
    return isfinite(semilatus_rectum) && semilatus_rectum > 0.0
        ? semilatus_rectum : NAN;
}


// ------------------------------------------------------------------------------------------------
// Semimajor-axis definitions
// ------------------------------------------------------------------------------------------------

// Return the Newtonian osculating semimajor axis a = -m/(2E) from the specific orbital energy.
double semimajor_axis_newtonian(const NewtonianRelativeOrbitState *state)
{
    if (state == NULL || !state->valid || !state->velocity_valid)
        return NAN;

    const double specific_energy =
        0.5*state->relative_speed_squared - state->total_mass/state->separation;
    if (specific_energy == 0.0)
        return INFINITY;

    const double semimajor_axis = -state->total_mass/(2.0*specific_energy);
    return isnan(semimajor_axis) ? NAN : semimajor_axis;
}


// Return the semimajor axis a_r = (r_a + r_p) / 2 based on the apocenter and pericenter distances.
double semimajor_axis_radial_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->radial_turning_points_valid)
        return NAN;
    return 0.5 * (state->pericenter + state->apocenter);
}


/**
 * @brief Return the PN-expanded ADM radial semimajor axis through the enabled PN order.
 *
 * This is Eq. (20a) of Memmesheimer, Gopakumar & Schaefer (2004), truncated after 2PN and
 * restored to the code's physical length units by multiplying the dimensionless result by M.
 */
double semimajor_axis_radial_analytic_pn(const PNOrbitState *state)
{
    long double energy_scale;
    long double energy_angular_momentum_squared;
    long double nu;
    if (!pn_orbit_bound_invariants(
            state, &energy_scale, &energy_angular_momentum_squared, &nu))
        return NAN;

    long double expansion = 1.0L;
    if (state->use_1pn)
        expansion += energy_scale / 4.0L * (-7.0L + nu);
    if (state->use_2pn)
        expansion += energy_scale*energy_scale / 16.0L
            * (1.0L + 10.0L*nu + nu*nu
               + (-68.0L + 44.0L*nu)/energy_angular_momentum_squared);

    const long double semimajor_axis = state->total_mass * expansion / energy_scale;
    return isfinite(semimajor_axis) && semimajor_axis > 0.0L
        ? (double)semimajor_axis : NAN;
}


// Return the energy-defined semimajor axis a_E = -M/(2 E_hat) with the PN energy.
double semimajor_axis_energy_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->valid || !isfinite(state->reduced_energy)
        || !(state->total_mass > 0.0) || !isfinite(state->total_mass))
        return NAN;
    if (state->reduced_energy == 0.0)
        return INFINITY;

    const double semimajor_axis = -state->total_mass / (2.0*state->reduced_energy);
    return isnan(semimajor_axis) ? NAN : semimajor_axis;
}


// Return a_n = (M/Omega_r^2)^(1/3) (applying Kepler's law to the PN radial frequency).
double semimajor_axis_mean_motion_pn(const PNOrbitFrequencyState *state)
{
    if (state == NULL || !state->valid || !(state->radial_frequency > 0.0))
        return NAN;

    const long double semimajor_axis = cbrtl(
        state->total_mass/(state->radial_frequency*state->radial_frequency));
    return isfinite(semimajor_axis) && semimajor_axis > 0.0L
        ? (double)semimajor_axis : NAN;
}


// Return a_phi = (M/Omega_phi^2)^(1/3) (applying Kepler's law to the PN azimuthal frequency).
double semimajor_axis_azimuthal_pn(const PNOrbitFrequencyState *state)
{
    if (state == NULL || !state->valid || !(state->azimuthal_frequency > 0.0))
        return NAN;

    const long double semimajor_axis = cbrtl(
        state->total_mass/(state->azimuthal_frequency*state->azimuthal_frequency));
    return isfinite(semimajor_axis) && semimajor_axis > 0.0L
        ? (double)semimajor_axis : NAN;
}


// ------------------------------------------------------------------------------------------------
// Eccentricity definitions
// ------------------------------------------------------------------------------------------------

// Return the magnitude of the Newtonian osculating eccentricity vector.
double eccentricity_newtonian(const NewtonianRelativeOrbitState *state)
{
    if (state == NULL || !state->valid || !state->velocity_valid)
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


// Return the eccentricity e_r=(r_a-r_p)/(r_a+r_p) based on the apocenter and pericenter distances.
double eccentricity_radial_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->radial_turning_points_valid)
        return NAN;

    const double radius_sum = state->pericenter + state->apocenter;
    if (!(radius_sum > 0.0) || !isfinite(radius_sum))
        return NAN;

    const double eccentricity = (state->apocenter - state->pericenter) / radius_sum;
    if (eccentricity < 0.0 && eccentricity > -64.0 * DBL_EPSILON)
        return 0.0;
    return eccentricity;
}


/**
 * @brief Return the PN-expanded ADM radial eccentricity through the enabled PN order.
 *
 * This is Eq. (20b) of Memmesheimer, Gopakumar & Schaefer (2004), truncated after 2PN. A negative
 * value of the truncated e_r^2 series produces NAN unless it is consistent with floating-point
 * roundoff.
 */
double eccentricity_radial_analytic_pn(const PNOrbitState *state)
{
    long double energy_scale;
    long double energy_angular_momentum_squared;
    long double nu;
    if (!pn_orbit_bound_invariants(
            state, &energy_scale, &energy_angular_momentum_squared, &nu))
        return NAN;

    long double eccentricity_squared = 1.0L - energy_angular_momentum_squared;
    long double scale = 1.0L + fabsl(energy_angular_momentum_squared);

    if (state->use_1pn) {
        const long double correction = energy_scale / 4.0L
            * (24.0L - 4.0L*nu
               + 5.0L*(-3.0L + nu)*energy_angular_momentum_squared);
        eccentricity_squared += correction;
        scale += fabsl(correction);
    }

    if (state->use_2pn) {
        const long double correction = energy_scale*energy_scale / 8.0L
            * (52.0L + 2.0L*nu + 2.0L*nu*nu
               - (80.0L - 55.0L*nu + 4.0L*nu*nu)
                 * energy_angular_momentum_squared
               - 8.0L*(-17.0L + 11.0L*nu)/energy_angular_momentum_squared);
        eccentricity_squared += correction;
        scale += fabsl(correction);
    }

    const long double roundoff_tolerance = 256.0L * DBL_EPSILON * fmaxl(1.0L, scale);
    if (eccentricity_squared < 0.0L) {
        if (eccentricity_squared >= -roundoff_tolerance)
            return 0.0;
        return NAN;
    }

    return isfinite(eccentricity_squared) ? (double)sqrtl(eccentricity_squared) : NAN;
}


/**
 * @brief Return the ADM quasi-Keplerian time eccentricity through the enabled PN order.
 *
 * This uses Eq. (21a) of Memmesheimer, Gopakumar & Schaefer (2004), truncated after 2PN and
 * evaluated in units G = c = 1. The turning-point-defined e_r is the common baseline.
 */
double eccentricity_time_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->radial_turning_points_valid)
        return NAN;

    const double radial_eccentricity = eccentricity_radial_pn(state);
    long double energy_scale;
    long double energy_angular_momentum_squared;
    long double nu;
    if (!isfinite(radial_eccentricity)
        || !pn_orbit_bound_invariants(
            state, &energy_scale, &energy_angular_momentum_squared, &nu))
        return NAN;

    long double ratio = 1.0L;

    if (state->use_1pn)
        ratio += energy_scale / 2.0L * (-8.0L + 3.0L*nu);

    if (state->use_2pn) {
        ratio += energy_scale*energy_scale
            / (8.0L*energy_angular_momentum_squared)
            * (-34.0L + 22.0L*nu
               + (-60.0L + 24.0L*nu)*sqrtl(energy_angular_momentum_squared)
               + (72.0L - 33.0L*nu + 12.0L*nu*nu)
                 * energy_angular_momentum_squared);
    }

    const long double eccentricity = radial_eccentricity * ratio;
    return isfinite(eccentricity) && eccentricity >= 0.0L ? (double)eccentricity : NAN;
}


/**
 * @brief Return the ADM quasi-Keplerian angular eccentricity through the enabled PN order.
 *
 * This uses Eq. (21b) of Memmesheimer, Gopakumar & Schaefer (2004), truncated after 2PN and
 * evaluated in units G = c = 1. As for e_t, the turning-point-defined e_r is the common baseline.
 */
double eccentricity_angular_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->radial_turning_points_valid)
        return NAN;

    const double radial_eccentricity = eccentricity_radial_pn(state);
    long double energy_scale;
    long double energy_angular_momentum_squared;
    long double nu;
    if (!isfinite(radial_eccentricity)
        || !pn_orbit_bound_invariants(
            state, &energy_scale, &energy_angular_momentum_squared, &nu))
        return NAN;

    long double ratio = 1.0L;

    if (state->use_1pn)
        ratio += energy_scale / 2.0L * nu;

    if (state->use_2pn)
        ratio += energy_scale*energy_scale
            / (32.0L*energy_angular_momentum_squared)
            * (136.0L - 56.0L*nu - 15.0L*nu*nu
               + nu*(20.0L + 11.0L*nu)*energy_angular_momentum_squared);

    const long double eccentricity = radial_eccentricity * ratio;
    return isfinite(eccentricity) && eccentricity >= 0.0L ? (double)eccentricity : NAN;
}


/**
 * @brief Return the eccentricity inferred from azimuthal frequencies at the radial turning points.
 *
 * Omega_p and Omega_a are obtained from Hamilton's equation Omega = (1/M)*dh/dj at p_r = 0.
 * Their common 1/M factor cancels from e_Omega. The logarithmic form is algebraically equivalent
 * to the square-root definition but avoids directly subtracting nearly equal square roots.
 */
double eccentricity_frequency_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->radial_turning_points_valid
        || !(state->total_mass > 0.0) || !(state->reduced_angular_momentum > 0.0))
        return NAN;

    const double x_pericenter = state->pericenter/state->total_mass;
    const double x_apocenter = state->apocenter/state->total_mass;
    const double omega_pericenter = pn_binary_reduced_hamiltonian_dj(
        x_pericenter, 0.0, state->reduced_angular_momentum,
        state->symmetric_mass_ratio, state->use_1pn, state->use_2pn);
    const double omega_apocenter = pn_binary_reduced_hamiltonian_dj(
        x_apocenter, 0.0, state->reduced_angular_momentum,
        state->symmetric_mass_ratio, state->use_1pn, state->use_2pn);
    if (!(omega_pericenter > 0.0) || !(omega_apocenter > 0.0))
        return NAN;

    long double frequency_difference =
        (long double)omega_pericenter - (long double)omega_apocenter;
    const long double frequency_scale =
        fmaxl((long double)omega_pericenter, (long double)omega_apocenter);
    if (frequency_difference < 0.0L) {
        if (frequency_difference >= -64.0L*DBL_EPSILON*frequency_scale)
            frequency_difference = 0.0L;
        else
            return NAN;
    }

    const long double log_frequency_ratio = log1pl(
        frequency_difference/(long double)omega_apocenter);
    const long double eccentricity = tanhl(0.25L*log_frequency_ratio);
    return isfinite(eccentricity) && eccentricity >= 0.0L
        ? (double)eccentricity : NAN;
}


// ------------------------------------------------------------------------------------------------
// Characteristic gravitational-wave frequencies
// ------------------------------------------------------------------------------------------------

// Characteristic peak gravitational-wave frequency of an eccentric binary according to Wen 2003.
static double peak_gravitational_wave_frequency_wen(
    double total_mass, double semilatus_rectum, double eccentricity)
{
    if (!(total_mass > 0.0) || !isfinite(total_mass)
        || !(semilatus_rectum > 0.0) || !isfinite(semilatus_rectum)
        || eccentricity < 0.0 || !(eccentricity < 1.0) || !isfinite(eccentricity))
        return NAN;

    const double frequency = pow(1.0 + eccentricity, 1.1954)
        * sqrt(total_mass/semilatus_rectum)
        / (acos(-1.0)*semilatus_rectum);
    return isfinite(frequency) ? frequency : NAN;
}


// Return Wen's peak-frequency fit evaluated with the Newtonian osculating a and e.
double peak_gravitational_wave_frequency_wen_newtonian(
    const NewtonianRelativeOrbitState *state)
{
    const double semimajor_axis = semimajor_axis_newtonian(state);
    const double eccentricity = eccentricity_newtonian(state);
    if (!(semimajor_axis > 0.0) || !isfinite(semimajor_axis)
        || eccentricity < 0.0 || !(eccentricity < 1.0) || !isfinite(eccentricity))
        return NAN;

    const double one_minus_eccentricity_squared = fma(-eccentricity, eccentricity, 1.0);
    const double semilatus_rectum = semimajor_axis*one_minus_eccentricity_squared;
    return peak_gravitational_wave_frequency_wen(
        state->total_mass, semilatus_rectum, eccentricity);
}


// Return Wen's peak-frequency fit evaluated with the numerical PN radial elements a_r and e_r.
double peak_gravitational_wave_frequency_wen_radial_pn(const PNOrbitState *state)
{
    if (state == NULL)
        return NAN;

    const double semilatus_rectum = semilatus_rectum_pn(state);
    const double eccentricity = eccentricity_radial_pn(state);
    return peak_gravitational_wave_frequency_wen(
        state->total_mass, semilatus_rectum, eccentricity);
}


// ------------------------------------------------------------------------------------------------
// Pair energy
// ------------------------------------------------------------------------------------------------

/**
 * @brief Return the conservative isolated-pair binding energy in the code's energy units.
 *
 * PNOrbitState stores the reduced Hamiltonian E_hat = H_pair / mu. Multiplication by the pair's
 * reduced mass mu = M*nu restores the physical pair Hamiltonian. Rest-mass energy, pair COM motion,
 * other bodies and dissipative terms are excluded.
 */
double energy_binary_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->valid)
        return NAN;

    const double reduced_mass = state->total_mass * state->symmetric_mass_ratio;
    const double energy = reduced_mass * state->reduced_energy;
    return isfinite(energy) ? energy : NAN;
}


double energy_rate_binary_pn(const PNOrbitRateState *state)
{
    return state != NULL && state->valid ? state->energy_rate : NAN;
}


// ------------------------------------------------------------------------------------------------
// Pair angular momentum
// ------------------------------------------------------------------------------------------------

/**
 * @brief Return the magnitude of the isolated pair's canonical orbital angular momentum.
 *
 * PNOrbitState stores j = J/(mu*M), where J = |r x p_rel|. Multiplication by mu*M restores
 * the physical angular momentum in the code's canonical units. The pair COM motion is excluded.
 */
double angular_momentum_binary_pn(const PNOrbitState *state)
{
    if (state == NULL || !state->valid)
        return NAN;

    const double reduced_mass = state->total_mass * state->symmetric_mass_ratio;
    const double angular_momentum =
        reduced_mass * state->total_mass * state->reduced_angular_momentum;
    return isfinite(angular_momentum) && angular_momentum >= 0.0
        ? angular_momentum : NAN;
}


double angular_momentum_rate_binary_pn(const PNOrbitRateState *state)
{
    return state != NULL && state->valid ? state->angular_momentum_rate : NAN;
}


// Return the vector L = r x p_rel for the isolated pair.
int angular_momentum_vector_binary_pn(const PNOrbitVectorState *state, double vector[3])
{
    return copy_binary_vector(
        state, state != NULL ? state->angular_momentum_vector : NULL, vector);
}


// ------------------------------------------------------------------------------------------------
// Other canonical pair vectors
// ------------------------------------------------------------------------------------------------

// Return the Newtonian canonical Laplace-Runge-Lenz vector of the isolated pair.
int laplace_runge_lenz_vector_binary_pn(const PNOrbitVectorState *state, double vector[3])
{
    return copy_binary_vector(
        state, state != NULL ? state->laplace_runge_lenz_vector : NULL, vector);
}


// Return the canonical osculating eccentricity vector A/(mu^2*M).
int eccentricity_vector_binary_pn(const PNOrbitVectorState *state, double vector[3])
{
    return copy_binary_vector(state, state != NULL ? state->eccentricity_vector : NULL, vector);
}


// ------------------------------------------------------------------------------------------------
// Orbit-averaged frequencies, periods and precession
// ------------------------------------------------------------------------------------------------

double frequency_radial_pn(const PNOrbitFrequencyState *state)
{
    return state != NULL && state->valid ? state->radial_frequency : NAN;
}


double frequency_azimuthal_pn(const PNOrbitFrequencyState *state)
{
    return state != NULL && state->valid ? state->azimuthal_frequency : NAN;
}


double period_radial_pn(const PNOrbitFrequencyState *state)
{
    if (state == NULL || !state->valid || !(state->radial_frequency > 0.0))
        return NAN;
    return 2.0*acos(-1.0)/state->radial_frequency;
}


double period_azimuthal_pn(const PNOrbitFrequencyState *state)
{
    if (state == NULL || !state->valid || !(state->azimuthal_frequency > 0.0))
        return NAN;
    return 2.0*acos(-1.0)/state->azimuthal_frequency;
}


double periastron_advance_factor_pn(const PNOrbitFrequencyState *state)
{
    return state != NULL && state->valid ? state->periastron_advance_factor : NAN;
}


double pericenter_advance_pn(const PNOrbitFrequencyState *state)
{
    if (state == NULL || !state->valid)
        return NAN;
    return 2.0*acos(-1.0)*(state->periastron_advance_factor - 1.0);
}


double precession_rate_pn(const PNOrbitFrequencyState *state)
{
    if (state == NULL || !state->valid)
        return NAN;
    return (state->periastron_advance_factor - 1.0)*state->radial_frequency;
}
