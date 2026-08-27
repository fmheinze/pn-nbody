/**
 * @file merger.c
 * @brief Routines for merging bodies in the simulations.
 *
 * Routines for merging bodies in the simulations, such as finding merger pairs, updating the
 * ode_params and the state vector after merger, as well as prescriptions for the properties
 * of the merger remnant.
 */

#include "eom.h"
#include "output.h"
#include "stdio.h"
#include "utils.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define C_KMS 299792.458

// Lousto-Zlochower recoil coefficients
#define LZ_KICK_A_KMS   1.2e4
#define LZ_KICK_B      -0.93
#define LZ_KICK_H_KMS   6.9e3
#define LZ_KICK_XI_RAD  (145.0 * M_PI / 180.0)

#define LZ_KICK_V11_KMS 3677.76
#define LZ_KICK_VA_KMS  2481.21
#define LZ_KICK_VB_KMS  1792.45
#define LZ_KICK_VC_KMS  1506.52

// Lousto-Zlochower remnant mass/spin coefficients
#define LZ_E2       0.341
#define LZ_E3       0.522
#define LZ_ES       0.673
#define LZ_EDELTA  -0.36
#define LZ_EA      -0.014
#define LZ_EB       0.045
#define LZ_EC       0.0
#define LZ_ED       0.2611
#define LZ_EE       0.0959
#define LZ_EF       0.0

#define LZ_J2      -2.81
#define LZ_J3       1.69
#define LZ_JA      -2.97
#define LZ_JB      -1.73


// ------------------------------------------------------------------------------------------------
// Structs
// ------------------------------------------------------------------------------------------------

// Binary struct for the Lousto-Zlochower remnant prescriptions.
typedef struct {
    int idx1;          // Smaller BH
    int idx2;          // Larger BH

    double m1;
    double m2;
    double M;
    double q;
    double eta;

    double Lhat[3];

    double alpha1[3];  // Dimensionless spin of smaller BH
    double alpha2[3];  // Dimensionless spin of larger BH

    double alpha1_par;
    double alpha2_par;

    double alpha1_perp[3];
    double alpha2_perp[3];

    double e1[3];      // In-plane separation direction
    double e2[3];      // In-plane direction rotated by 90 deg
} lz_binary_data;


// Remnant parameter struct for the Lousto-Zlochower remnant prescriptions.
typedef struct {
    double mass;
    double erad_frac;
    double chi[3];
    double v_kick[3];
    double v_kick_kms[3];
} lz_remnant;


// ------------------------------------------------------------------------------------------------
// Helper functions
// ------------------------------------------------------------------------------------------------

// Initializes the remnant parameter struct with a remnant mass and else zeros.
static void init_simple_remnant(lz_remnant *rem, double mass)
{
    rem->mass = mass;
    rem->erad_frac = 0.0;

    for (int k = 0; k < 3; k++) {
        rem->chi[k] = 0.0;
        rem->v_kick[k] = 0.0;
        rem->v_kick_kms[k] = 0.0;
    }
}

/**
 * @brief Computes the orbital angular momentum of a pair of bodies.
 *
 * @param[in]   params       Parameter struct containing general information about the system
 * @param[in]   w            State vector w = [positions, momenta, spins]
 * @param[in]   i            Index of the first body
 * @param[in]   j            Index of the second body
 * @param[out]  L            Orbital angular momentum vector
 */
void compute_pair_orbital_angular_momentum(struct ode_params *params, double *w, int i, int j,
    double L[3])
{
    const int num_dim = params->num_dim;
    const int num_bodies = params->num_bodies_initial;
    const int array_half = num_dim * num_bodies;

    if (num_dim != 3)
        errorexit("Orbital angular momentum vector requires num_dim = 3");

    const double mi = params->masses[i];
    const double mj = params->masses[j];
    const double mu = mi * mj / (mi + mj);

    double r[3];
    double p_rel[3];

    for (int k = 0; k < 3; k++) {
        const double xi = w[num_dim * i + k];
        const double xj = w[num_dim * j + k];

        const double pi = w[array_half + num_dim * i + k];
        const double pj = w[array_half + num_dim * j + k];

        r[k] = xi - xj;

        p_rel[k] = mu * (pi / mi - pj / mj);
    }

    L[0] = r[1] * p_rel[2] - r[2] * p_rel[1];
    L[1] = r[2] * p_rel[0] - r[0] * p_rel[2];
    L[2] = r[0] * p_rel[1] - r[1] * p_rel[0];
}


/**
 * @brief Computes the Newtonian binding energy of a pair of bodies.
 *
 * @param[in]   params       Parameter struct containing general information about the system
 * @param[in]   w            State vector w = [positions, momenta, spins]
 * @param[in]   i            Index of the first body
 * @param[in]   j            Index of the second body
 * @returns Newtonian binding energy for the bodies `i` and `j`
 */
double pair_newtonian_binding_energy(struct ode_params *params, double *w, int i, int j)
{
    const int num_dim = params->num_dim;
    const int num_bodies = params->num_bodies_initial;
    const int array_half = num_dim * num_bodies;

    if (!params->active[i] || !params->active[j] || i == j)
        return DBL_MAX;

    const double mi = params->masses[i];
    const double mj = params->masses[j];
    const double M  = mi + mj;

    if (!(mi > 0.0) || !(mj > 0.0) || !(M > 0.0))
        return DBL_MAX;

    double r2 = 0.0;
    double vrel2 = 0.0;

    for (int k = 0; k < num_dim; k++) {
        const double xi = w[num_dim * i + k];
        const double xj = w[num_dim * j + k];

        const double pi = w[array_half + num_dim * i + k];
        const double pj = w[array_half + num_dim * j + k];

        const double dx = xi - xj;
        const double dv = pi / mi - pj / mj;

        r2 += dx * dx;
        vrel2 += dv * dv;
    }

    if (!(r2 > 0.0) || !isfinite(r2) || !isfinite(vrel2))
        return DBL_MAX;

    const double r = sqrt(r2);
    const double mu = mi * mj / M;

    return 0.5 * mu * vrel2 - mi * mj / r;
}


/**
 * @brief Evaluates whether a pair of bodies is gravitationally bound.
 *
 * To test whether the pair is bound the function checks whether the Newtonian binding energy is
 * negative. In relativistic settings this is only an approximation.
 *
 * @param[in]   params       Parameter struct containing general information about the system
 * @param[in]   w            State vector w = [positions, momenta, spins]
 * @param[in]   i            Index of the first body
 * @param[in]   j            Index of the second body
 * @returns 1 if the two bodies are gravitationally bound, else 0
 */
int pair_is_bound(struct ode_params *params, double *w, int i, int j)
{
    const double E = pair_newtonian_binding_energy(params, w, i, j);
    if (!isfinite(E)) return 0;
    return E < 0.0;
}


// ------------------------------------------------------------------------------------------------
// Barausse remnant mass (see Barausse et al. 2012)
// ------------------------------------------------------------------------------------------------

// Innermost stable circular orbit for the remnant mass description of Barausse et al. 2012
double barausse_E_ISCO(double a)
{
    if (a >  1.0) a =  1.0;
    if (a < -1.0) a = -1.0;

    const double Z1 = 1.0 + pow(1.0-a*a, 1.0/3.0) * (pow(1.0+a, 1.0/3.0) + pow(1.0-a, 1.0/3.0));
    const double Z2 = sqrt(3.0 * a*a + Z1 * Z1);

    const double r_eq_ISCO = 3.0 + Z2 - sign_double(a) * sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));

    return sqrt(1.0 - 2.0 / (3.0 * r_eq_ISCO));
}


/**
 * @brief Computes NR-informed remnant mass for the merger of bodies i and j.
 *
 * Based on Barausse et al. 2012.
 *
 * @param[in]   params     Parameter struct containing general information about the system
 * @param[in]   w          State vector w = [positions, momenta, spins]
 * @param[in]   i          Index of the first merger object
 * @param[in]   j          Index of the second merger object
 * @return Remnant mass for the merger of bodies `i` and `j`
 */
double barausse_remnant_mass(struct ode_params *ode_params, double *w, int i, int j)
{
    const int num_dim = ode_params->num_dim;
    const int array_half = num_dim * ode_params->num_bodies_initial;
    const int spin_offset = 2 * array_half;

    if (num_dim != 3)
        errorexit("Remnant-mass fit requires num_dim = 3");

    const double m1 = ode_params->masses[i];
    const double m2 = ode_params->masses[j];
    const double M = m1 + m2;

    if (M <= 0.0)
        return 0.0;

    const double eta = (m1 * m2) / (M * M);

    double Lhat[3];
    compute_pair_orbital_angular_momentum(ode_params, w, i, j, Lhat);

    const double Lnorm = norm(Lhat, num_dim);

    if (Lnorm <= 0.0) {
        // Degenerate case: nonspinning projection
        Lhat[0] = 0.0;
        Lhat[1] = 0.0;
        Lhat[2] = 0.0;
    } else {
        for (int k = 0; k < 3; k++)
            Lhat[k] /= Lnorm;
    }

    const double chi1_dot_L =
        w[spin_offset + 3*i + 0] * Lhat[0]
      + w[spin_offset + 3*i + 1] * Lhat[1]
      + w[spin_offset + 3*i + 2] * Lhat[2];

    const double chi2_dot_L =
        w[spin_offset + 3*j + 0] * Lhat[0]
      + w[spin_offset + 3*j + 1] * Lhat[1]
      + w[spin_offset + 3*j + 2] * Lhat[2];

    double a_tilde = (m1*m1 * chi1_dot_L + m2*m2 * chi2_dot_L) / (M * M);

    if (a_tilde >  1.0) a_tilde =  1.0;
    if (a_tilde < -1.0) a_tilde = -1.0;

    const double p0 = 0.04827;
    const double p1 = 0.01707;

    const double e_isco = barausse_E_ISCO(a_tilde);

    const double e_rad_frac = (1.0 - e_isco) * eta + 4.0 * eta * eta * (4.0 * p0
        + 16.0 * p1 * a_tilde * (a_tilde + 1.0) + e_isco - 1.0);

    return M * (1.0 - e_rad_frac);
}


// ------------------------------------------------------------------------------------------------
// Lousto - Zlochower remnant fits (see Lousto et al. 2010/2012 or Blecha et al. 2016)
// ------------------------------------------------------------------------------------------------

/**
 * @brief Computes the binary data necessary for the Lousto - Zlochower remnant prescription.
 *
 * The function fills the binary data struct from information in the ode_params and state vector w.
 *
 * @param[in]   params     Parameter struct containing general information about the system
 * @param[in]   w          State vector w = [positions, momenta, spins]
 * @param[in]   i          Index of the first merger object
 * @param[in]   j          Index of the second merger object
 * @param[out]  bd         Merging binary data struct
 * @return 1 if the binary data could be computed, else 0
 */
int get_lz_binary_data(struct ode_params *params, double *w, int i, int j, lz_binary_data *bd)
{
    const int num_dim = params->num_dim;
    const int num_bodies = params->num_bodies_initial;
    const int array_half = num_dim * num_bodies;
    const int spin_offset = 2 * array_half;

    if (num_dim != 3)
        errorexit("Lousto/Zlochower remnant prescription requires num_dim = 3");

    if (!params->active[i] || !params->active[j] || i == j)
        return 0;

    int idx1 = i;
    int idx2 = j;

    // Relabel so idx1 is the smaller BH and q <= 1
    if (params->masses[idx1] > params->masses[idx2]) {
        idx1 = j;
        idx2 = i;
    }

    const double m1 = params->masses[idx1];
    const double m2 = params->masses[idx2];
    const double M = m1 + m2;

    if (!(m1 > 0.0) || !(m2 > 0.0) || !(M > 0.0))
        return 0;

    bd->idx1 = idx1;
    bd->idx2 = idx2;

    bd->m1 = m1;
    bd->m2 = m2;
    bd->M = M;
    bd->q = m1 / m2;
    bd->eta = bd->q / ((1.0 + bd->q) * (1.0 + bd->q));

    double x1[3], x2[3];
    double p1[3], p2[3];

    for (int k = 0; k < 3; k++) {
        x1[k] = w[num_dim * idx1 + k];
        x2[k] = w[num_dim * idx2 + k];

        p1[k] = w[array_half + num_dim * idx1 + k];
        p2[k] = w[array_half + num_dim * idx2 + k];

        bd->alpha1[k] = w[spin_offset + 3 * idx1 + k];
        bd->alpha2[k] = w[spin_offset + 3 * idx2 + k];
    }

    // e1: separation direction from larger BH to smaller BH
    double e1_raw[3];

    for (int k = 0; k < 3; k++)
        e1_raw[k] = x1[k] - x2[k];

    if (!safe_normalize(e1_raw, bd->e1))
        return 0;

    // Relative orbital angular momentum
    const double mu = m1 * m2 / M;

    double r12[3];
    double p_rel[3];

    for (int k = 0; k < 3; k++) {
        r12[k] = x1[k] - x2[k];
        p_rel[k] = mu * (p1[k] / m1 - p2[k] / m2);
    }

    double L_raw[3];
    cross_product(r12, p_rel, L_raw);
    if (!safe_normalize(L_raw, bd->Lhat))
        return 0;

    // e2 = Lhat x e1
    double e2_raw[3];
    cross_product(bd->Lhat, bd->e1, e2_raw);
    if (!safe_normalize(e2_raw, bd->e2))
        return 0;

    bd->alpha1_par = dot_product(bd->alpha1, bd->Lhat, 3);
    bd->alpha2_par = dot_product(bd->alpha2, bd->Lhat, 3);

    vec_project_perp(bd->alpha1, bd->Lhat, bd->alpha1_perp, 3);
    vec_project_perp(bd->alpha2, bd->Lhat, bd->alpha2_perp, 3);

    return 1;
}


// Innermost stable circular orbit for the remnant parameter description of Lousto et al. 2010/2012
double lz_E_ISCO_tilde(const lz_binary_data *bd)
{
    const double q = bd->q;
    const double eta = bd->eta;

    const double a1_par = bd->alpha1_par;
    const double a2_par = bd->alpha2_par;

    const double a1_dot_a1 = dot_product((double *)bd->alpha1, (double *)bd->alpha1, 3);
    const double a2_dot_a2 = dot_product((double *)bd->alpha2, (double *)bd->alpha2, 3);
    const double a1_dot_a2 = dot_product((double *)bd->alpha1, (double *)bd->alpha2, 3);

    const double term0 = 1.0 - sqrt(8.0 / 9.0);
    const double term_eta = 0.103803 * eta;

    const double term_spin_linear = 1.0 / (36.0 * sqrt(3.0) * (1.0 + q) * (1.0 + q))
        * (q * (1.0 + 2.0*q) * a1_par + (2.0 + q) * a2_par);

    const double spin_quad = a2_dot_a2 - 3.0 * a2_par * a2_par
        - 2.0*q * (a1_dot_a2 - 3.0 * a1_par * a2_par)
        + q*q * (a1_dot_a1 - 3.0 * a1_par * a1_par);

    const double term_spin_quad = -5.0 / (324.0 * sqrt(2.0) * (1.0 + q) * (1.0 + q))
        * spin_quad;

    return term0 + term_eta + term_spin_linear + term_spin_quad;
}


// J_ISCO for the remnant parameter description of Lousto et al. 2010/2012
void lz_J_ISCO_tilde_vector(lz_binary_data *bd, double Jvec[3])
{
    const double q = bd->q;
    const double eta = bd->eta;

    const double a1_par = bd->alpha1_par;
    const double a2_par = bd->alpha2_par;

    const double a1_dot_a1 = dot_product((double *)bd->alpha1, (double *)bd->alpha1, 3);
    const double a2_dot_a2 = dot_product((double *)bd->alpha2, (double *)bd->alpha2, 3);
    const double a1_dot_a2 = dot_product((double *)bd->alpha1, (double *)bd->alpha2, 3);

    const double base = 2.0 * sqrt(3.0) - 1.5255862 * eta;

    const double spin_linear_L =
        -1.0 / (9.0 * sqrt(2.0) * (1.0 + q) * (1.0 + q)) * (q * (7.0 + 8.0*q) * a1_par
        + (8.0 + 7.0*q) * a2_par);

    const double spin_quad =
        a2_dot_a2 - 3.0 * a2_par * a2_par - 2.0*q * (a1_dot_a2 - 3.0 * a1_par * a2_par)
        + q*q * (a1_dot_a1 - 3.0 * a1_par * a1_par);

    const double spin_quad_L = 2.0 / (9.0 * sqrt(3.0) * (1.0 + q) * (1.0 + q)) * spin_quad;

    const double L_coeff = base + spin_linear_L + spin_quad_L;

    for (int k = 0; k < 3; k++)
        Jvec[k] = L_coeff * bd->Lhat[k];

    for (int k = 0; k < 3; k++) {
        Jvec[k] += -1.0 / (9.0 * sqrt(2.0) * (1.0 + q) * (1.0 + q))
            * (q * (1.0 + 4.0*q) * bd->alpha1[k] + (4.0 + q) * bd->alpha2[k]);

        Jvec[k] += (bd->alpha2[k] + q*q * bd->alpha1[k]) / (eta * (1.0 + q) * (1.0 + q));
    }
}


/**
 * @brief Computes NR-informed radiated energy fraction for the merger of a binary.
 *
 * Based on Lousto et al. 2012.
 *
 * @param[in]   bd         Data of the merging binary
 * @returns The ene
 */
double lz_radiated_energy_fraction(const lz_binary_data *bd)
{
    const double q = bd->q;
    const double eta = bd->eta;

    const double a1_par = bd->alpha1_par;
    const double a2_par = bd->alpha2_par;

    double alpha_plus[3];
    double alpha_minus[3];

    double alpha_plus_perp[3];
    double alpha_minus_perp[3];

    for (int k = 0; k < 3; k++) {
        alpha_plus[k]  = bd->alpha2[k] + q * bd->alpha1[k];
        alpha_minus[k] = bd->alpha2[k] - q * bd->alpha1[k];

        alpha_plus_perp[k]  = bd->alpha2_perp[k] + q * bd->alpha1_perp[k];
        alpha_minus_perp[k] = bd->alpha2_perp[k] - q * bd->alpha1_perp[k];
    }

    const double alpha_plus2 = dot_product(alpha_plus, alpha_plus, 3);
    const double alpha_minus2 = dot_product(alpha_minus, alpha_minus, 3);

    const double alpha_plus_perp2 = dot_product(alpha_plus_perp, alpha_plus_perp, 3);
    const double alpha_minus_perp2 = dot_product(alpha_minus_perp, alpha_minus_perp, 3);

    // Phase angles of the in-plane spin combinations relative to the infall direction near merger
    // Approximated with phase average
    // TODO: Random phase and deterministic phase based on instantaneous separation direction e1
    const double cos2_theta_plus = 0.5;
    const double cos2_theta_minus = 0.5;

    const double Eisco = lz_E_ISCO_tilde(bd);

    const double spin_bracket = LZ_ES * (a2_par + q*q * a1_par)
    + LZ_EDELTA * (1.0 - q) * (a2_par - q * a1_par)
    + LZ_EA * alpha_plus2
    + LZ_EB * alpha_plus_perp2 * (cos2_theta_plus + LZ_EC)
    + LZ_ED * alpha_minus2
    + LZ_EE * alpha_minus_perp2 * (cos2_theta_minus + LZ_EF);

    double erad_frac = eta * Eisco
        + LZ_E2 * eta * eta
        + LZ_E3 * eta * eta * eta
        + eta * eta / ((1.0 + q) * (1.0 + q)) * spin_bracket;

    if (!isfinite(erad_frac))
        erad_frac = 0.0;

    erad_frac = clamp_double(erad_frac, 0.0, 0.3);

    return erad_frac;
}


/**
 * @brief Computes NR-informed spin for the merger remnant of a binary.
 *
 * Based on Lousto et al. 2012.
 *
 * @param[in]   bd         Data of the merging binary
 * @param[in]   erad_frac  Radiated energy fraction
 * @param[out]  chi_final  Remnant spin vector
 */
void lz_final_spin(lz_binary_data *bd, double erad_frac, double chi_final[3])
{
    const double q = bd->q;
    const double eta = bd->eta;

    const double a1_par = bd->alpha1_par;
    const double a2_par = bd->alpha2_par;

    double Jisco[3];
    lz_J_ISCO_tilde_vector(bd, Jisco);

    for (int k = 0; k < 3; k++)
        chi_final[k] = eta * Jisco[k];

    const double J_nonspin = LZ_J2 * eta * eta + LZ_J3 * eta * eta * eta;

    for (int k = 0; k < 3; k++)
        chi_final[k] += J_nonspin * bd->Lhat[k];

    const double J_spin = eta * eta / ((1.0 + q) * (1.0 + q))
        * (LZ_JA * (a2_par + q*q * a1_par) + LZ_JB * (1.0 - q) * (a2_par - q * a1_par));

    for (int k = 0; k < 3; k++)
        chi_final[k] += J_spin * bd->Lhat[k];

    const double mass_factor = 1.0 - erad_frac;

    if (mass_factor > 0.0) {
        const double scale = 1.0 / (mass_factor * mass_factor);

        for (int k = 0; k < 3; k++)
            chi_final[k] *= scale;
    }

    // Safety cap
    const double chi_mag = norm(chi_final, 3);
    if (chi_mag > 0.999) {
        for (int k = 0; k < 3; k++)
            chi_final[k] *= 0.999 / chi_mag;
    }
}


/**
 * @brief Computes NR-informed recoil/kick velocity for the merger remnant of a binary.
 *
 * Based on Lousto et al. 2012 (see Blecha et al. 2016 for formulas closer to this implementation).
 * The phase of the out-of-plane superkick term is not determined robustly by this simple
 * prescription. Here we approximate it using the instantaneous direction of Delta_perp relative
 * to e1. For statistical studies one often randomizes this phase instead.
 *
 * @param[in]   bd         Data of the merging binary
 * @param[out]  v_kick     Computed kick velocity
 * @param[out]  v_kick_kms Computed kick velocity in kms
 */
void lz_kick_velocity(lz_binary_data *bd, double v_kick[3], double v_kick_kms[3])
{
    const double q = bd->q;
    const double eta = bd->eta;

    for (int k = 0; k < 3; k++) {
        v_kick[k] = 0.0;
        v_kick_kms[k] = 0.0;
    }

    const double alpha1_parallel = bd->alpha1_par;
    const double alpha2_parallel = bd->alpha2_par;

    double Delta_perp[3];
    for (int k = 0; k < 3; k++)
        Delta_perp[k] = bd->alpha2_perp[k] - q * bd->alpha1_perp[k];
    const double Delta_perp_mag = norm(Delta_perp, 3);

    const double S_tilde_parallel = 2.0 * (alpha2_parallel + q*q * alpha1_parallel)
        / ((1.0 + q) * (1.0 + q));

    // Mass-asymmetry kick
    const double v_m_kms = LZ_KICK_A_KMS * eta * eta * (1.0 - q) / (1.0 + q)
        * (1.0 + LZ_KICK_B * eta);

    // In-plane spin kick from spin components parallel to L
    const double v_perp_kms = LZ_KICK_H_KMS * eta * eta / (1.0 + q)
        * (alpha2_parallel - q * alpha1_parallel);

    // Phase convention (the most model-dependent part)
    // Approximate cos(phi_Delta - phi_1) using the instantaneous in-plane separation direction e1
    double cos_phase = 0.0;

    if (Delta_perp_mag > 0.0) {
        cos_phase = dot_product(Delta_perp, (double *)bd->e1, 3) / Delta_perp_mag;
        cos_phase = clamp_double(cos_phase, -1.0, 1.0);
    }

    const double V_kms = LZ_KICK_V11_KMS
        + LZ_KICK_VA_KMS * S_tilde_parallel
        + LZ_KICK_VB_KMS * S_tilde_parallel * S_tilde_parallel
        + LZ_KICK_VC_KMS * S_tilde_parallel * S_tilde_parallel * S_tilde_parallel;

    const double v_parallel_kms = 16.0 * eta*eta / (1.0 + q) * V_kms * Delta_perp_mag * cos_phase;

    for (int k = 0; k < 3; k++) {
        v_kick_kms[k] = v_m_kms * bd->e1[k]
            + v_perp_kms * (cos(LZ_KICK_XI_RAD) * bd->e1[k] + sin(LZ_KICK_XI_RAD) * bd->e2[k])
            + v_parallel_kms * bd->Lhat[k];

        v_kick[k] = v_kick_kms[k] / C_KMS;
    }
}


/**
 * @brief Computes NR-informed remnant parameters for the merger of bodies i and j.
 *
 * Based on Lousto et al. 2012.
 *
 * @param[in]   params     Parameter struct containing general information about the system
 * @param[in]   w          State vector w = [positions, momenta, spins]
 * @param[in]   i          Index of the first merger object
 * @param[in]   j          Index of the second merger object
 * @param[out]  rem        Remnant struct
 * @return 1 if remnant parameters could be computed, 0 otherwise.
 */
int compute_lz_remnant(struct ode_params *params, double *w, int i, int j, lz_remnant *rem)
{
    if (rem == NULL)
        return 0;

    rem->mass = 0.0;
    rem->erad_frac = 0.0;

    for (int k = 0; k < 3; k++) {
        rem->chi[k] = 0.0;
        rem->v_kick[k] = 0.0;
        rem->v_kick_kms[k] = 0.0;
    }

    lz_binary_data bd;

    if (!get_lz_binary_data(params, w, i, j, &bd))
        return 0;

    rem->erad_frac = lz_radiated_energy_fraction(&bd);
    rem->mass = bd.M * (1.0 - rem->erad_frac);

    lz_final_spin(&bd, rem->erad_frac, rem->chi);
    lz_kick_velocity(&bd, rem->v_kick, rem->v_kick_kms);

    return 1;
}


// ------------------------------------------------------------------------------------------------
// Merger functions
// ------------------------------------------------------------------------------------------------

/**
 * @brief Finds the closest merging pair.
 *
 * Finds the closest active pair satisfying r_ij < MERGE_FACTOR * (m_i + m_j) which is
 * gravitationally bound according to the Newtonian binding energy.
 *
 * @param[in]   params     Parameter struct containing general information about the system
 * @param[in]   w          State vector w = [positions, momenta, spins]
 * @param[out]  i_merge    Index of the first merger object
 * @param[out]  j_merge    Index of the second merger object
 * @return 1 if such a pair is found, 0 otherwise. The pair is written to *i_merge and *j_merge.
 */
int find_merger_pair(struct ode_params *params, double *w, int *i_merge, int *j_merge)
{
    int found = 0;
    double best_dist2 = DBL_MAX;

    const int n_bodies = params->num_bodies_initial;
    const int num_dim = params->num_dim;

    for (int i = 0; i < n_bodies; i++) {
        if (!params->active[i])
            continue;

        for (int j = i + 1; j < n_bodies; j++) {
            if (!params->active[j])
                continue;

            const double mi = params->masses[i];
            const double mj = params->masses[j];
            const double M_tot = mi + mj;

            const double r_merge = params->merge_factor * M_tot;
            const double r_merge2 = r_merge * r_merge;

            double dist2 = 0.0;

            for (int n = 0; n < num_dim; n++) {
                const double dx = w[num_dim * i + n] - w[num_dim * j + n];
                dist2 += dx * dx;
            }

            if (dist2 < r_merge2 && dist2 < best_dist2 && pair_is_bound(params, w, i, j)) {
                best_dist2 = dist2;
                *i_merge = i;
                *j_merge = j;
                found = 1;
            }
        }
    }
    return found;
}


/**
 * @brief Merges active bodies i and j.
 *
 * @param   params        Parameter struct containing general information about the system
 * @param   w             State vector w = [positions, momenta, spins]
 * @param   i             Index of the first merger object
 * @param   j             Index of the second merger object
 * @param   t             Current time
 * @param   file_merger   Pointer to the merger output file
 */
void merge_pair(struct ode_params *params, double *w, int i, int j, double t, FILE *file_merger)
{
    const int num_dim = params->num_dim;
    const int array_half = num_dim * params->num_bodies_initial;
    const int spin_offset = 2 * array_half;

    if (!params->active[i] || !params->active[j])
        return;

    if (i == j)
        return;

    const double mi = params->masses[i];
    const double mj = params->masses[j];

    // Compute the remnant mass, spin and kick velocity based on a selected prescription
    lz_remnant rem;
    init_simple_remnant(&rem, mi + mj);

    if (strcmp(params->remnant_prescription, "lz") == 0) {
        if (!compute_lz_remnant(params, w, i, j, &rem)) {
            printf("Warning: Lousto-Zlochower remnant parameters could not be computed!\n");
            printf("Using simple mass addition without remnant kick and spin instead!\n");
            init_simple_remnant(&rem, mi + mj);
        }
    }
    else if (strcmp(params->remnant_prescription, "barausse") == 0) {
        rem.mass = barausse_remnant_mass(params, w, i, j);

        if (!(rem.mass > 0.0) || !isfinite(rem.mass)) {
            printf("Warning: Barausse remnant mass could not be computed!\n");
            printf("Using simple mass addition without remnant kick and spin instead!\n");
            init_simple_remnant(&rem, mi + mj);
        } else {
            rem.erad_frac = 1.0 - rem.mass / (mi + mj);

            lz_binary_data bd;
            if (get_lz_binary_data(params, w, i, j, &bd)) {
                lz_final_spin(&bd, rem.erad_frac, rem.chi);
                lz_kick_velocity(&bd, rem.v_kick, rem.v_kick_kms);
            }
        }
    }

    const long long id_i = params->body_id[i];
    const long long id_j = params->body_id[j];
    const long long id_remnant = params->next_body_id++;

    const int gen_i = params->generation[i];
    const int gen_j = params->generation[j];
    const int gen_remnant = 1 + (gen_i > gen_j ? gen_i : gen_j);

    double x_i_old[num_dim];
    double x_j_old[num_dim];
    double x_rem[num_dim];

    double p_i_old[num_dim];
    double p_j_old[num_dim];
    double p_rem[num_dim];

    double s_i_old[3];
    double s_j_old[3];
    double s_rem[3];

    double r2 = 0.0;

    // Save old parent states and compute remnant state
    for (int n = 0; n < num_dim; n++) {
        const int xi = num_dim * i + n;
        const int xj = num_dim * j + n;

        const int pi = array_half + num_dim * i + n;
        const int pj = array_half + num_dim * j + n;

        x_i_old[n] = w[xi];
        x_j_old[n] = w[xj];

        p_i_old[n] = w[pi];
        p_j_old[n] = w[pj];

        const double dx = x_i_old[n] - x_j_old[n];
        r2 += dx * dx;

        const double M_initial = mi + mj;
        x_rem[n] = (mi * x_i_old[n] + mj * x_j_old[n]) / M_initial;
        p_rem[n] = p_i_old[n] + p_j_old[n] + rem.mass * rem.v_kick[n];
    }
    for (int n = 0; n < 3; n++) {
        const int si = spin_offset + 3 * i + n;
        const int sj = spin_offset + 3 * j + n;

        s_i_old[n] = w[si];
        s_j_old[n] = w[sj];
        s_rem[n] = rem.chi[n];
    }

    const double r_ij = sqrt(r2);

    // Write merger event before overwriting bookkeeping
    if (file_merger != NULL) {
        output_write_merger_event(
            file_merger,
            t,
            params,
            i,
            j,
            i,
            id_i,
            id_j,
            id_remnant,
            gen_i,
            gen_j,
            gen_remnant,
            mi,
            mj,
            rem.mass,
            x_i_old,
            x_j_old,
            x_rem,
            p_i_old,
            p_j_old,
            p_rem,
            s_i_old,
            s_j_old,
            s_rem,
            rem.v_kick_kms,
            r_ij
        );
    }

    // Update state vector
    for (int n = 0; n < num_dim; n++) {
        const int xi = num_dim * i + n;
        const int xj = num_dim * j + n;

        const int pi = array_half + num_dim * i + n;
        const int pj = array_half + num_dim * j + n;

        w[xi] = x_rem[n];
        w[pi] = p_rem[n];

        w[xj] = x_rem[n];
        w[pj] = 0.0;
    }
    for (int n = 0; n < 3; n++) {
        const int si = spin_offset + 3 * i + n;
        const int sj = spin_offset + 3 * j + n;

        w[si] = s_rem[n];
        w[sj] = 0.0;
    }

    // Update merger history params
    params->masses[i] = rem.mass;
    params->masses[j] = 0.0;

    params->active[i] = 1;
    params->active[j] = 0;
    params->num_active--;

    params->body_id[i] = id_remnant;
    params->body_id[j] = -1;

    params->generation[i] = gen_remnant;
    params->generation[j] = -1;
}


/**
 * @brief Repeatedly searches for merger pairs and merges them.
 *
 * This handles triple/quadruple situations by merging one pair, then re-running the search on the
 * updated system.
 *
 * @param   params       Parameter struct containing general information about the system
 * @param   w            State vector w = [positions, momenta, spins]
 * @param   t            Current time
 * @param   file_merger  Pointer to the merger output file
 * @return 1 if at least one merger happened, 0 otherwise.
 */
int test_and_merge_bodies(struct ode_params *params, double *w, double t, FILE *file_merger)
{
    int merged_any = 0;

    while (1) {
        int i_merge = -1;
        int j_merge = -1;

        const int found = find_merger_pair(params, w, &i_merge, &j_merge);

        if (!found)
            break;

        merge_pair(params, w, i_merge, j_merge, t, file_merger);
        merged_any = 1;
    }

    return merged_any;
}
