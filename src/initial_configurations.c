/**
 * @file initial_configurations.c
 * @brief Functions for returning the intial positions and momenta for specific systems.
 *
 * Functions that can be used to input high-level properties of a system, from which
 * the corresponding initial positions and momenta of all the bodies are computed.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "initial_configurations.h"
#include "utils.h"
#include "eom.h"
#include "pn_binary.h"


// ------------------------------------------------------------------------------------------------
// PN initial configuration helper functions
// ------------------------------------------------------------------------------------------------

/**
 * @brief Check whether all bodies satisfy a minimum pair separation.
 *
 * Tests every unordered body pair against the mass-dependent separation criterion
 *
 *     r_ij >= min_sep_factor * (m_i + m_j).
 *
 * The positions are read from the flattened array as x[i * D + k]. The check includes all
 * num_bodies_initial bodies and returns immediately when a pair violating the criterion is found.
 * A nonpositive min_sep_factor disables the check.
 *
 * @param[in] params         Simulation parameters.
 * @param[in] x              Flattened position array with N * D elements.
 * @param[in] min_sep_factor Factor relating the minimum separation to the combined mass of a pair.
 *
 * @return 1 if the check is disabled or every pair satisfies the minimum separation;
 *         0 if at least one pair is too close.
 */
static int min_sep_ok(struct ode_params* params, double* x, double min_sep_factor)
{
    const int N = params->num_bodies_initial;
    const int D = params->num_dim;

    if (min_sep_factor <= 0.0)
        return 1;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double rij2 = 0.0;

            for (int k = 0; k < D; k++) {
                const double dx = x[i * D + k] - x[j * D + k];
                rij2 += dx * dx;
            }

            const double rij = sqrt(rij2);
            const double min_sep = min_sep_factor * (params->masses[i] + params->masses[j]);

            if (rij < min_sep)
                return 0;
        }
    }

    return 1;
}


/**
 * @brief Changes the position and full 3D orientation of a binary.
 *
 * The binary is initially generated in the canonical orbital frame:
 *
 *     x-axis: periapsis/reference direction
 *     y-axis: positive true-anomaly direction
 *     z-axis: angular-momentum direction
 *
 * This function maps that basis to:
 *
 *     e_hat: periapsis/reference direction
 *     q_hat = h_hat x e_hat
 *     h_hat: orbital angular-momentum direction
 *
 * For circular binaries, e_hat is simply a phase-reference direction.
 */
void position_binary(double com_pos[3], double h_input[3], double e_input[3], double w0[12])
{
    double h_hat[3] = {h_input[0], h_input[1], h_input[2]};
    double e_hat[3] = {e_input[0], e_input[1], e_input[2]};
    double q_hat[3];

    double pos1[3] = {w0[0], w0[1], w0[2]};
    double pos2[3] = {w0[3], w0[4], w0[5]};
    double p1[3]   = {w0[6], w0[7], w0[8]};
    double p2[3]   = {w0[9], w0[10], w0[11]};

    double pos1_new[3], pos2_new[3], p1_new[3], p2_new[3];

    // q_hat completes a right-handed orbital basis
    cross_product(h_hat, e_hat, q_hat);
    normalize(q_hat, q_hat);

    // Rotate positions and momenta from canonical frame to target frame
    map_from_orbital_basis(pos1, e_hat, q_hat, h_hat, pos1_new);
    map_from_orbital_basis(pos2, e_hat, q_hat, h_hat, pos2_new);
    map_from_orbital_basis(p1,   e_hat, q_hat, h_hat, p1_new);
    map_from_orbital_basis(p2,   e_hat, q_hat, h_hat, p2_new);

    // Shift positions to target COM
    w0[0] = pos1_new[0] + com_pos[0];
    w0[1] = pos1_new[1] + com_pos[1];
    w0[2] = pos1_new[2] + com_pos[2];

    w0[3] = pos2_new[0] + com_pos[0];
    w0[4] = pos2_new[1] + com_pos[1];
    w0[5] = pos2_new[2] + com_pos[2];

    // Momenta are rotated but not shifted
    w0[6] = p1_new[0];
    w0[7] = p1_new[1];
    w0[8] = p1_new[2];

    w0[9]  = p2_new[0];
    w0[10] = p2_new[1];
    w0[11] = p2_new[2];
}


/**
 * @brief Construct the incoming canonical momentum of the left scattering subsystem.
 *
 * The magnitude p0_rel is interpreted as the total canonical COM momentum of either subsystem,
 * not as a relative velocity or as a per-body momentum. The right subsystem receives the opposite
 * momentum. For positive b, the left subsystem initially moves in the +x, -y direction.
 */
static void scattering_com_momentum(double d0, double p0_rel, double b, double P[3])
{
    const double b_over_d0 = b / d0;

    P[0] = p0_rel * sqrt(fmax(0.0, 1.0 - b_over_d0 * b_over_d0));
    P[1] = -p0_rel * b_over_d0;
    P[2] = 0.0;
}


/**
 * @brief Place two subsystem centers a distance d0 apart in their Newtonian total-COM frame.
 */
static void place_scattering_coms(double d0, double mass_left, double mass_right,
    double left_com[3], double right_com[3])
{
    const double total_mass = mass_left + mass_right;

    left_com[0] = -(mass_right / total_mass) * d0;
    left_com[1] = 0.0;
    left_com[2] = 0.0;

    right_com[0] = (mass_left / total_mass) * d0;
    right_com[1] = 0.0;
    right_com[2] = 0.0;
}


/**
 * @brief Add a binary COM momentum while preserving its Newtonian relative canonical momentum.
 *
 * Splitting P in proportion to the component masses is the linear canonical Jacobi decomposition.
 * It is the appropriate leading-order superposition for these scattering presets, but it is not
 * an exact finite-momentum PN Poincare boost.
 */
static void add_binary_com_momentum(double w_binary[12], double m1, double m2,
    const double P[3])
{
    const double M = m1 + m2;

    for (int k = 0; k < 3; k++) {
        w_binary[6 + k] += (m1 / M) * P[k];
        w_binary[9 + k] += (m2 / M) * P[k];
    }
}


// ------------------------------------------------------------------------------------------------
// Initial configurations
// ------------------------------------------------------------------------------------------------


/**
 * @brief Computes initial positions and momenta for a binary.
 *
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, p1, p2]) of two
 * bodies with masses m1 and m2 in a binary system with relative semi-major axis a,
 * eccentricity e, and true anomaly f0 in the center of mass frame and in the xy-plane.
 *
 * @param[in]   params          Parameter struct containing general information about the system
 * @param[in]   binary_params   Struct containing the binary parameters (initialized elsewhere)
 * @param[out]  w0              Initial positions and momenta, w0 = [pos1, pos2, p1, p2]
 */
void ic_binary(struct ode_params* ode_params, struct binary_params* binary_params,
    double m1, double m2, double* w0)
{
    // Unpack needed parameters
    const double a = binary_params->a;
    const double e = binary_params->e;
    const double f0 = binary_params->f0;
    const double p = binary_params->p;
    const int use_1pn = (ode_params->pn_terms[1] == 1);
    const int use_2pn = (ode_params->pn_terms[2] == 1);

    // Compute total mass, reduced mass and symmetric mass ratio
    const double M = m1 + m2;
    const double mu = m1 * m2 / M;
    const double nu = mu / M;

    // Geometric quantities
    const double cosf0 = cos(f0);
    const double sinf0 = sin(f0);

    const double denom = 1.0 + e * cosf0;
    if (denom <= 0.0)
        errorexit("invalid phase/eccentricity combination");

    const double r0 = p / denom;
    const double x0 = r0 / M;

    const double rp = a * (1.0 - e);
    const double ra = a * (1.0 + e);
    const double xp = rp / M;
    const double xa = ra / M;

    if (rp <= 0.0 || ra <= 0.0 || x0 <= 0.0)
        errorexit("non-positive PN radius");

    double j;
    if (e < 1e-12) {
        j = pn_binary_solve_circular_j(a / M, nu, use_1pn, use_2pn);
        if (!isfinite(j))
            errorexit("could not bracket circular angular momentum");
    }
    else {
        j = pn_binary_solve_eccentric_j(xp, xa, nu, use_1pn, use_2pn);
        if (!isfinite(j))
            errorexit("could not bracket eccentric angular momentum");
    }

    const double energy = (e < 1e-12)
        ? pn_binary_turning_hamiltonian(a / M, j, nu, use_1pn, use_2pn)
        : pn_binary_turning_hamiltonian(xp, j, nu, use_1pn, use_2pn);

    double pr_hat_abs = pn_binary_solve_pr_hat_abs(
        x0, j, energy, nu, use_1pn, use_2pn);
    if (!isfinite(pr_hat_abs))
        errorexit("could not solve radial momentum");
    double pr_hat;

    if (e < 1e-12 && ode_params->pn_terms[3] == 1) {
        // Quasi-circular inspiral radial momentum
        const double xdot = -(64.0 / 5.0) * nu / pow(x0, 3);

        const double A = pn_binary_dh_dpr_coeff_at_zero(
            x0, j, nu, use_1pn, use_2pn);

        if (!isfinite(A) || fabs(A) < 1e-14)
            errorexit("could not compute circular 2.5PN radial momentum");

        pr_hat = xdot / A;
    }
    else {
        // Conservative eccentric-orbit radial momentum
        double pr_sign = 0.0;

        if (sinf0 > 1e-14)
            pr_sign = 1.0;
        else if (sinf0 < -1e-14)
            pr_sign = -1.0;

        pr_hat = pr_sign * pr_hat_abs;
    }
    const double pt_hat = j / x0;

    // Unit vectors in the orbital plane
    // n points radially outward, lambda is the positive-phi tangential unit vector
    const double nx = cosf0;
    const double ny = sinf0;
    const double lx = -sinf0;
    const double ly = cosf0;

    const double rvec_x = r0 * nx;
    const double rvec_y = r0 * ny;

    // Physical relative canonical momentum P_rel = mu * p_hat
    const double prel_x = mu * (pr_hat * nx + pt_hat * lx);
    const double prel_y = mu * (pr_hat * ny + pt_hat * ly);

    const double pos1_x =  m2 / M * rvec_x;
    const double pos1_y =  m2 / M * rvec_y;
    const double pos2_x = -m1 / M * rvec_x;
    const double pos2_y = -m1 / M * rvec_y;

    if (ode_params->num_dim == 2) {
        w0[0] = pos1_x;
        w0[1] = pos1_y;

        w0[2] = pos2_x;
        w0[3] = pos2_y;

        w0[4] =  prel_x;
        w0[5] =  prel_y;

        w0[6] = -prel_x;
        w0[7] = -prel_y;
    }
    else {
        w0[0] = pos1_x;
        w0[1] = pos1_y;
        w0[2] = 0.0;

        w0[3] = pos2_x;
        w0[4] = pos2_y;
        w0[5] = 0.0;

        w0[6] =  prel_x;
        w0[7] =  prel_y;
        w0[8] = 0.0;

        w0[9]  = -prel_x;
        w0[10] = -prel_y;
        w0[11] = 0.0;
    }
}


/**
 * @brief Computes initial positions and momenta for a hierarchical triple.
 *
 * The system is treated as an inner binary plus an outer binary:
 *      inner binary:  body 1 + body 2
 *      outer binary:  inner-COM + body 3
 *
 * @param[in]   ode_params                General ODE/system parameters
 * @param[in]   inner_binary_params       Orbital parameters of bodies 1 and 2
 * @param[in]   outer_binary_params       Orbital parameters of inner-COM and body 3
 * @param[out]  w0                        Initial positions and momenta
 */
void ic_hierarchical_triple(struct ode_params* ode_params,
    struct binary_params* inner_binary_params, struct binary_params* outer_binary_params,
    double* w0)
{
    double m1 = ode_params->masses[0];
    double m2 = ode_params->masses[1];
    double m3 = ode_params->masses[2];
    double m_inner = m1 + m2;

    double w_inner[12];
    double w_outer[12];

    // Initialize the inner binary in its own COM frame
    ic_binary(ode_params, inner_binary_params, m1, m2, w_inner);

    // Initialize the effective outer binary
    ic_binary(ode_params, outer_binary_params, m_inner, m3, w_outer);

    // Orient the outer binary around the total triple COM
    double triple_com_pos[3] = {0.0, 0.0, 0.0};
    position_binary(triple_com_pos, outer_binary_params->h_hat, outer_binary_params->e_hat, w_outer);

    double inner_com_pos[3] = {w_outer[0], w_outer[1], w_outer[2]};
    double inner_com_mom[3] = {w_outer[6], w_outer[7], w_outer[8]};
    double tertiary_pos[3] = {w_outer[3], w_outer[4], w_outer[5]};
    double tertiary_mom[3] = {w_outer[9], w_outer[10], w_outer[11]};

    // Orient the inner binary and place its COM on the outer orbit
    position_binary(inner_com_pos, inner_binary_params->h_hat, inner_binary_params->e_hat, w_inner);

    // Add the outer COM momentum to the two inner bodies
    w_inner[6]  += (m1 / m_inner) * inner_com_mom[0];
    w_inner[7]  += (m1 / m_inner) * inner_com_mom[1];
    w_inner[8]  += (m1 / m_inner) * inner_com_mom[2];

    w_inner[9]  += (m2 / m_inner) * inner_com_mom[0];
    w_inner[10] += (m2 / m_inner) * inner_com_mom[1];
    w_inner[11] += (m2 / m_inner) * inner_com_mom[2];

    // Fill final state vector
    if (ode_params->num_dim == 2) {
        w0[0] = w_inner[0];
        w0[1] = w_inner[1];

        w0[2] = w_inner[3];
        w0[3] = w_inner[4];

        w0[4] = tertiary_pos[0];
        w0[5] = tertiary_pos[1];

        w0[6] = w_inner[6];
        w0[7] = w_inner[7];

        w0[8] = w_inner[9];
        w0[9] = w_inner[10];

        w0[10] = tertiary_mom[0];
        w0[11] = tertiary_mom[1];
    }
    else if (ode_params->num_dim == 3) {
        w0[0] = w_inner[0];
        w0[1] = w_inner[1];
        w0[2] = w_inner[2];

        w0[3] = w_inner[3];
        w0[4] = w_inner[4];
        w0[5] = w_inner[5];

        w0[6] = tertiary_pos[0];
        w0[7] = tertiary_pos[1];
        w0[8] = tertiary_pos[2];

        w0[9]  = w_inner[6];
        w0[10] = w_inner[7];
        w0[11] = w_inner[8];

        w0[12] = w_inner[9];
        w0[13] = w_inner[10];
        w0[14] = w_inner[11];

        w0[15] = tertiary_mom[0];
        w0[16] = tertiary_mom[1];
        w0[17] = tertiary_mom[2];
    }
}


/**
 * @brief Computes initial positions and momenta for a binary-single scattering.
 *
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, pos3, p1, p2, p3])
 * of a binary-single scattering with specified binary and scattering parameters. The
 * scattering takes place in the xy-plane.
 *
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[in]   binary_params   Struct containing the binary parameters (initialized elsewhere)
 * @param[in]   d0              Initial distance of the binary's center of mass and the single
 * @param[in]   p0_rel          Total canonical COM momentum magnitude of either scattering subsystem
 * @param[in]   b               Scattering impact parameter
 * @param[out]  w0              Initial positions and momenta,
 *                              w0 = [pos_b1, pos_b2, pos_s, p_b1, p_b2, p_s]
 */
void ic_binary_single_scattering(struct ode_params* ode_params,
    struct binary_params* binary_params, double d0, double p0_rel, double b, double* w0)
{
    const double m1 = ode_params->masses[0];
    const double m2 = ode_params->masses[1];
    const double m3 = ode_params->masses[2];
    const double binary_mass = m1 + m2;

    double w0_binary[12];
    double P[3];
    double binary_com[3];
    double single_pos[3];

    // Create the binary initial values
    ic_binary(ode_params, binary_params, m1, m2, w0_binary);

    scattering_com_momentum(d0, p0_rel, b, P);
    place_scattering_coms(d0, binary_mass, m3, binary_com, single_pos);

    // Position the binary at its right place with the specified orientation
    position_binary(binary_com, binary_params->h_hat, binary_params->e_hat, w0_binary);

    // Give the binary total momentum P without changing its relative canonical momentum
    add_binary_com_momentum(w0_binary, m1, m2, P);

    // Setting up the binary
    for(int i = 0; i < 6; i++)
        w0[i] = w0_binary[i];
    for(int i = 9; i < 15; i++)
        w0[i] = w0_binary[i-3];

    // Setting up the single body
    w0[6] = single_pos[0];
    w0[7] = single_pos[1];
    w0[8] = single_pos[2];
    w0[15] = -P[0];
    w0[16] = -P[1];
    w0[17] = -P[2];
}


/**
 * @brief Computes initial positions and momenta for a circular binary-single scattering with
 * low-eccentricity tangential and radial momenta.
 *
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, pos3, p1, p2, p3])
 * of a circular relativistic binary-single scattering with specified binary and scattering
 * parameters. The scattering takes place in the xy-plane. The function uses r0, pt0 and pr0
 * (the initial separation, tangential and radial momentum), which are commonly quoted as
 * relativistic low-eccentricity binary initial parameters.
 *
 * @param[in]   d0              Initial distance of the binary's center of mass and the single
 * @param[in]   p0_rel          Total canonical COM momentum magnitude of either scattering subsystem
 * @param[in]   b               Scattering impact parameter
 * @param[in]   binary_phi0     Binary phase shift
 * @param[in]   binary_r0       Initial binary radius
 * @param[in]   binary_pt0      Tangential component of the initial momentum of the binary members
 * @param[in]   binary_pr0      Radial component of the initial momentum of the binary members
 * @param[in]   orientation     Binary orientation vector (NULL -> oriented in direction of motion)
 * @param[out]  w0              Initial positions and momenta,
 *                              w0 = [pos_b1, pos_b2, pos_s, p_b1, p_b2, p_s]
 */
void ic_binary_single_scattering_circ(double d0, double p0_rel, double b, double binary_phi0,
    double binary_r0, double binary_pt0, double binary_pr0, double* orientation, double* w0)
{
    double binary_orientation[3];
    double P[3];

    const double b_over_d0 = b / d0;
    const double direction_x = sqrt(fmax(0.0, 1.0 - b_over_d0 * b_over_d0));
    const double direction_y = -b_over_d0;

    scattering_com_momentum(d0, p0_rel, b, P);

    // If no orientation is given, the binary will be oriented along its momentum vector
    if (orientation == NULL) {
        binary_orientation[0] = direction_x;
        binary_orientation[1] = direction_y;
        binary_orientation[2] = 0.0;
    }
    else {
        binary_orientation[0] = orientation[0];
        binary_orientation[1] = orientation[1];
        binary_orientation[2] = orientation[2];
    }

    // As a starting point set up the binary in the xy-plane at the origin
    double pos1[3] = {binary_r0, 0.0, 0.0};
    double pos2[3] = {-binary_r0, 0.0, 0.0};
    double p1[3] = {binary_pr0, -binary_pt0, 0.0};
    double p2[3] = {-binary_pr0, binary_pt0, 0.0};

    // Perform a rotation to account for the specified phase shift
    double R[3][3];
    double axis[3] = {0.0, 0.0, 1.0};
    create_rotation_matrix(axis, binary_phi0, R);
    rotate_vector(pos1, R, pos1);
    rotate_vector(pos2, R, pos2);
    rotate_vector(p1, R, p1);
    rotate_vector(p2, R, p2);

    // Perform a rotation to align the orbital plane with the orientation vector
    align_vectors_rotation_matrix(axis, binary_orientation, R);
    rotate_vector(pos1, R, pos1);
    rotate_vector(pos2, R, pos2);
    rotate_vector(p1, R, p1);
    rotate_vector(p2, R, p2);

    // Put the binary at the specified position and add the momentum vector for the scattering
    pos1[0] -= d0/2;
    pos2[0] -= d0/2;
    p1[0] += 0.5 * P[0];
    p1[1] += 0.5 * P[1];
    p2[0] += 0.5 * P[0];
    p2[1] += 0.5 * P[1];

    // Fill the state vector
    w0[0] = pos1[0];
    w0[1] = pos1[1];
    w0[2] = pos1[2];
    w0[3] = pos2[0];
    w0[4] = pos2[1];
    w0[5] = pos2[2];
    w0[6] = d0/2;
    w0[7] = 0.0;
    w0[8] = 0.0;
    w0[9] = p1[0];
    w0[10] = p1[1];
    w0[11] = p1[2];
    w0[12] = p2[0];
    w0[13] = p2[1];
    w0[14] = p2[2];
    w0[15] = -P[0];
    w0[16] = -P[1];
    w0[17] = -P[2];
}


/**
 * @brief Computes initial positions and momenta for a binary-binary scattering.
 *
 * Returns the initial positions and momenta (as w0 = [pos1, pos2, pos3, pos4, p1, p2, p3, p4])
 * of a binary-binary scattering with specified binary and scattering parameters. The
 * scattering takes place in the xy-plane.
 *
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[in]   binary1_params  Struct with the parameters of binary 1 (initialized elsewhere)
 * @param[in]   binary2_params  Struct with the parameters of binary 2 (initialized elsewhere)
 * @param[in]   d0              Initial distance of the centers of mass of the binaries
 * @param[in]   p0_rel          Total canonical COM momentum magnitude of either binary
 * @param[in]   b               Scattering impact parameter
 * @param[out]  w0              Initial positions and momenta,
 *                              w0 = [pos1_1, pos2_1, pos1_2, pos2_2, p1_1, p2_1, p1_2, p2_2]
 */
void ic_binary_binary_scattering(struct ode_params* ode_params,
    struct binary_params* binary1_params, struct binary_params* binary2_params, double d0,
    double p0_rel, double b, double* w0)
{
    const double m1 = ode_params->masses[0];
    const double m2 = ode_params->masses[1];
    const double m3 = ode_params->masses[2];
    const double m4 = ode_params->masses[3];
    const double binary1_mass = m1 + m2;
    const double binary2_mass = m3 + m4;

    double w0_binary1[12];
    double w0_binary2[12];
    double P[3];
    double binary1_com[3];
    double binary2_com[3];

    // Create the binary initial values
    ic_binary(ode_params, binary1_params, m1, m2, w0_binary1);
    ic_binary(ode_params, binary2_params, m3, m4, w0_binary2);

    scattering_com_momentum(d0, p0_rel, b, P);
    place_scattering_coms(d0, binary1_mass, binary2_mass, binary1_com, binary2_com);

    // Position the binaries at their right place with the specified orientations
    position_binary(binary1_com, binary1_params->h_hat, binary1_params->e_hat, w0_binary1);
    position_binary(binary2_com, binary2_params->h_hat, binary2_params->e_hat, w0_binary2);

    // Give the binaries equal and opposite total canonical COM momenta
    add_binary_com_momentum(w0_binary1, m1, m2, P);

    const double minus_P[3] = {-P[0], -P[1], -P[2]};
    add_binary_com_momentum(w0_binary2, m3, m4, minus_P);

    // Setting up the binaries
    for(int i = 0; i < 6; i++)
        w0[i] = w0_binary1[i];
    for(int i = 6; i < 12; i++)
        w0[i] = w0_binary2[i-6];
    for(int i = 12; i < 18; i++)
        w0[i] = w0_binary1[i-6];
    for(int i = 18; i < 24; i++)
        w0[i] = w0_binary2[i-12];
}


/**
 * @brief Computes initial positions and momenta for a relativistic binary-binary scattering.
 *
 * Returns the initial positions and momenta (as w0 = [pos1, pos2, pos3, pos4, p1, p2, p3, p4])
 * of a relativistic binary-binary scattering with specified binary and scattering parameters. The
 * scattering takes place in the xy-plane. The function currently only supports equal-mass binaries
 * for which it is common to report the radial and tangential component of the binary members that
 * for example produce a quasi-circular orbit.
 *
 * @param[in]   d0              Initial distance of the binary's center of mass and the single
 * @param[in]   p0_rel          Total canonical COM momentum magnitude of either binary
 * @param[in]   b               Scattering impact parameter
 * @param[in]   binary_phi0_X   Phase shift of binary X
 * @param[in]   binary_r0_X     Initial radius of binary X
 * @param[in]   binary_pt0_X    Tangential component of the initial momentum of binary X members
 * @param[in]   binary_pr0_X    Radial component of the initial momentum of binary X members
 * @param[in]   orientation_X   Orientation of binary X (NULL -> oriented in direction of motion)
 * @param[out]  w0              Initial positions and momenta,
 *                              w0 = [pos1_1, pos2_1, pos1_2, pos2_2, p1_1, p2_1, p1_2, p2_2]
 */
void ic_binary_binary_scattering_circ(double d0, double p0_rel, double b, double binary1_phi0,
    double binary1_r0, double binary1_pt0, double binary1_pr0, double* orientation_1,
    double binary_phi0_2, double binary2_r0, double binary2_pt0, double binary2_pr0,
    double* orientation_2, double* w0)
{
    double binary_orientation_1[3];
    double binary_orientation_2[3];
    double P[3];

    // Compute the x- and y- components of the momentum vector of the two systems
    const double b_over_d0 = b / d0;
    const double direction_x = sqrt(fmax(0.0, 1.0 - b_over_d0 * b_over_d0));
    const double direction_y = -b_over_d0;

    scattering_com_momentum(d0, p0_rel, b, P);

    if (orientation_1 == NULL) {
        binary_orientation_1[0] = direction_x;
        binary_orientation_1[1] = direction_y;
        binary_orientation_1[2] = 0.0;
    }
    else {
        binary_orientation_1[0] = orientation_1[0];
        binary_orientation_1[1] = orientation_1[1];
        binary_orientation_1[2] = orientation_1[2];
    }

    if (orientation_2 == NULL) {
        binary_orientation_2[0] = -direction_x;
        binary_orientation_2[1] = -direction_y;
        binary_orientation_2[2] = 0.0;
    }
    else {
        binary_orientation_2[0] = orientation_2[0];
        binary_orientation_2[1] = orientation_2[1];
        binary_orientation_2[2] = orientation_2[2];
    }

    // As a starting point set up the binaries in the xy-plane at the origin
    double pos1_1[3] = {binary1_r0, 0.0, 0.0};
    double pos2_1[3] = {-binary1_r0, 0.0, 0.0};
    double p1_1[3] = {binary1_pr0, -binary1_pt0, 0.0};
    double p2_1[3] = {-binary1_pr0, binary1_pt0, 0.0};

    double pos1_2[3] = {binary2_r0, 0.0, 0.0};
    double pos2_2[3] = {-binary2_r0, 0.0, 0.0};
    double p1_2[3] = {binary2_pr0, -binary2_pt0, 0.0};
    double p2_2[3] = {-binary2_pr0, binary2_pt0, 0.0};

    // Perform a rotation to account for the specified phase shift
    double R[3][3];
    double axis[3] = {0.0, 0.0, 1.0};
    create_rotation_matrix(axis, binary1_phi0, R);
    rotate_vector(pos1_1, R, pos1_1);
    rotate_vector(pos2_1, R, pos2_1);
    rotate_vector(p1_1, R, p1_1);
    rotate_vector(p2_1, R, p2_1);
    create_rotation_matrix(axis, binary_phi0_2, R);
    rotate_vector(pos1_2, R, pos1_2);
    rotate_vector(pos2_2, R, pos2_2);
    rotate_vector(p1_2, R, p1_2);
    rotate_vector(p2_2, R, p2_2);

    // Perform a rotation to align the orbital plane with the orientation vector
    align_vectors_rotation_matrix(axis, binary_orientation_1, R);
    rotate_vector(pos1_1, R, pos1_1);
    rotate_vector(pos2_1, R, pos2_1);
    rotate_vector(p1_1, R, p1_1);
    rotate_vector(p2_1, R, p2_1);
    align_vectors_rotation_matrix(axis, binary_orientation_2, R);
    rotate_vector(pos1_2, R, pos1_2);
    rotate_vector(pos2_2, R, pos2_2);
    rotate_vector(p1_2, R, p1_2);
    rotate_vector(p2_2, R, p2_2);

    // Put the binaries at the specified position and add the momentum vector for the scattering
    pos1_1[0] -= d0/2;
    pos2_1[0] -= d0/2;
    p1_1[0] += 0.5 * P[0];
    p1_1[1] += 0.5 * P[1];
    p2_1[0] += 0.5 * P[0];
    p2_1[1] += 0.5 * P[1];
    pos1_2[0] += d0/2;
    pos2_2[0] += d0/2;
    p1_2[0] -= 0.5 * P[0];
    p1_2[1] -= 0.5 * P[1];
    p2_2[0] -= 0.5 * P[0];
    p2_2[1] -= 0.5 * P[1];

    // Fill the state vector
    w0[0] = pos1_1[0];
    w0[1] = pos1_1[1];
    w0[2] = pos1_1[2];
    w0[3] = pos2_1[0];
    w0[4] = pos2_1[1];
    w0[5] = pos2_1[2];
    w0[6] = pos1_2[0];
    w0[7] = pos1_2[1];
    w0[8] = pos1_2[2];
    w0[9] = pos2_2[0];
    w0[10] = pos2_2[1];
    w0[11] = pos2_2[2];
    w0[12] = p1_1[0];
    w0[13] = p1_1[1];
    w0[14] = p1_1[2];
    w0[15] = p2_1[0];
    w0[16] = p2_1[1];
    w0[17] = p2_1[2];
    w0[18] = p1_2[0];
    w0[19] = p1_2[1];
    w0[20] = p1_2[2];
    w0[21] = p2_2[0];
    w0[22] = p2_2[1];
    w0[23] = p2_2[2];
}


/**
 * @brief Computes initial positions and momenta for a three-body figure eight orbit.
 *
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, pos3, p1, p2, p3])
 * for a stable periodic three-body figure eight orbit in the xy-plane. This can be achieved for
 * post-Newtonian approximations up to 2PN using initial parameters for px and py fitted from
 * Lousto & Nakano (2008).
 *
 * @param[in]   params      Parameter struct containing general information about the system
 * @param[in]   width       Width of the figure eight orbit
 * @param[out]  w0          Initial positions and momenta, w0 = [pos1, pos2, pos3, p1, p2, p3]
 */
void ic_figure_eight_orbit(struct ode_params* params, double width, double* w0)
{
    double pos_x, pos_y, lambda;
    double px = 0.0;
    double py = 0.0;

    // Scaling
    lambda = width/108.1;

    // Initial positions
    pos_x = 97.0 * lambda;
    pos_y = -24.31 * lambda;

    // Newtonian momenta
    if (params->pn_terms[0] == 1 && params->pn_terms[1] == 0 &&
        params->pn_terms[2] == 0 && params->pn_terms[3] == 0)
    {
        px = -0.09324/sqrt(lambda);
        py = -0.08647/sqrt(lambda);
    }
    // 1PN momenta
    else if (params->pn_terms[0] == 1 && params->pn_terms[1] == 1 && params->pn_terms[2] == 0) {
        px = -sqrt(0.008693032833827606/lambda + 0.000798860400642637/pow(lambda, 2)
                                               + 0.00013381114672890315/pow(lambda, 3));
        py = -sqrt(0.007480061222224325/lambda + 0.001241927410006741/pow(lambda, 2)
                                               + 0.00028564641727617235/pow(lambda, 3));
    }
    // 2PN momenta
    else if (params->pn_terms[0] == 1 && params->pn_terms[1] == 1 && params->pn_terms[2] == 1) {
        px = -sqrt(0.008692910686038705/lambda + 0.0007977722653187864/pow(lambda, 2)
             + 6.351332174711012e-05/pow(lambda, 3) - 3.0103527312470005e-05/pow(lambda, 4));
        py = -sqrt(0.007477759360235814/lambda + 0.0012723359375445093/pow(lambda, 2)
             + 6.583447113705563e-05/pow(lambda, 3) - 5.3457362474338285e-06/pow(lambda, 4));
    }

    // Set initial parameters (scale the system by the masses in case they differ from 1.0)
    if (params->num_dim == 2) {
        w0[0] = pos_x * params->masses[0];
        w0[1] = pos_y * params->masses[0];

        w0[2] = -pos_x * params->masses[1];
        w0[3] = -pos_y * params->masses[1];

        w0[4] = 0.0;
        w0[5] = 0.0;

        w0[6] = -0.5*px * params->masses[0];
        w0[7] = -0.5*py * params->masses[0];

        w0[8] = -0.5*px * params->masses[1];
        w0[9] = -0.5*py * params->masses[1];

        w0[10] = px * params->masses[2];
        w0[11] = py * params->masses[2];
    }
    else if (params->num_dim == 3) {
        w0[0] = pos_x * params->masses[0];
        w0[1] = pos_y * params->masses[0];
        w0[2] = 0.0;

        w0[3] = -pos_x * params->masses[1];
        w0[4] = -pos_y * params->masses[1];
        w0[5] = 0.0;

        w0[6] = 0.0;
        w0[7] = 0.0;
        w0[8] = 0.0;

        w0[9] = -0.5*px * params->masses[0];
        w0[10] = -0.5*py * params->masses[0];
        w0[11] = 0.0;

        w0[12] = -0.5*px * params->masses[1];
        w0[13] = -0.5*py * params->masses[1];
        w0[14] = 0.0;

        w0[15] = px * params->masses[2];
        w0[16] = py * params->masses[2];
        w0[17] = 0.0;
    }
}


// ------------------------------------------------------------------------------------------------
// Newtonian Plummer cluster initial data
// ------------------------------------------------------------------------------------------------

/**
 * @brief Sample a radius from a Plummer distribution truncated at a maximum radius.
 *
 * Draws a uniform variate and applies the inverse Plummer cumulative-mass distribution,
 *
 *     M(<r) / M = r^3 / (r^2 + a^2)^(3/2),
 *     r = a / sqrt(u^(-2/3) - 1).
 *
 * Samples larger than rmax are rejected, producing the Plummer radial distribution conditioned on
 * r <= rmax. The uniform variate is kept away from zero and one to avoid numerical singularities.
 * The function terminates through errorexit if no valid radius is obtained after 100000 attempts.
 *
 * @param[in]     a      Plummer scale radius; expected to be positive.
 * @param[in]     rmax   Maximum accepted radius; expected to be positive.
 * @param[in,out] rng    State of the pseudo-random number generator.
 *
 * @return A sampled radius satisfying 0 < r <= rmax.
 */
static double ic_sample_truncated_plummer_radius(double a, double rmax, unsigned long long *rng)
{
    for (int attempt = 0; attempt < 100000; attempt++) {
        double u = rng_uniform(rng);

        if (u < 1e-14)
            u = 1e-14;
        if (u > 1.0 - 1e-14)
            u = 1.0 - 1e-14;

        const double denom = pow(u, -2.0 / 3.0) - 1.0;
        if (denom <= 0.0)
            continue;

        const double r = a / sqrt(denom);

        if (r <= rmax)
            return r;
    }

    errorexit("Could not sample truncated Plummer radius");
}


/**
 * @brief Sample a Plummer-model speed as a fraction of the local escape speed.
 *
 * Uses Aarseth-style rejection sampling for the dimensionless speed fraction
 *
 *     q = v / v_escape,
 *
 * whose Plummer distribution is proportional to
 *
 *     g(q) = q^2 (1 - q^2)^(7/2),    0 <= q <= 1.
 *
 * Candidate values of q are drawn uniformly from [0,1), while the comparison value is drawn
 * uniformly from [0,0.1). The value 0.1 bounds the maximum of g(q), so accepting candidates for
 * which y <= g(q) produces the required speed-fraction distribution. Sampling continues until a
 * candidate is accepted.
 *
 * @param[in,out] rng    State of the pseudo-random number generator.
 *
 * @return The sampled dimensionless speed fraction q in the range [0,1).
 */
static double ic_sample_plummer_speed_fraction(unsigned long long *rng)
{
    while (1) {
        const double q = rng_uniform(rng);
        const double y = 0.1 * rng_uniform(rng);
        const double g = q*q * pow(1.0 - q*q, 3.5);

        if (y <= g)
            return q;
    }
}


// Compute the Newtonian gravitational potential energy of a cluster.
static double ic_cluster_newtonian_potential_energy(struct ode_params* params, double* x)
{
    const int N = params->num_bodies_initial;
    const int D = params->num_dim;

    double U = 0.0;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double rij2 = 0.0;

            for (int k = 0; k < D; k++) {
                const double dx = x[i * D + k] - x[j * D + k];
                rij2 += dx * dx;
            }

            const double rij = sqrt(rij2);

            if (rij <= 0.0)
                errorexit("Zero pair separation in virialized cluster");

            U -= params->masses[i] * params->masses[j] / rij;
        }
    }

    return U;
}


/**
 * @brief Generate a finite-N Newtonian Plummer cluster.
 *
 * Samples positions from a spherically symmetric Plummer distribution truncated at a finite outer
 * radius. Realizations containing pairs closer than the requested mass-scaled minimum separation
 * are rejected. The mass-weighted center of mass is then shifted to the origin.
 *
 * Velocities are sampled isotropically from the Newtonian Plummer distribution function using the
 * local Newtonian escape speed. The center-of-mass velocity is removed, after which all velocities
 * are rescaled by a common factor so that the finite-N realization satisfies the requested
 * Newtonian virial ratio
 *
 *     K / |U| = virial_ratio,
 *
 * where
 *
 *     K = sum_i m_i v_i^2 / 2,
 *     U = -sum_{i<j} m_i m_j / r_ij.
 *
 * The cluster scale is specified through the nominal Newtonian virial radius,
 *
 *     R_v = compactness * M_tot,
 *
 * and the Plummer scale radius is obtained from the untruncated Plummer-model relation
 *
 *     a = 3 pi R_v / 16.
 *
 * The output canonical momenta are assigned using the Newtonian relation p_i = m_i v_i.
 * Enabled PN terms do not affect the construction. Consequently, this routine provides approximate
 * initial data for PN evolution but does not construct a PN or relativistic equilibrium.
 *
 * @param[in]  params           System parameters
 * @param[in]  compactness      Ratio R_v / M_tot describing the cluster's compactness.
 * @param[in]  virial_ratio     Requested Newtonian ratio K / |U|.
 * @param[in]  rmax_factor      Truncation radius in units of R_v.
 * @param[in]  min_sep_factor   Minimum pair separation in units of m_i + m_j; zero disables it.
 * @param[in]  seed             Seed for deterministic random sampling.
 * @param[out] w0               Initial state vector containing positions and momenta.
 */
void ic_newtonian_plummer_cluster(struct ode_params* params, double compactness,
    double virial_ratio, double rmax_factor, double min_sep_factor, unsigned long long seed,
    double* w0)
{
    const int N = params->num_bodies_initial;
    const int D = params->num_dim;

    if (D != 3)
        errorexit("virialized_cluster currently requires num_dim = 3");

    if (N < 3)
        errorexit("virialized_cluster requires num_bodies >= 3");

    double Mtot = 0.0;
    for (int i = 0; i < N; i++) {
        if (params->masses[i] <= 0.0)
            errorexit("virialized_cluster requires all masses to be > 0");

        Mtot += params->masses[i];
    }

    if (compactness <= 0.0)
        errorexit("cluster_compactness must be > 0");

    if (virial_ratio <= 0.0)
        errorexit("cluster_virial_ratio must be > 0");

    if (rmax_factor <= 0.0)
        errorexit("cluster_rmax_factor must be > 0");

    if (min_sep_factor < 0.0)
        errorexit("cluster_min_separation_factor must be >= 0");

    if (seed == 0ULL)
        seed = 1ULL;

    // Interpret compactness as virial radius / total mass.
    const double Rv = compactness * Mtot;
    const double a = (3.0 * M_PI / 16.0) * Rv;
    const double rmax = rmax_factor * Rv;

    double *x;
    double *v;
    allocate_vector(&x, N * D);
    allocate_vector(&v, N * D);

    unsigned long long rng = seed;

    /*
     * Sample positions. If the initial finite-N realization contains very close
     * pairs, resample the entire cluster. This avoids immediate artificial mergers.
     */
    int accepted = 0;

    for (int cluster_attempt = 0; cluster_attempt < 10000; cluster_attempt++) {
        for (int i = 0; i < N; i++) {
            double n[3];
            random_unit_vector(&rng, n);

            const double r = ic_sample_truncated_plummer_radius(a, rmax, &rng);

            for (int k = 0; k < D; k++)
                x[i * D + k] = r * n[k];
        }

        //Shift mass-weighted center of mass to the origin
        double xcom[3] = {0.0, 0.0, 0.0};

        for (int i = 0; i < N; i++) {
            for (int k = 0; k < D; k++)
                xcom[k] += params->masses[i] * x[i * D + k];
        }

        for (int k = 0; k < D; k++)
            xcom[k] /= Mtot;

        for (int i = 0; i < N; i++) {
            for (int k = 0; k < D; k++)
                x[i * D + k] -= xcom[k];
        }

        if (min_sep_ok(params, x, min_sep_factor)) {
            accepted = 1;
            break;
        }
    }

    if (!accepted)
        errorexit("Could not sample virialized cluster satisfying min-separation criterion");

    // Sample velocities from the isotropic Plummer distribution function
    for (int i = 0; i < N; i++) {
        double r2 = 0.0;
        for (int k = 0; k < D; k++)
            r2 += x[i * D + k] * x[i * D + k];

        const double r = sqrt(r2);
        const double vesc = sqrt(2.0 * Mtot / sqrt(r*r + a*a));
        const double q = ic_sample_plummer_speed_fraction(&rng);
        const double speed = q * vesc;

        double n[3];
        random_unit_vector(&rng, n);

        for (int k = 0; k < D; k++)
            v[i * D + k] = speed * n[k];
    }

    // Remove center-of-mass velocity
    double pcom[3] = {0.0, 0.0, 0.0};

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < D; k++)
            pcom[k] += params->masses[i] * v[i * D + k];
    }

    for (int k = 0; k < D; k++)
        pcom[k] /= Mtot;

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < D; k++)
            v[i * D + k] -= pcom[k];
    }

    // Rescale velocities to exact requested virial ratio K/|U|
    double K = 0.0;
    for (int i = 0; i < N; i++) {
        double v2 = 0.0;

        for (int k = 0; k < D; k++)
            v2 += v[i * D + k] * v[i * D + k];

        K += 0.5 * params->masses[i] * v2;
    }

    const double U = ic_cluster_newtonian_potential_energy(params, x);

    if (K <= 0.0)
        errorexit("Initial kinetic energy is zero in virialized cluster");

    if (U >= 0.0)
        errorexit("Initial potential energy is non-negative in virialized cluster");

    const double velocity_scale = sqrt(virial_ratio * fabs(U) / K);

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < D; k++)
            v[i * D + k] *= velocity_scale;
    }

    // Fill the state vector
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < D; k++) {
            w0[i * D + k] = x[i * D + k];
            w0[N * D + i * D + k] = params->masses[i] * v[i * D + k];
        }
    }

    printf("Virialized cluster diagnostics:\n");
    printf("Mtot                     = % .16e\n", Mtot);
    printf("cluster_compactness Rv/M = % .16e\n", compactness);
    printf("Plummer scale a          = % .16e\n", a);
    printf("Truncation radius rmax   = % .16e\n", rmax);
    printf("Initial potential U      = % .16e\n", U);
    printf("Requested K/|U|          = % .16e\n", virial_ratio);
    print_divider();

    free_vector(x);
    free_vector(v);
}


// ------------------------------------------------------------------------------------------------
// Relativistic monoenergetic cluster, close to Bamber et al. 2025 Appendix A
// ------------------------------------------------------------------------------------------------

// Compute the normalized energy density rho/rho_c of the monoenergetic model.
static double rc_rho_bar(double z, double zc)
{
    if (z <= 1.0)
        return 0.0;

    const double denom = pow(zc, 3.0) * sqrt(zc*zc - 1.0);
    return pow(z, 3.0) * sqrt(z*z - 1.0) / denom;
}


// Compute the normalized isotropic pressure P/rho_c of the monoenergetic model.
static double rc_p_bar(double z, double zc)
{
    if (z <= 1.0)
        return 0.0;

    const double denom = pow(zc, 3.0) * sqrt(zc*zc - 1.0);
    return (1.0 / 3.0) * z * pow(z*z - 1.0, 1.5) / denom;
}


// Compute the rest-mass density ratio rho_0/rho_c of the monoenergetic model.
static double rc_rho0_over_rhoc(double z, double zc)
{
    if (z <= 1.0)
        return 0.0;

    const double denom = pow(zc, 3.0) * sqrt(zc*zc - 1.0);
    return z*z * sqrt(z*z - 1.0) / denom;
}


// Append one radial solution point to the model's preallocated profile arrays.
static int rc_append_model_point(struct rel_mono_model *m, double x, double mu, double z)
{
    if (m->n >= m->cap)
        return 0;

    m->x[m->n] = x;
    m->mu[m->n] = mu;
    m->z[m->n] = z;
    m->n++;

    return 1;
}


// Evaluate the dimensionless relativistic equilibrium equations at one radius.
static void rc_rhs(double x, double mu, double z, double zc, double *dmu_dx, double *dz_dx,
    double *dmu0_dx)
{
    const double rho = rc_rho_bar(z, zc);
    const double P = rc_p_bar(z, zc);
    const double rho0 = rc_rho0_over_rhoc(z, zc);

    *dmu_dx = x*x * rho;

    if (x <= 0.0) {
        *dz_dx = 0.0;
        *dmu0_dx = 0.0;
        return;
    }

    const double A = 1.0 - 2.0 * mu / x;

    if (A <= 0.0)
        errorexit("relativistic monoenergetic model crossed x = 2 mu");

    *dz_dx = -z * (mu + x*x*x * P) / (x * (x - 2.0 * mu));

    // Proper rest-mass integral
    *dmu0_dx = x*x * rho0 / sqrt(A);
}


// Advance the relativistic equilibrium solution by one fourth-order Runge–Kutta step.
static void rc_rk4_step(double x, double h, double *mu, double *z, double *mu0, double zc)
{
    double k1m, k1z, k1m0;
    double k2m, k2z, k2m0;
    double k3m, k3z, k3m0;
    double k4m, k4z, k4m0;

    rc_rhs(x, *mu, *z, zc, &k1m, &k1z, &k1m0);
    rc_rhs(x + 0.5*h, *mu + 0.5*h*k1m, *z + 0.5*h*k1z, zc,
        &k2m, &k2z, &k2m0);
    rc_rhs(x + 0.5*h, *mu + 0.5*h*k2m, *z + 0.5*h*k2z, zc,
        &k3m, &k3z, &k3m0);
    rc_rhs(x + h, *mu + h*k3m, *z + h*k3z, zc,
        &k4m, &k4z, &k4m0);

    *mu  += h * (k1m  + 2.0*k2m  + 2.0*k3m  + k4m)  / 6.0;
    *z   += h * (k1z  + 2.0*k2z  + 2.0*k3z  + k4z)  / 6.0;
    *mu0 += h * (k1m0 + 2.0*k2m0 + 2.0*k3m0 + k4m0) / 6.0;
}


// Release the dynamically allocated profile arrays owned by a relativistic cluster model
static void rc_free_model(struct rel_mono_model *m)
{
    free_vector(m->x);
    free_vector(m->mu);
    free_vector(m->z);
    free_vector(m->xiso);
    free_vector(m->cdf);

    m->x = NULL;
    m->mu = NULL;
    m->z = NULL;
    m->xiso = NULL;
    m->cdf = NULL;
}


/**
 * @brief Integrate a dimensionless relativistic monoenergetic cluster model.
 *
 * Constructs the smooth, spherically symmetric Shapiro-Teukolsky/Bamber monoenergetic equilibrium
 * associated with the central local Lorentz factor zc. The dimensionless variables are
 *
 *     x  = r sqrt(4 pi rho_c),
 *     mu = m sqrt(4 pi rho_c),
 *
 * where rho_c is the central energy density. The local Lorentz factor z decreases outward from zc,
 * and the matter surface is defined by z = 1, where the density and pressure vanish.
 *
 * The function integrates the coupled equations for the enclosed gravitational mass mu,
 * local Lorentz factor z, and proper rest mass mu0. Integration starts close to the regular center
 * using the leading-order behavior
 *
 *     mu(x) = x^3 / 3,
 *
 * which follows from the normalization rho(zc) / rho_c = 1. A fourth-order Runge-Kutta method is
 * used with a radial step that grows with radius but is reduced near the surface to limit the
 * fractional change in z. When a step crosses z = 1, the final solution point is linearly
 * interpolated to the surface.
 *
 * A solution is rejected if it becomes nonfinite, exhausts the available grid before reaching the
 * surface, contains too few radial points, or approaches the condition x <= 2 mu. For a successful
 * solution, the surface quantities x_surf, mu_surf, mu0_surf and the inverse compactness
 * R / M = x_surf / mu_surf are stored in the output model.
 *
 * After the equilibrium integration, the areal-radius profile is converted to an isotropic-radius
 * profile. The exterior Schwarzschild relation fixes the isotropic radius at the surface, and the
 * coordinate transformation
 *
 *     d ln(xiso) / dx = 1 / [x sqrt(1 - 2 mu/x)]
 *
 * is integrated inward through the numerical mass profile.
 *
 * Finally, a radial sampling CDF is constructed from the proper rest-mass element
 *
 *     dmu0 = x^2 (rho0/rho_c) / sqrt(1 - 2 mu/x) dx.
 *
 * The CDF is evaluated by trapezoidal integration and normalized to the range [0,1]. It is
 * subsequently used to sample particle radii according to proper rest mass rather than coordinate
 * volume or gravitational mass.
 *
 * @param[in]  zc      Central local Lorentz factor; must be greater than one.
 * @param[in]  ngrid   Maximum number of radial profile points to allocate and
 *                     store.
 * @param[out] m       Model structure receiving the integrated profiles,
 *                     surface quantities, compactness, and normalized CDF.
 *
 * @return 1 if the equilibrium model, isotropic-radius mapping, and sampling CDF are successful;
 *         0 if the input or numerical solution was invalid.
 */
static int rc_integrate_model(double zc, int ngrid, struct rel_mono_model *m)
{
    // Dimensionless integration of the monoenergetic spherical equilibrium sequence
    if (zc <= 1.0)
        return 0;

    m->n = 0;
    m->cap = ngrid;

    allocate_vector(&m->x, ngrid);
    allocate_vector(&m->mu, ngrid);
    allocate_vector(&m->z, ngrid);
    allocate_vector(&m->xiso, ngrid);
    allocate_vector(&m->cdf, ngrid);

    m->zc = zc;

    // Start close to the origin. Since rho_bar(zc)=1 by construction, mu ~ x^3 / 3
    double x = 1e-6;
    double mu = x*x*x / 3.0;
    double z = zc;
    double mu0 = x*x*x * rc_rho0_over_rhoc(zc, zc) / 3.0;

    if (!rc_append_model_point(m, x, mu, z)) {
        rc_free_model(m);
        return 0;
    }

    while (z > 1.0 && m->n < ngrid) {
        const double x_old = x;
        const double mu_old = mu;
        const double z_old = z;
        const double mu0_old = mu0;

        double h = 2e-3 * fmax(1.0, x);

        // Limit fractional change in z per step
        double dm_tmp, dz_tmp, dm0_tmp;
        rc_rhs(x, mu, z, zc, &dm_tmp, &dz_tmp, &dm0_tmp);

        if (fabs(dz_tmp) > 0.0) {
            const double h_z = 0.01 * fabs(z - 1.0) / fabs(dz_tmp);
            if (h_z > 0.0 && h_z < h)
                h = h_z;
        }

        if (h < 1e-7)
            h = 1e-7;

        rc_rk4_step(x, h, &mu, &z, &mu0, zc);
        x += h;

        if (!isfinite(x) || !isfinite(mu) || !isfinite(z) || !isfinite(mu0)) {
            rc_free_model(m);
            return 0;
        }

        if (x <= 2.0 * mu) {
            rc_free_model(m);
            return 0;
        }

        if (z <= 1.0) {
            // Interpolate last step to surface z = 1
            const double a = (z_old - 1.0) / (z_old - z);

            const double x_s = x_old + a * (x - x_old);
            const double mu_s = mu_old + a * (mu - mu_old);
            const double mu0_s = mu0_old + a * (mu0 - mu0_old);

            if (!rc_append_model_point(m, x_s, mu_s, 1.0)) {
                rc_free_model(m);
                return 0;
            }

            m->x_surf = x_s;
            m->mu_surf = mu_s;
            m->mu0_surf = mu0_s;
            m->compactness = x_s / mu_s;
            break;
        }

        if (!rc_append_model_point(m, x, mu, z)) {
            rc_free_model(m);
            return 0;
        }
    }

    if (m->n < 10 || m->z[m->n - 1] > 1.0) {
        rc_free_model(m);
        return 0;
    }

    if (m->mu_surf <= 0.0 || m->x_surf <= 2.0 * m->mu_surf) {
        rc_free_model(m);
        return 0;
    }

    // Isotropic-radius mapping
    const int n = m->n;
    const double Xs = m->x_surf;
    const double Ms = m->mu_surf;

    const double xiso_s =
        0.5 * (Xs - Ms + sqrt(Xs * (Xs - 2.0 * Ms)));

    m->xiso[n - 1] = xiso_s;

    double g_next = log(xiso_s / Xs);

    for (int i = n - 2; i >= 0; i--) {
        const double x1 = m->x[i];
        const double x2 = m->x[i + 1];

        const double A1 = 1.0 - 2.0 * m->mu[i] / x1;
        const double A2 = 1.0 - 2.0 * m->mu[i + 1] / x2;

        if (A1 <= 0.0 || A2 <= 0.0) {
            rc_free_model(m);
            return 0;
        }

        const double f1 = (1.0 / sqrt(A1) - 1.0) / x1;
        const double f2 = (1.0 / sqrt(A2) - 1.0) / x2;

        const double dg = 0.5 * (f1 + f2) * (x2 - x1);
        const double g_i = g_next - dg;

        m->xiso[i] = x1 * exp(g_i);
        g_next = g_i;
    }

    // CDF from proper rest-mass density
    m->cdf[0] = 0.0;

    for (int i = 1; i < n; i++) {
        const double x1 = m->x[i - 1];
        const double x2 = m->x[i];

        const double A1 = 1.0 - 2.0 * m->mu[i - 1] / x1;
        const double A2 = 1.0 - 2.0 * m->mu[i] / x2;

        if (A1 <= 0.0 || A2 <= 0.0) {
            rc_free_model(m);
            return 0;
        }

        const double w1 = x1*x1 * rc_rho0_over_rhoc(m->z[i - 1], zc) / sqrt(A1);
        const double w2 = x2*x2 * rc_rho0_over_rhoc(m->z[i], zc) / sqrt(A2);

        m->cdf[i] = m->cdf[i - 1] + 0.5 * (w1 + w2) * (x2 - x1);
    }

    const double total = m->cdf[n - 1];

    if (total <= 0.0 || !isfinite(total)) {
        rc_free_model(m);
        return 0;
    }

    for (int i = 0; i < n; i++)
        m->cdf[i] /= total;

    return 1;
}


// Compute the model radius-to-mass ratio R/M for central Lorentz factor zc.
static double rc_model_compactness(double zc, int ngrid)
{
    struct rel_mono_model m;

    if (!rc_integrate_model(zc, ngrid, &m))
        return -1.0;

    const double C = m.compactness;
    rc_free_model(&m);
    return C;
}


// Return whether two valid compactness values bracket the requested target.
static int rc_compactness_bracketed(double C1, double C2, double target)
{
    const double f1 = C1 - target;
    const double f2 = C2 - target;

    return (f1 <= 0.0 && f2 >= 0.0) ||
           (f2 <= 0.0 && f1 >= 0.0);
}


/**
 * @brief Refine a valid central-Lorentz-factor bracket by bisection.
 *
 * Both endpoints must correspond to successful model integrations and their compactness values
 * must bracket the requested target. Every midpoint is evaluated before it replaces an endpoint,
 * so a failed integration never becomes part of the active bracket. The refinement succeeds only
 * when the returned model matches the requested compactness within the supplied tolerance.
 *
 * @return 1 on success and 0 if the bracket or any midpoint is invalid.
 */
static int rc_refine_zc_bracket(double zc_lo, double C_lo, double zc_hi, double C_hi,
    double target_compactness, int ngrid, double compactness_tol, double *zc_solution)
{
    const double zc_rtol = 1e-12;

    if (zc_solution == NULL || zc_lo >= zc_hi ||
        C_lo <= 0.0 || !isfinite(C_lo) ||
        C_hi <= 0.0 || !isfinite(C_hi) ||
        !rc_compactness_bracketed(C_lo, C_hi, target_compactness))
        return 0;

    if (fabs(C_lo - target_compactness) <= compactness_tol) {
        *zc_solution = zc_lo;
        return 1;
    }

    if (fabs(C_hi - target_compactness) <= compactness_tol) {
        *zc_solution = zc_hi;
        return 1;
    }

    for (int it = 0; it < 80; it++) {
        const double zc_mid = zc_lo + 0.5 * (zc_hi - zc_lo);

        if (zc_mid == zc_lo || zc_mid == zc_hi)
            break;

        const double C_mid = rc_model_compactness(zc_mid, ngrid);

        if (C_mid <= 0.0 || !isfinite(C_mid))
            return 0;

        if (fabs(C_mid - target_compactness) <= compactness_tol) {
            *zc_solution = zc_mid;
            return 1;
        }

        if (rc_compactness_bracketed(C_lo, C_mid, target_compactness)) {
            zc_hi = zc_mid;
            C_hi = C_mid;
        }
        else {
            zc_lo = zc_mid;
            C_lo = C_mid;
        }

        if (zc_hi - zc_lo <=
            zc_rtol * fmax(1.0, fabs(zc_lo + 0.5 * (zc_hi - zc_lo))))
            break;
    }

    if (fabs(C_lo - target_compactness) <=
        fabs(C_hi - target_compactness)) {
        if (fabs(C_lo - target_compactness) > compactness_tol)
            return 0;

        *zc_solution = zc_lo;
    }
    else {
        if (fabs(C_hi - target_compactness) > compactness_tol)
            return 0;

        *zc_solution = zc_hi;
    }

    return 1;
}

/**
 * @brief Find the central Lorentz factor that produces a requested inverse compactness.
 *
 * The solver scans the central Lorentz factor `zc` logarithmically over a fixed range above unity
 * and evaluates the corresponding model radius-to-mass ratio `R/M`. Failed or non-finite model
 * integrations are skipped. Once two successive valid models bracket the requested value, the
 * solution is refined with 80 bisection iterations. A failed midpoint evaluation contracts the
 * search toward the valid lower endpoint. Valid scan results are printed for diagnostic purposes,
 * and failure to find a bracket terminates the program.
 *
 * @param   target_compactness  Desired inverse compactness `R/M`.
 * @param   ngrid               Maximum number of radial grid points used to integrate each model.
 * @return Central Lorentz factor `zc` corresponding to the requested `R/M`.
 */
static double rc_solve_zc_for_compactness(double target_compactness, int ngrid)
{
    const double zc_min = 1.0005;
    const double zc_max = 50.0;
    const int nscan = 400;

    if (target_compactness <= 0.0 || !isfinite(target_compactness))
        errorexit("cluster_compactness must be finite and positive");

    if (ngrid < 10)
        errorexit("cluster_ngrid must be at least 10");

    const double compactness_tol =
        1e-8 * fmax(1.0, fabs(target_compactness));

    const double log_zc_min = log(zc_min);
    const double log_zc_max = log(zc_max);

    int have_prev = 0;
    double zc_prev = 0.0;
    double C_prev = 0.0;

    for (int s = 0; s <= nscan; s++) {
        const double fraction = (double)s / (double)nscan;
        const double zc = exp(log_zc_min + fraction * (log_zc_max - log_zc_min));

        const double C = rc_model_compactness(zc, ngrid);

        if (C <= 0.0 || !isfinite(C)) {
            have_prev = 0;
            continue;
        }

        printf("rel. cluster scan: zc = %.12e, R/M = %.12e\n", zc, C);

        if (fabs(C - target_compactness) <= compactness_tol)
            return zc;

        if (have_prev &&
            rc_compactness_bracketed(C_prev, C, target_compactness)) {
            double zc_solution;

            if (rc_refine_zc_bracket(
                    zc_prev, C_prev, zc, C, target_compactness,
                    ngrid, compactness_tol, &zc_solution))
                return zc_solution;
        }

        zc_prev = zc;
        C_prev = C;
        have_prev = 1;
    }

    errorexit("Could not find a valid zc for requested cluster_compactness. "
              "Try increasing cluster_ngrid, using a less extreme compactness, "
              "or specifying cluster_central_z manually.");

    return -1.0;
}


// Sample proper-rest-mass CDF and interpolate areal radius, isotropic radius and Lorentz factor.
static void rc_sample_radius(struct rel_mono_model *m, unsigned long long *rng,
    double *x_areal, double *x_iso, double *z_local)
{
    const double u = rng_uniform(rng);

    int lo = 0;
    int hi = m->n - 1;

    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;

        if (m->cdf[mid] < u)
            lo = mid;
        else
            hi = mid;
    }

    double a = 0.0;

    if (m->cdf[hi] > m->cdf[lo])
        a = (u - m->cdf[lo]) / (m->cdf[hi] - m->cdf[lo]);

    *x_areal = m->x[lo] + a * (m->x[hi] - m->x[lo]);
    *x_iso = m->xiso[lo] + a * (m->xiso[hi] - m->xiso[lo]);
    *z_local = m->z[lo] + a * (m->z[hi] - m->z[lo]);
}


/**
 * @brief Generate initial data for a relativistic monoenergetic spherical cluster.
 *
 * Constructs the smooth equilibrium model for a specified central redshift parameter, scales it to
 * the total mass and independently samples isotropic positions and momenta for all bodies. If
 * solve_central_z is nonzero, the central redshift parameter is determined from
 * target_compactness; otherwise, central_z is used directly.
 *
 * Samples are rejected until every body pair satisfies
 *
 *     r_ij >= min_sep_factor * (m_i + m_j),
 *
 * or until the maximum number of sampling attempts is reached. Optionally, the finite-N
 * center-of-mass position and total momentum are removed. The resulting state is stored in w0.
 *
 * The underlying equilibrium model assumes equal-mass particles.
 *
 * @param[in]  params              Simulation parameters.
 * @param[in]  target_compactness  Target smooth-model compactness R/M.
 * @param[in]  central_z           Central redshift parameter.
 * @param[in]  solve_central_z     If nonzero, solve for the central redshift.
 * @param[in]  ngrid               Number of grid points used to integrate the equilibrium model.
 * @param[in]  seed                Seed for random sampling.
 * @param[in]  min_sep_factor      Minimum pair-separation factor.
 * @param[in]  remove_com          If nonzero, remove the center-of-mass displacement and momentum.
 * @param[out] w0                  Initial state array containing 2 * N * D elements.
 */
void ic_relativistic_monoenergetic_cluster(struct ode_params* params, double target_compactness,
    double central_z, int solve_central_z, int ngrid, unsigned long long seed,
    double min_sep_factor, int remove_com, double* w0)
{
    const int N = params->num_bodies_initial;
    const int D = params->num_dim;

    if (D != 3)
        errorexit("relativistic_monoenergetic_cluster requires num_dim = 3");

    double M0_target = 0.0;

    for (int i = 0; i < N; i++) {
        if (params->masses[i] <= 0.0)
            errorexit("relativistic_monoenergetic_cluster requires all masses > 0");

        M0_target += params->masses[i];
    }

    for (int i = 1; i < N; i++) {
        if (fabs(params->masses[i] - params->masses[0]) > 1e-12 * M0_target) {
            printf("Warning: the Bamber et al. monoenergetic cluster construction "
                   "assumes equal rest-mass particles. Unequal masses are only a "
                   "pragmatic PN extension.\n");
            break;
        }
    }

    if (seed == 0ULL)
        seed = 1ULL;

    double zc = central_z;

    if (solve_central_z)
        zc = rc_solve_zc_for_compactness(target_compactness, ngrid);

    if (zc <= 1.0)
        errorexit("central z must be > 1");

    struct rel_mono_model model;

    if (!rc_integrate_model(zc, ngrid, &model))
        errorexit("Could not integrate relativistic monoenergetic cluster model");

    // Scale dimensionless model to the total mass used by the PN code
    const double length_scale = M0_target / model.mu0_surf;
    const double M_smooth = length_scale * model.mu_surf;

    double *x;
    double *p;

    allocate_vector(&x, N * D);
    allocate_vector(&p, N * D);

    unsigned long long rng = seed;
    int accepted = 0;

    for (int attempt = 0; attempt < 10000; attempt++) {
        for (int i = 0; i < N; i++) {
            double n_pos[3];
            double n_p[3];

            random_unit_vector(&rng, n_pos);
            random_unit_vector(&rng, n_p);

            double x_areal, x_iso, z_local;
            rc_sample_radius(&model, &rng, &x_areal, &x_iso, &z_local);

            const double r_areal = length_scale * x_areal;
            const double r_iso = length_scale * x_iso;

            for (int k = 0; k < D; k++)
                x[i*D + k] = r_iso * n_pos[k];

            // Local orthonormal momentum magnitude
            double p_mag = params->masses[i] *
                sqrt(fmax(0.0, z_local*z_local - 1.0));

            if (r_iso <= 0.0)
                errorexit("zero isotropic radius in relativistic cluster");

            const double psi2 = r_areal / r_iso;
            p_mag *= psi2;

            for (int k = 0; k < D; k++)
                p[i*D + k] = p_mag * n_p[k];
        }

        if (min_sep_ok(params, x, min_sep_factor)) {
            accepted = 1;
            break;
        }
    }

    if (!accepted)
        errorexit("Could not sample relativistic monoenergetic cluster satisfying "
                  "cluster_min_separation_factor");

    if (remove_com) {
        double xcom[3] = {0.0, 0.0, 0.0};
        double ptot[3] = {0.0, 0.0, 0.0};

        for (int i = 0; i < N; i++) {
            for (int k = 0; k < D; k++) {
                xcom[k] += params->masses[i] * x[i*D + k];
                ptot[k] += p[i*D + k];
            }
        }

        for (int k = 0; k < D; k++)
            xcom[k] /= M0_target;

        for (int i = 0; i < N; i++) {
            for (int k = 0; k < D; k++) {
                x[i*D + k] -= xcom[k];
                p[i*D + k] -= (params->masses[i] / M0_target) * ptot[k];
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < D; k++) {
            w0[i*D + k] = x[i*D + k];
            w0[N*D + i*D + k] = p[i*D + k];
        }
    }

    printf("Relativistic monoenergetic cluster diagnostics:\n");
    printf("central_z zc                 = % .16e\n", zc);
    printf("smooth R/M                   = % .16e\n", model.compactness);
    printf("smooth M0/M                  = % .16e\n", model.mu0_surf / model.mu_surf);
    printf("dimensionless R              = % .16e\n", model.x_surf);
    printf("dimensionless M              = % .16e\n", model.mu_surf);
    printf("dimensionless M0             = % .16e\n", model.mu0_surf);
    printf("physical M0 scale            = % .16e\n", M0_target);
    printf("physical gravitational M     = % .16e\n", M_smooth);
    printf("physical areal R             = % .16e\n", length_scale * model.x_surf);
    printf("physical isotropic Rbar      = % .16e\n", length_scale * model.xiso[model.n - 1]);
    printf("finite-N COM removed         = %d\n", remove_com);
    print_divider();

    free_vector(x);
    free_vector(p);
    rc_free_model(&model);
}
