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


// ------------------------------------------------------------------------------------------------
// PN initial configuration helper functions
// ------------------------------------------------------------------------------------------------


static double pn_hhat_rel(double x, double pr_hat, double j, double nu, int use_1pn, int use_2pn)
{
    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;
    const double inv_x3 = inv_x2 * inv_x;

    const double pr2 = pr_hat * pr_hat;
    const double j2 = j * j;
    const double p2 = pr2 + j2 * inv_x2;

    const double p4 = p2 * p2;
    const double p6 = p4 * p2;

    double h = 0.5 * p2 - inv_x;

    if (use_1pn) {
        h += 0.125 * (3.0 * nu - 1.0) * p4
           - 0.5 * inv_x * ((3.0 + nu) * p2 + nu * pr2)
           + 0.5 * inv_x2;
    }

    if (use_2pn) {
        h += (1.0 - 5.0 * nu + 5.0 * nu * nu) * p6 / 16.0

           + inv_x / 8.0 *
             ((5.0 - 20.0 * nu - 3.0 * nu * nu) * p4
              - 2.0 * nu * nu * pr2 * p2
              - 3.0 * nu * nu * pr2 * pr2)

           + 0.5 * inv_x2 *
             ((5.0 + 8.0 * nu) * p2 + 3.0 * nu * pr2)

           - 0.25 * (1.0 + 3.0 * nu) * inv_x3;
    }

    return h;
}


static double pn_hhat_turning_point(double x, double j, double nu, int use_1pn, int use_2pn)
{
    return pn_hhat_rel(x, 0.0, j, nu, use_1pn, use_2pn);
}


static double pn_dh_dx_turning_numeric(double x, double j, double nu, int use_1pn, int use_2pn)
{
    double dx = 1e-6 * fabs(x);

    if (dx < 1e-8) dx = 1e-8;
    if (x - dx <= 0.0) dx = 0.5 * x;

    const double hp = pn_hhat_turning_point(x + dx, j, nu, use_1pn, use_2pn);
    const double hm = pn_hhat_turning_point(x - dx, j, nu, use_1pn, use_2pn);

    return (hp - hm) / (2.0 * dx);
}


static double pn_solve_circular_j(double x, double nu, int use_1pn, int use_2pn)
{
    double lo = 0.0;
    double hi = sqrt(x) * 2.0 + 1.0;

    double flo = pn_dh_dx_turning_numeric(x, lo, nu, use_1pn, use_2pn);
    double fhi = pn_dh_dx_turning_numeric(x, hi, nu, use_1pn, use_2pn);

    int expand_count = 0;
    while (flo * fhi > 0.0 && expand_count < 80) {
        hi *= 2.0;
        fhi = pn_dh_dx_turning_numeric(x, hi, nu, use_1pn, use_2pn);
        expand_count++;
    }

    if (flo * fhi > 0.0)
        errorexit("could not bracket circular angular momentum");

    for (int i = 0; i < 160; i++) {
        const double mid = 0.5 * (lo + hi);
        const double fmid = pn_dh_dx_turning_numeric(x, mid, nu, use_1pn, use_2pn);

        if (flo * fmid <= 0.0) {
            hi = mid;
            fhi = fmid;
        }
        else {
            lo = mid;
            flo = fmid;
        }
    }

    return 0.5 * (lo + hi);
}


static double pn_solve_eccentric_j(double xp, double xa, double nu, int use_1pn, int use_2pn)
{
    double lo = 0.0;
    double hi = sqrt(0.5 * (xp + xa)) * 2.0 + 1.0;

    double flo = pn_hhat_turning_point(xp, lo, nu, use_1pn, use_2pn)
               - pn_hhat_turning_point(xa, lo, nu, use_1pn, use_2pn);

    double fhi = pn_hhat_turning_point(xp, hi, nu, use_1pn, use_2pn)
               - pn_hhat_turning_point(xa, hi, nu, use_1pn, use_2pn);

    int expand_count = 0;
    while (flo * fhi > 0.0 && expand_count < 80) {
        hi *= 2.0;
        fhi = pn_hhat_turning_point(xp, hi, nu, use_1pn, use_2pn)
            - pn_hhat_turning_point(xa, hi, nu, use_1pn, use_2pn);
        expand_count++;
    }

    if (flo * fhi > 0.0)
        errorexit("could not bracket eccentric angular momentum");

    for (int i = 0; i < 160; i++) {
        const double mid = 0.5 * (lo + hi);
        const double fmid = pn_hhat_turning_point(xp, mid, nu, use_1pn, use_2pn)
                          - pn_hhat_turning_point(xa, mid, nu, use_1pn, use_2pn);

        if (flo * fmid <= 0.0) {
            hi = mid;
            fhi = fmid;
        }
        else {
            lo = mid;
            flo = fmid;
        }
    }

    return 0.5 * (lo + hi);
}


static double pn_solve_pr_hat_abs(double x, double j, double energy, double nu, 
    int use_1pn, int use_2pn)
{
    /*
     * Solve H(x, pr, j) = energy for |pr| at fixed x and j.
     * The Hamiltonian is even in pr, so solve in q = pr^2.
     */

    double qlo = 0.0;
    double flo = pn_hhat_rel(x, 0.0, j, nu, use_1pn, use_2pn) - energy;

    if (fabs(flo) < 1e-13)
        return 0.0;

    if (flo > 0.0) {
        if (flo < 1e-10)
            return 0.0;
        errorexit("requested phase gives no real radial momentum");
    }

    double qhi = 1e-12;
    double fhi = pn_hhat_rel(x, sqrt(qhi), j, nu, use_1pn, use_2pn) - energy;

    int expand_count = 0;
    while (fhi <= 0.0 && expand_count < 120) {
        qhi *= 2.0;
        fhi = pn_hhat_rel(x, sqrt(qhi), j, nu, use_1pn, use_2pn) - energy;
        expand_count++;
    }

    if (fhi <= 0.0)
        errorexit("could not bracket radial momentum");

    for (int i = 0; i < 160; i++) {
        const double qmid = 0.5 * (qlo + qhi);
        const double fmid = pn_hhat_rel(x, sqrt(qmid), j, nu, use_1pn, use_2pn) - energy;

        if (fmid <= 0.0) {
            qlo = qmid;
            flo = fmid;
        }
        else {
            qhi = qmid;
            fhi = fmid;
        }
    }

    return sqrt(0.5 * (qlo + qhi));
}


static double pn_dh_dpr_coeff_at_zero(double x, double j, double nu, int use_1pn, int use_2pn)
{
    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;

    const double j2 = j * j;
    const double u = j2 * inv_x2;  // tangential p^2 at pr = 0

    /*
     * We compute A such that, near pr = 0,
     *
     *      dH/dpr = A * pr + O(pr^3).
     *
     * Newtonian:
     *      A_N = 1.
     */
    double A = 1.0;

    if (use_1pn) {
        /*
         * 1PN contribution:
         *
         * H_1PN =
         *      1/8 (3 nu - 1) p^4
         *    - 1/2 /x [ (3 + nu) p^2 + nu pr^2 ]
         *    + 1/2 /x^2
         *
         * with p^2 = pr^2 + j^2/x^2.
         */
        A += 0.5 * (3.0 * nu - 1.0) * u
           - (3.0 + 2.0 * nu) * inv_x;
    }

    if (use_2pn) {
        /*
         * 2PN contribution to the linear-in-pr coefficient.
         *
         * H_2PN =
         *      c6 p^6
         *    + 1/x/8 [ a p^4 - 2 nu^2 pr^2 p^2 - 3 nu^2 pr^4 ]
         *    + 1/2/x^2 [ (5 + 8nu) p^2 + 3nu pr^2 ]
         *    - 1/4(1+3nu)/x^3
         *
         * At pr = 0, only terms linear in pr in dH/dpr contribute.
         */
        const double c6 = (1.0 - 5.0 * nu + 5.0 * nu * nu) / 16.0;
        const double a4 = 5.0 - 20.0 * nu - 3.0 * nu * nu;

        A += 6.0 * c6 * u * u
           + inv_x / 8.0 * (4.0 * a4 * u - 4.0 * nu * nu * u)
           + inv_x2 * (5.0 + 11.0 * nu);
    }

    return A;
}


// ------------------------------------------------------------------------------------------------
// Initial configurations
// ------------------------------------------------------------------------------------------------


/**
 * @brief Computes initial positions and momenta for a binary.
 * 
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, p1, p2]) of two
 * bodies with masses m1 and m2 in a binary system with relative semi-major axis a, 
 * eccentricity e, and initial phase phi0 in the center of mass frame and in the xy-plane.
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
    const double phi0 = binary_params->phi0;
    const double p = binary_params->p;
    const int use_1pn = (ode_params->pn_terms[1] == 1);
    const int use_2pn = (ode_params->pn_terms[2] == 1);

    // Compute total mass, reduced mass and symmetric mass ratio
    const double M = m1 + m2; 
    const double mu = m1 * m2 / M;
    const double nu = mu / M;

    // Geometric quantities
    const double cosphi0 = cos(phi0);
    const double sinphi0 = sin(phi0);

    const double denom = 1.0 + e * cosphi0;
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
    if (e < 1e-12) 
        j = pn_solve_circular_j(a / M, nu, use_1pn, use_2pn);
    else
        j = pn_solve_eccentric_j(xp, xa, nu, use_1pn, use_2pn);

    const double energy = (e < 1e-12)
        ? pn_hhat_turning_point(a / M, j, nu, use_1pn, use_2pn)
        : pn_hhat_turning_point(xp, j, nu, use_1pn, use_2pn);

    double pr_hat_abs = pn_solve_pr_hat_abs(x0, j, energy, nu, use_1pn, use_2pn);
    double pr_hat;

    if (e < 1e-12 && ode_params->pn_terms[3] == 1) {
        /*
        * Quasi-circular inspiral radial momentum.
        *
        * Leading quadrupole circular inspiral:
        *
        *      dx/d(t/M) = -64/5 * nu / x^3
        *
        * Hamilton equation:
        *
        *      dx/d(t/M) = dHhat/dpr_hat
        *                = A * pr_hat
        *
        * so
        *
        *      pr_hat = xdot / A.
        */
        const double xdot = -(64.0 / 5.0) * nu / pow(x0, 3);

        const double A = pn_dh_dpr_coeff_at_zero(
            x0,
            j,
            nu,
            use_1pn,
            use_2pn
        );

        if (fabs(A) < 1e-14)
            errorexit("could not compute circular 2.5PN radial momentum");

        pr_hat = xdot / A;
    }
    else {
        /*
        * Conservative eccentric-orbit radial momentum.
        */
        double pr_sign = 0.0;

        if (sinphi0 > 1e-14)
            pr_sign = 1.0;
        else if (sinphi0 < -1e-14)
            pr_sign = -1.0;

        pr_hat = pr_sign * pr_hat_abs;
    }
    const double pt_hat = j / x0;

    /*
     * Unit vectors in the orbital plane.
     * n points radially outward; lambda is the positive-phi tangential unit vector.
     */
    const double nx = cosphi0;
    const double ny = sinphi0;
    const double lx = -sinphi0;
    const double ly = cosphi0;

    const double rvec_x = r0 * nx;
    const double rvec_y = r0 * ny;

    /*
     * Physical relative canonical momentum P_rel = mu * p_hat.
     * Body momenta in the COM frame are p1 = +P_rel, p2 = -P_rel.
     */
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
 * @brief Changes the position and orientation of a binary.
 * 
 * Changes the center-of-mass position and orientation of a binary, assuming the initial
 * orientation is in the z-direction (i.e. given by the vector [0, 0, 1]).
 * 
 * @param[in]       com_pos     Target center-of-mass position
 * @param[in]       orientation Target orientation vector (does not have to be normalized)
 * @param[in,out]   w0          Initial and target positions and momenta, w = [pos1, pos2, p1, p2]
 */
static void position_binary(double com_pos[3], double orientation[3], double w0[12]) 
{
    double initial_orientation[3] = {0, 0, 1};
    double pos1[3] = {w0[0], w0[1], w0[2]};
    double pos2[3] = {w0[3], w0[4], w0[5]};
    double p1[3] = {w0[6], w0[7], w0[8]};
    double p2[3] = {w0[9], w0[10], w0[11]};
    double pos1_new[3], pos2_new[3], p1_new[3], p2_new[3];
    double axis[3];
    double angle;
    double rotation_matrix[3][3];

    normalize(orientation, orientation);

    // Special case if the specified orientation equals the initial orientation
    if(orientation[0] == 0 && orientation[1] == 0 && orientation[2] == 1) {
        axis[0] = 0;
        axis[1] = 0;
        axis[2] = 1;
    }
    else {
        // Compute axis of rotation (cross product of initial and target orientation)
        cross_product(initial_orientation, orientation, axis);
    }
    
    // Compute angle of rotation (dot product) and create the rotation matrix
    angle = acos(dot_product(initial_orientation, orientation, 3));
    create_rotation_matrix(axis, angle, rotation_matrix);

    // Rotate positions and momenta of both masses
    rotate_vector(pos1, rotation_matrix, pos1_new);
    rotate_vector(pos2, rotation_matrix, pos2_new);
    rotate_vector(p1, rotation_matrix, p1_new);
    rotate_vector(p2, rotation_matrix, p2_new);

    // Update w0 with rotated and shifted positions and momenta
    w0[0] = pos1_new[0] + com_pos[0];
    w0[1] = pos1_new[1] + com_pos[1];
    w0[2] = pos1_new[2] + com_pos[2];
    w0[3] = pos2_new[0] + com_pos[0];
    w0[4] = pos2_new[1] + com_pos[1];
    w0[5] = pos2_new[2] + com_pos[2];
    w0[6] = p1_new[0];
    w0[7] = p1_new[1];
    w0[8] = p1_new[2];
    w0[9] = p2_new[0];
    w0[10] = p2_new[1];
    w0[11] = p2_new[2];
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
 * @param[in]   inner_binary_orientation  Orientation of inner binary angular momentum
 *                                        NULL -> [0, 0, 1]
 * @param[in]   outer_binary_orientation  Orientation of outer binary angular momentum
 *                                        NULL -> [0, 0, 1]
 * @param[out]  w0                        Initial positions and momenta
 */
void ic_hierarchical_triple(struct ode_params* ode_params, 
    struct binary_params* inner_binary_params, struct binary_params* outer_binary_params,
    double* inner_binary_orientation, double* outer_binary_orientation, double* w0)
{
    double m1 = ode_params->masses[0];
    double m2 = ode_params->masses[1];
    double m3 = ode_params->masses[2];
    double m_inner = m1 + m2;

    double w_inner[12];
    double w_outer[12];

    if (inner_binary_orientation == NULL) {
        inner_binary_orientation[0] = 0.0;
        inner_binary_orientation[1] = 0.0;
        inner_binary_orientation[2] = 1.0;
    }

    if (outer_binary_orientation == NULL) {
        outer_binary_orientation[0] = 0.0;
        outer_binary_orientation[1] = 0.0;
        outer_binary_orientation[2] = 1.0;
    }

    // Initialize the inner binary in its own COM frame
    ic_binary(ode_params, inner_binary_params, m1, m2, w_inner);

    // Initialize the effective outer binary
    ic_binary(ode_params, outer_binary_params, m_inner, m3, w_outer);

    // Orient the outer binary around the total triple COM
    double triple_com_pos[3] = {0.0, 0.0, 0.0};
    position_binary(triple_com_pos, outer_binary_orientation, w_outer);

    double inner_com_pos[3] = {w_outer[0], w_outer[1], w_outer[2]};
    double inner_com_mom[3] = {w_outer[6], w_outer[7], w_outer[8]};
    double tertiary_pos[3] = {w_outer[3], w_outer[4], w_outer[5]};
    double tertiary_mom[3] = {w_outer[9], w_outer[10], w_outer[11]};

    // Orient the inner binary and place its COM on the outer orbit
    position_binary(inner_com_pos, inner_binary_orientation, w_inner);

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
 * @param[in]   v0_rel          Initial relative approach velocity of the binary and the single 
 * @param[in]   b               Scattering impact parameter
 * @param[in]   orientation     Binary orientation vector (NULL -> oriented in direction of motion)
 * @param[out]  w0              Initial positions and momenta, 
 *                              w0 = [pos_b1, pos_b2, pos_s, p_b1, p_b2, p_s]
 */
void ic_binary_single_scattering(struct ode_params* ode_params, 
    struct binary_params* binary_params, double d0, double v0_rel, double b, double* orientation, 
    double* w0)
{
    double v_x, v_y;
    double w0_binary[12];
    double binary_orientation[3];
    double binary_pos[3] = {0.0, 0.0, 0.0};

    // Create the binary initial values
    ic_binary(ode_params, binary_params, ode_params->masses[0], ode_params->masses[1], w0_binary);

    // Compute the absolute values of v_x and v_y
    v_x = v0_rel/2 * sqrt(1 - b*b/(d0*d0));
    v_y = v0_rel/2 * b/d0;

    // Position the binary at its right place with the specified orientation
    binary_pos[0] = -d0/2;
    if(orientation == NULL){
        binary_orientation[0] = v_x;
        binary_orientation[1] = -v_y;
        binary_orientation[2] = 0.0;
    }
    else{
        binary_orientation[0] = orientation[0];
        binary_orientation[1] = orientation[1];
        binary_orientation[2] = orientation[2];
    }
    position_binary(binary_pos, binary_orientation, w0_binary);

    // Add p_x and p_y momentum to the binary
    w0_binary[6] += ode_params->masses[0] * v_x;
    w0_binary[7] -= ode_params->masses[0] * v_y;
    w0_binary[9] += ode_params->masses[1] * v_x;
    w0_binary[10] -= ode_params->masses[1] * v_y;

    // Setting up the binary
    for(int i = 0; i < 6; i++)
        w0[i] = w0_binary[i];
    for(int i = 9; i < 15; i++)
        w0[i] = w0_binary[i-3];

    // Setting up the single body
    w0[6] = d0/2;
    w0[7] = 0.0;
    w0[8] = 0.0;
    w0[15] = -ode_params->masses[2] * v_x;
    w0[16] = ode_params->masses[2] * v_y;
    w0[17] = 0.0;
}


/**
 * @brief Computes initial positions and momenta for a relativistic binary-single scattering.
 * 
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, pos3, p1, p2, p3])
 * of a relativistic binary-single scattering with specified binary and scattering parameters. The
 * scattering takes place in the xy-plane. The function currently only supports equal-mass binaries
 * for which it is common to report the radial and tangential component of the binary members that
 * for example produce a quasi-circular orbit.
 * 
 * @param[in]   d0              Initial distance of the binary's center of mass and the single
 * @param[in]   p0_rel          Initial relative approach momentum of the binary and the single
 * @param[in]   b               Scattering impact parameter
 * @param[in]   binary_phi0     Binary phase shift
 * @param[in]   binary_r0       Initial binary radius
 * @param[in]   binary_pt0      Tangential component of the initial momentum of the binary members
 * @param[in]   binary_pr0      Radial component of the initial momentum of the binary members
 * @param[in]   orientation     Binary orientation vector (NULL -> oriented in direction of motion)
 * @param[out]  w0              Initial positions and momenta,
 *                              w0 = [pos_b1, pos_b2, pos_s, p_b1, p_b2, p_s]
 */
void ic_binary_single_scattering_rel(double d0, double p0_rel, double b, double binary_phi0, 
    double binary_r0, double binary_pt0, double binary_pr0, double* orientation, double* w0)
{
    double px0_rel, py0_rel;
    double binary_orientation[3];

    // Compute the x- and y- components of the momentum vector of the two systems
    px0_rel = p0_rel * sqrt(1 - b*b/(d0*d0));
    py0_rel = -p0_rel * b/d0;

    // If no orientation is given, the binary will be oriented along its momentum vector
    if(orientation == NULL){
        binary_orientation[0] = px0_rel;
        binary_orientation[1] = py0_rel;
        binary_orientation[2] = 0.0;
    }
    else{
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
    p1[0] += px0_rel;
    p1[1] += py0_rel;
    p2[0] += px0_rel;
    p2[1] += py0_rel;

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
    w0[15] = -px0_rel;
    w0[16] = -py0_rel;
    w0[17] = 0.0;
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
 * @param[in]   v0_rel          Initial relative approach velocity of the binaries
 * @param[in]   b               Scattering impact parameter
 * @param[in]   orientation_1   Orientation of binary 1 (NULL -> oriented in direction of motion)
 * @param[in]   orientation_2   Orientation of binary 2 (NULL -> oriented in direction of motion)
 * @param[out]  w0              Initial positions and momenta, 
 *                              w0 = [pos1_1, pos2_1, pos1_2, pos2_2, p1_1, p2_1, p1_2, p2_2]
 */
void ic_binary_binary_scattering(struct ode_params* ode_params, 
    struct binary_params* binary1_params, struct binary_params* binary2_params, double d0, 
    double v0_rel, double b, double* orientation_1, double* orientation_2, double* w0)
{
    double v_x, v_y;
    double w0_binary1[12];
    double w0_binary2[12];
    double binary1_orientation[3];
    double binary2_orientation[3];
    double binary1_pos[3] = {0.0, 0.0, 0.0};
    double binary2_pos[3] = {0.0, 0.0, 0.0};

    // Create the binary initial values
    ic_binary(ode_params, binary1_params, ode_params->masses[0], ode_params->masses[1], w0_binary1);
    ic_binary(ode_params, binary2_params, ode_params->masses[2], ode_params->masses[3], w0_binary2);

    // Compute the absolute values of v_x and v_y
    v_x = v0_rel/2 * sqrt(1 - pow(b/d0, 2));
    v_y = v0_rel/2 * b/d0;

    // Position the binaries at their right place with the specified orientations
    binary1_pos[0] = -d0/2;
    binary2_pos[0] = d0/2;
    if(orientation_1 == NULL){
        binary1_orientation[0] = v_x;
        binary1_orientation[1] = -v_y;
        binary1_orientation[2] = 0.0;
    }
    else{
        binary1_orientation[0] = orientation_1[0];
        binary1_orientation[1] = orientation_1[1];
        binary1_orientation[2] = orientation_1[2];
    }
    if(orientation_2 == NULL){
        binary2_orientation[0] = -v_x;
        binary2_orientation[1] = v_y;
        binary2_orientation[2] = 0.0;
    }
    else{
        binary2_orientation[0] = orientation_2[0];
        binary2_orientation[1] = orientation_2[1];
        binary2_orientation[2] = orientation_2[2];
    }
    position_binary(binary1_pos, binary1_orientation, w0_binary1);
    position_binary(binary2_pos, binary2_orientation, w0_binary2);

    // Add p_x and p_y momentum to the binaries
    w0_binary1[6] += ode_params->masses[0] * v_x;
    w0_binary1[7] -= ode_params->masses[0] * v_y;
    w0_binary1[9] += ode_params->masses[1] * v_x;
    w0_binary1[10] -= ode_params->masses[1] * v_y;
    w0_binary2[6] -= ode_params->masses[2] * v_x;
    w0_binary2[7] += ode_params->masses[2] * v_y;
    w0_binary2[9] -= ode_params->masses[3] * v_x;
    w0_binary2[10] += ode_params->masses[3] * v_y;

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
 * @param[in]   p0_rel          Initial relative approach momentum of the binary and the single
 * @param[in]   b               Scattering impact parameter
 * @param[in]   binary_phi0_X   Phase shift of binary X
 * @param[in]   binary_r0_X     Initial radius of binary X
 * @param[in]   binary_pt0_X    Tangential component of the initial momentum of binary X members
 * @param[in]   binary_pr0_X    Radial component of the initial momentum of binary X members
 * @param[in]   orientation_X   Orientation of binary X (NULL -> oriented in direction of motion)
 * @param[out]  w0              Initial positions and momenta,
 *                              w0 = [pos1_1, pos2_1, pos1_2, pos2_2, p1_1, p2_1, p1_2, p2_2]
 */
void ic_binary_binary_scattering_rel(double d0, double p0_rel, double b, double binary1_phi0, 
    double binary1_r0, double binary1_pt0, double binary1_pr0, double* orientation_1, 
    double binary_phi0_2, double binary2_r0, double binary2_pt0, double binary2_pr0, 
    double* orientation_2, double* w0)
{
    double scatter_px0, scatter_py0;
    double binary_orientation_1[3];
    double binary_orientation_2[3];

    // Compute the x- and y- components of the momentum vector of the two systems
    scatter_px0 = p0_rel * sqrt(1 - pow(b/d0, 2));
    scatter_py0 = -p0_rel * b/d0;

    if(orientation_1 == NULL){
        binary_orientation_1[0] = scatter_px0;
        binary_orientation_1[1] = scatter_py0;
        binary_orientation_1[2] = 0.0;
    }
    else{
        binary_orientation_1[0] = orientation_1[0];
        binary_orientation_1[1] = orientation_1[1];
        binary_orientation_1[2] = orientation_1[2];
    }
    if(orientation_2 == NULL){
        binary_orientation_2[0] = -scatter_px0;
        binary_orientation_2[1] = -scatter_py0;
        binary_orientation_2[2] = 0.0;
    }
    else{
        binary_orientation_2[0] = orientation_1[0];
        binary_orientation_2[1] = orientation_1[1];
        binary_orientation_2[2] = orientation_1[2];
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
    p1_1[0] += scatter_px0;
    p1_1[1] += scatter_py0;
    p2_1[0] += scatter_px0;
    p2_1[1] += scatter_py0;
    pos1_2[0] += d0/2;
    pos2_2[0] += d0/2;
    p1_2[0] -= scatter_px0;
    p1_2[1] -= scatter_py0;
    p2_2[0] -= scatter_px0;
    p2_2[1] -= scatter_py0;

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
