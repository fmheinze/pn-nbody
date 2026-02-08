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
#include "parameters.h"


/**
 * @brief Computes initial positions and momenta for a Newtonian binary.
 * 
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, p1, p2]) of two
 * bodies with masses m1 and m2 in a Newtonian binary system with relative semi-major axis a, 
 * eccentricity e, and initial phase phi0 in the center of mass frame and in the xy-plane.
 * 
 * @param[in]   params          Parameter struct containing general information about the system
 * @param[in]   binary_params   Struct containing the binary parameters (initialized elsewhere)
 * @param[out]  w0              Initial positions and momenta, w0 = [pos1, pos2, p1, p2]
 */
void ic_newtonian_binary(struct ode_params* ode_params, struct binary_params* binary_params,
    double* w0)
{
    // Unpack needed parameters
    double a = binary_params->a;
    double e = binary_params->e;
    double phi0 = binary_params->phi0;
    double p = binary_params->p;
    double m1 = ode_params->masses[0];
    double m2 = ode_params->masses[1];

    // Compute total mass, reduced mass and angular momentum
    double M = m1 + m2; 
    double mu = m1 * m2 / M;
    double L = mu * sqrt(a * M * (1 - e * e));

    // Relative distance and its derivative at phi0
    double cosphi0 = cos(phi0);
    double sinphi0 = sin(phi0);
    double denom_inv = 1 / (1 + e * cosphi0);
    double r0 = p * denom_inv;
    double dr_dphi0 = e * p * sinphi0 * denom_inv * denom_inv;

    // Positions and velocities of the individual masses in the center of mass frame
    double r0cosphi0 = r0 * cosphi0;
    double r0sinphi0 = r0 * sinphi0;
    double p_factor1 = dr_dphi0 * cosphi0 - r0sinphi0;
    double p_factor2 = dr_dphi0 * sinphi0 + r0cosphi0;
    double p1_factor = L / (r0 * r0);
    double p2_factor = -L / (r0 * r0);

    // Set initial parameters
    if (ode_params->num_dim == 2){
        w0[0] = m2 / M * r0cosphi0;
        w0[1] = m2 / M * r0sinphi0;

        w0[2] = -m1 / M * r0cosphi0;
        w0[3] = -m1 / M * r0sinphi0;

        w0[4] = p1_factor * p_factor1;
        w0[5] = p1_factor * p_factor2;

        w0[6] = p2_factor * p_factor1;
        w0[7] = p2_factor * p_factor2;
    }
    if (ode_params->num_dim == 3) {
        w0[0] = m2 / M * r0cosphi0;
        w0[1] = m2 / M * r0sinphi0;
        w0[2] = 0.0;

        w0[3] = -m1 / M * r0cosphi0;
        w0[4] = -m1 / M * r0sinphi0;
        w0[5] = 0.0;

        w0[6] = p1_factor * p_factor1;
        w0[7] = p1_factor * p_factor2;
        w0[8] = 0.0;

        w0[9] = p2_factor * p_factor1;
        w0[10] = p2_factor * p_factor2;
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
 * @brief Computes initial positions and momenta for a Newtonian binary-single scattering.
 * 
 * Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, pos3, p1, p2, p3])
 * of a Newtonian binary-single scattering with specified binary and scattering parameters. The
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
    ic_newtonian_binary(ode_params, binary_params, w0_binary);

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
 * @brief Computes initial positions and momenta for a Newtonian binary-binary scattering.
 * 
 * Returns the initial positions and momenta (as w0 = [pos1, pos2, pos3, pos4, p1, p2, p3, p4])
 * of a Newtonian binary-binary scattering with specified binary and scattering parameters. The
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
    ic_newtonian_binary(ode_params, binary1_params, w0_binary1);
    ic_newtonian_binary(ode_params, binary2_params, w0_binary2);

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
    double pos_x, pos_y, px, py, lambda;

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
    else if (params->pn_terms[0] == 1 & params->pn_terms[1] == 1 & params->pn_terms[2] == 0) {
        px = -sqrt(0.008693032833827606/lambda + 0.000798860400642637/pow(lambda, 2) 
                                               + 0.00013381114672890315/pow(lambda, 3));
        py = -sqrt(0.007480061222224325/lambda + 0.001241927410006741/pow(lambda, 2) 
                                               + 0.00028564641727617235/pow(lambda, 3));
    }
    // 2PN momenta
    else if (params->pn_terms[0] == 1 & params->pn_terms[1] == 1 & params->pn_terms[2] == 1) {
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
