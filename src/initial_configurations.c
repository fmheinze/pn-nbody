#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "initial_configurations.h"
#include "utils.h"
#include "pn_eom.h"


void newtonian_binary(struct ode_params* params, double* w0, double a, double e, double phi0)
/*Returns the initial positions and the initial momenta (as w0 = [pos1, pos2, p1, p2]) of the two bodies with 
masses m1 and m2 in a binary system with relative semi-major axis a, eccentricity e and initial phase phi0
in the center of mass system and in the xy-plane. dim sets the number of dimensions (2 or 3).*/
{
    // Basic binary parameters
    double m1 = params->masses[0];
    double m2 = params->masses[1];
    double M = m1 + m2; 
    double mu = m1 * m2 / M;                                    // Reduced mass
    double L = mu * sqrt(a * M * (1 - e * e));                  // Angular momentum
    double p = L * L / (mu * mu * M);                           // Semi-latus rectum

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

    if(params->dim == 2){
        w0[0] = m2 / M * r0cosphi0;
        w0[1] = m2 / M * r0sinphi0;

        w0[2] = -m1 / M * r0cosphi0;
        w0[3] = -m1 / M * r0sinphi0;

        w0[4] = p1_factor * p_factor1;
        w0[5] = p1_factor * p_factor2;

        w0[6] = p2_factor * p_factor1;
        w0[7] = p2_factor * p_factor2;
    }
    else if(params->dim == 3){
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
    else{
        fprintf(stderr, "Error: dim must be either 2 or 3!\n");
        exit(EXIT_FAILURE);
    }
    return;
}


void position_binary(double com_pos[3], double orientation[3], double w0[12])
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

    // Rotate positions and velocities of both masses
    rotate_vector(pos1, rotation_matrix, pos1_new);
    rotate_vector(pos2, rotation_matrix, pos2_new);
    rotate_vector(p1, rotation_matrix, p1_new);
    rotate_vector(p2, rotation_matrix, p2_new);

    // Update w0 with rotated positions and velocities
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


void binary_single_scattering_symmetric(struct ode_params* params, double* w0, double a, double e, double phi0,
                                        double d0, double v_rel, double b, double* orientation)
/*Computes the inital values for a symmetric Newtonian binary-single scattering setup. The arguments are:

params          contains masses = [m_b1, m_b2, m_s]
a               semi-major axis of the binary
e               eccentricity of the binary
phi0            initial phase of the binary
d0              initial separation
v_rel           initial relative approach velocity   
b               scattering impact parameter           
orientation     binary orientation vector (if NULL the binary is oriented in its direction of motion)
w0              pointer to the array in which the initial positions and velocities of the bodies will be stored
                w0 = [pos_b1, pos_b2, pos_s, p_b1, p_b2, p_s]*/
{
    double v_x, v_y;
    double w0_binary[12];
    double binary_orientation[3];
    double binary_pos[3] = {0.0, 0.0, 0.0};

    // Create the binary initial values
    newtonian_binary(params, w0_binary, a, e, phi0);

    // Compute the absolute values of v_x and v_y
    v_x = v_rel/2 * sqrt(1 - pow(b/d0, 2));
    v_y = v_rel/2 * b/d0;

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
    w0_binary[6] += params->masses[0] * v_x;
    w0_binary[7] -= params->masses[0] * v_y;
    w0_binary[9] += params->masses[1] * v_x;
    w0_binary[10] -= params->masses[1] * v_y;

    // Setting up the binary
    for(int i = 0; i < 6; i++)
        w0[i] = w0_binary[i];
    for(int i = 9; i < 15; i++)
        w0[i] = w0_binary[i-3];

    // Setting up the single body
    w0[6] = d0/2;
    w0[7] = 0.0;
    w0[8] = 0.0;
    w0[15] = -params->masses[2] * v_x;
    w0[16] = params->masses[2] * v_y;
    w0[17] = 0.0;

    return;
}


void binary_single_scattering(struct ode_params* params, double* w0, double d0, double scatter_p0, double b, double phi0, 
                              double binary_r0, double binary_pt0, double binary_pr0, double* orientation, int invert_rotation)
{
    double scatter_px0, scatter_py0;
    double binary_orientation[3];

    // Compute the x- and y- components of the momentum vector of the two systems
    scatter_px0 = scatter_p0 * sqrt(1 - pow(b/d0, 2));
    scatter_py0 = -scatter_p0 * b/d0;

    if(orientation == NULL){
        binary_orientation[0] = scatter_px0;
        binary_orientation[1] = scatter_py0;
        binary_orientation[2] = 0.0;
    }
    else{
        binary_orientation[0] = orientation[0];
        binary_orientation[1] = orientation[1];
        binary_orientation[2] = orientation[2];
    }

    // As a starting point set up the binary in the xy-plane at the origin
    double p1[3] = {binary_r0, 0.0, 0.0};
    double p2[3] = {-binary_r0, 0.0, 0.0};
    double m1[3] = {binary_pr0, -binary_pt0, 0.0};
    double m2[3] = {-binary_pr0, binary_pt0, 0.0};
    if (invert_rotation) {
        m1[0] *= -1;
        m1[1] *= -1;
        m2[0] *= -1;
        m2[1] *= -1;
    }
    
    // Perform a rotation to account for the specified phase shift
    double R[3][3];
    double axis[3] = {0.0, 0.0, 1.0};
    create_rotation_matrix(axis, phi0, R);
    rotate_vector(p1, R, p1);
    rotate_vector(p2, R, p2);
    rotate_vector(m1, R, m1);
    rotate_vector(m2, R, m2);

    // Perform a rotation to align the normal vector of the orbital plane with the orientation vector
    align_vectors_rotation_matrix(axis, binary_orientation, R);
    rotate_vector(p1, R, p1);
    rotate_vector(p2, R, p2);
    rotate_vector(m1, R, m1);
    rotate_vector(m2, R, m2);

    // Put the binary at the specified position and add the momentum vector for the scattering
    p1[0] -= d0/2;
    p2[0] -= d0/2;
    m1[0] += scatter_px0;
    m1[1] += scatter_py0;
    m2[0] += scatter_px0;
    m2[1] += scatter_py0;

    // Fill the state vector
    w0[0] = p1[0];
    w0[1] = p1[1];
    w0[2] = p1[2];
    w0[3] = p2[0];
    w0[4] = p2[1];
    w0[5] = p2[2];
    w0[6] = d0/2;
    w0[7] = 0.0;
    w0[8] = 0.0;
    w0[9] = m1[0];
    w0[10] = m1[1];
    w0[11] = m1[2];
    w0[12] = m2[0];
    w0[13] = m2[1];
    w0[14] = m2[2];
    w0[15] = -scatter_px0;
    w0[16] = -scatter_py0;
    w0[17] = 0.0;
}


void binary_binary_scattering(struct ode_params* params, double* w0, double d0, double scatter_p0, double b, 
                              double phi0_1, double phi0_2, double binary_r0_1, double binary_pt0_1, double binary_pr0_1,
                              double binary_r0_2, double binary_pt0_2, double binary_pr0_2, 
                              double* orientation_1, double* orientation_2, int invert_rotation_1, int invert_rotation_2)
{
    double scatter_px0, scatter_py0;
    double binary_orientation_1[3];
    double binary_orientation_2[3];

    // Compute the x- and y- components of the momentum vector of the two systems
    scatter_px0 = scatter_p0 * sqrt(1 - pow(b/d0, 2));
    scatter_py0 = -scatter_p0 * b/d0;

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
    double p1_1[3] = {binary_r0_1, 0.0, 0.0};
    double p2_1[3] = {-binary_r0_1, 0.0, 0.0};
    double m1_1[3] = {binary_pr0_1, -binary_pt0_1, 0.0};
    double m2_1[3] = {-binary_pr0_1, binary_pt0_1, 0.0};
    if (invert_rotation_1) {
        m1_1[0] *= -1;
        m1_1[1] *= -1;
        m2_1[0] *= -1;
        m2_1[1] *= -1;
    }
    double p1_2[3] = {binary_r0_2, 0.0, 0.0};
    double p2_2[3] = {-binary_r0_2, 0.0, 0.0};
    double m1_2[3] = {binary_pr0_2, -binary_pt0_2, 0.0};
    double m2_2[3] = {-binary_pr0_2, binary_pt0_2, 0.0};
    if (invert_rotation_1) {
        m1_2[0] *= -1;
        m1_2[1] *= -1;
        m2_2[0] *= -1;
        m2_2[1] *= -1;
    }
    
    // Perform a rotation to account for the specified phase shift
    double R[3][3];
    double axis[3] = {0.0, 0.0, 1.0};
    create_rotation_matrix(axis, phi0_1, R);
    rotate_vector(p1_1, R, p1_1);
    rotate_vector(p2_1, R, p2_1);
    rotate_vector(m1_1, R, m1_1);
    rotate_vector(m2_1, R, m2_1);
    create_rotation_matrix(axis, phi0_2, R);
    rotate_vector(p1_2, R, p1_2);
    rotate_vector(p2_2, R, p2_2);
    rotate_vector(m1_2, R, m1_2);
    rotate_vector(m2_2, R, m2_2);

    // Perform a rotation to align the normal vector of the orbital plane with the orientation vector
    align_vectors_rotation_matrix(axis, binary_orientation_1, R);
    rotate_vector(p1_1, R, p1_1);
    rotate_vector(p2_1, R, p2_1);
    rotate_vector(m1_1, R, m1_1);
    rotate_vector(m2_1, R, m2_1);
    align_vectors_rotation_matrix(axis, binary_orientation_2, R);
    rotate_vector(p1_2, R, p1_2);
    rotate_vector(p2_2, R, p2_2);
    rotate_vector(m1_2, R, m1_2);
    rotate_vector(m2_2, R, m2_2);

    // Put the binaries at the specified position and add the momentum vector for the scattering
    p1_1[0] -= d0/2;
    p2_1[0] -= d0/2;
    m1_1[0] += scatter_px0;
    m1_1[1] += scatter_py0;
    m2_1[0] += scatter_px0;
    m2_1[1] += scatter_py0;
    p1_2[0] += d0/2;
    p2_2[0] += d0/2;
    m1_2[0] -= scatter_px0;
    m1_2[1] -= scatter_py0;
    m2_2[0] -= scatter_px0;
    m2_2[1] -= scatter_py0;

    // Fill the state vector
    w0[0] = p1_1[0];
    w0[1] = p1_1[1];
    w0[2] = p1_1[2];
    w0[3] = p2_1[0];
    w0[4] = p2_1[1];
    w0[5] = p2_1[2];
    w0[6] = p1_2[0];
    w0[7] = p1_2[1];
    w0[8] = p1_2[2];
    w0[9] = p2_2[0];
    w0[10] = p2_2[1];
    w0[11] = p2_2[2];
    w0[12] = m1_1[0];
    w0[13] = m1_1[1];
    w0[14] = m1_1[2];
    w0[15] = m2_1[0];
    w0[16] = m2_1[1];
    w0[17] = m2_1[2];
    w0[18] = m1_2[0];
    w0[19] = m1_2[1];
    w0[20] = m1_2[2];
    w0[21] = m2_2[0];
    w0[22] = m2_2[1];
    w0[23] = m2_2[2];
}


void binary_binary_scattering_symmetric(struct ode_params* params, double* w0,
                                        double a1, double a2, double e1, double e2, double phi01, double phi02,
                                        double d0, double v_rel, double b, double* orientation1, double* orientation2)
/*Computes the inital values for a symmetric binary-binary scattering setup. The arguments are:

params          contains masses = [m_a1, m_a2, m_b1, m_b2]
aX              semi-major axis of the binary X
eX              eccentricity of the binary X
phi0X           initial phase of the binary X
d0              initial separation of the centers of masses of the binaries
v_rel           initial relative approach velocity
b               scattering impact parameter           
orientationX    binary X orientation vectors (if NULL the binary is oriented in its direction of motion)
w0              pointer to the array in which the initial positions and velocities of the bodies will be stored
                w0 = [pos_a1, pos_a2, pos_b1, pos_b2, p_a1, p_a2, p_b1, p_b2]*/
{
    double v_x, v_y;
    double w0_binary1[12];
    double w0_binary2[12];
    double binary1_orientation[3];
    double binary2_orientation[3];
    double binary1_pos[3] = {0.0, 0.0, 0.0};
    double binary2_pos[3] = {0.0, 0.0, 0.0};

    // Create the binary initial values
    newtonian_binary(params, w0_binary1, a1, e1, phi01);
    newtonian_binary(params, w0_binary2, a2, e2, phi02);

    // Compute the absolute values of v_x and v_y
    v_x = v_rel/2 * sqrt(1 - pow(b/d0, 2));
    v_y = v_rel/2 * b/d0;

    // Position the binaries at their right place with the specified orientations
    binary1_pos[0] = -d0/2;
    binary2_pos[0] = d0/2;
    if(orientation1 == NULL){
        binary1_orientation[0] = v_x;
        binary1_orientation[1] = -v_y;
        binary1_orientation[2] = 0.0;
    }
    else{
        binary1_orientation[0] = orientation1[0];
        binary1_orientation[1] = orientation1[1];
        binary1_orientation[2] = orientation1[2];
    }
    if(orientation2 == NULL){
        binary2_orientation[0] = -v_x;
        binary2_orientation[1] = v_y;
        binary2_orientation[2] = 0.0;
    }
    else{
        binary2_orientation[0] = orientation2[0];
        binary2_orientation[1] = orientation2[1];
        binary2_orientation[2] = orientation2[2];
    }
    position_binary(binary1_pos, binary1_orientation, w0_binary1);
    position_binary(binary2_pos, binary2_orientation, w0_binary2);

    // Add p_x and p_y momentum to the binaries
    w0_binary1[6] += params->masses[0] * v_x;
    w0_binary1[7] -= params->masses[0] * v_y;
    w0_binary1[9] += params->masses[1] * v_x;
    w0_binary1[10] -= params->masses[1] * v_y;
    w0_binary2[6] -= params->masses[2] * v_x;
    w0_binary2[7] += params->masses[2] * v_y;
    w0_binary2[9] -= params->masses[3] * v_x;
    w0_binary2[10] += params->masses[3] * v_y;

    // Setting up the binaries
    for(int i = 0; i < 6; i++)
        w0[i] = w0_binary1[i];
    for(int i = 6; i < 12; i++)
        w0[i] = w0_binary2[i-6];
    for(int i = 12; i < 18; i++)
        w0[i] = w0_binary1[i-6];
    for(int i = 18; i < 24; i++)
        w0[i] = w0_binary2[i-12];

    return;
}


void figure_eight_orbit(struct ode_params* params, double* w0, double width){
    double pos_x, pos_y, px, py, lambda;

    lambda = width/108.1;
    if (width > 10000 || width < 100)
        printf("Warning: width = %lf, figure-eight orbit only accurate for 100 < width < 100000!\n", width);
    
    pos_x = 97.0 * lambda;
    pos_y = -24.31 * lambda;

    if(params->pn_terms[0] == 1 & params->pn_terms[1] == 0 & params->pn_terms[2] == 0 & params->pn_terms[3] == 0){
        px = -0.09324/sqrt(lambda);
        py = -0.08647/sqrt(lambda);
    }
    else if(params->pn_terms[0] == 1 & params->pn_terms[1] == 1 & params->pn_terms[2] == 0 & params->pn_terms[3] == 0){
        px = -sqrt(0.008693032833827606/lambda + 0.000798860400642637/pow(lambda, 2) + 0.00013381114672890315/pow(lambda, 3));
        py = -sqrt(0.007480061222224325/lambda + 0.001241927410006741/pow(lambda, 2) + 0.00028564641727617235/pow(lambda, 3));
    }
    else if(params->pn_terms[0] == 1 & params->pn_terms[1] == 1 & params->pn_terms[2] == 1 & params->pn_terms[3] == 0){
        px = -sqrt(0.008692910686038705/lambda + 0.0007977722653187864/pow(lambda, 2) + 6.351332174711012e-05/pow(lambda, 3) - 3.0103527312470005e-05/pow(lambda, 4));
        py = -sqrt(0.007477759360235814/lambda + 0.0012723359375445093/pow(lambda, 2) + 6.583447113705563e-05/pow(lambda, 3) - 5.3457362474338285e-06/pow(lambda, 4));
    }
    else if(params->pn_terms[0] == 1 & params->pn_terms[1] == 1 & params->pn_terms[2] == 1 & params->pn_terms[3] == 1){
        px = -sqrt(0.008692910686038705/lambda + 0.0007977722653187864/pow(lambda, 2) + 6.351332174711012e-05/pow(lambda, 3) - 3.0103527312470005e-05/pow(lambda, 4));
        py = -sqrt(0.007477759360235814/lambda + 0.0012723359375445093/pow(lambda, 2) + 6.583447113705563e-05/pow(lambda, 3) - 5.3457362474338285e-06/pow(lambda, 4));
    }

    if (params->dim == 2){
        w0[0] = pos_x;
        w0[1] = pos_y;
        w0[2] = -pos_x;
        w0[3] = -pos_y;
        w0[4] = 0.0;
        w0[5] = 0.0;
        w0[6] = -0.5*px;
        w0[7] = -0.5*py;
        w0[8] = -0.5*px;
        w0[9] = -0.5*py;
        w0[10] = px;
        w0[11] = py;
    }
    else if (params->dim == 3){
        w0[0] = pos_x;
        w0[1] = pos_y;
        w0[2] = 0.0;
        w0[3] = -pos_x;
        w0[4] = -pos_y;
        w0[5] = 0.0;
        w0[6] = 0.0;
        w0[7] = 0.0;
        w0[8] = 0.0;
        w0[9] = -0.5*px;
        w0[10] = -0.5*py;
        w0[11] = 0.0;
        w0[12] = -0.5*px;
        w0[13] = -0.5*py;
        w0[14] = 0.0;
        w0[15] = px;
        w0[16] = py;
        w0[17] = 0.0;
    }
}