#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "initial_configurations.h"
#include "utils.h"
#include "pn_eom.h"


void newtonian_binary(double m1, double m2, double a, double e, double phi0, int dim, double* w0)
/*Returns the initial positions and the initial velocities (as w0 = [pos1, pos2, v1, v2]) of the two bodies with 
masses m1 and m2 in a binary system with relative semi-major axis a, eccentricity e and initial phase phi0
in the center of mass system and in the xy-plane. dim sets the number of dimensions (2 or 3).*/
{
    // Basic binary parameters
    double M = m1 + m2;                                             // Total mass
    double L = m1 * m2 * sqrt(2 * a / M * (1 - pow(e, 2)));         // Angular momentum
    double p = pow(L, 2) / (pow(m1 * m2, 2) / M);                   // Semi-latus rectum

    // Relative distance and its derivative at phi0
    double r0 = p / (1 + e * cos(phi0));
    double dr_dphi0 = e * p * sin(phi0) / pow(1 + e * cos(phi0), 2);

    // Positions and velocities of the individual masses in the center of mass frame
    double r0cosphi0 = r0 * cos(phi0);
    double r0sinphi0 = r0 * sin(phi0);
    double v_factor1 = dr_dphi0 * cos(phi0) - r0sinphi0;
    double v_factor2 = dr_dphi0 * sin(phi0) + r0cosphi0;
    double v1_factor = L / (m1 * pow(r0, 2));
    double v2_factor = -L / (m2 * pow(r0, 2));

    if(dim == 2){
        w0[0] = m2 / M * r0cosphi0;
        w0[1] = m2 / M * r0sinphi0;

        w0[2] = -m1 / M * r0cosphi0;
        w0[3] = -m1 / M * r0sinphi0;

        w0[4] = v1_factor * v_factor1;
        w0[5] = v1_factor * v_factor2;

        w0[6] = v2_factor * v_factor1;
        w0[7] = v2_factor * v_factor2;
    }
    else if(dim == 3){
        w0[0] = m2 / M * r0cosphi0;
        w0[1] = m2 / M * r0sinphi0;
        w0[2] = 0.0;

        w0[3] = -m1 / M * r0cosphi0;
        w0[4] = -m1 / M * r0sinphi0;
        w0[5] = 0.0;

        w0[6] = v1_factor * v_factor1;
        w0[7] = v1_factor * v_factor2;
        w0[8] = 0.0;

        w0[9] = v2_factor * v_factor1;
        w0[10] = v2_factor * v_factor2;
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
    double v1[3] = {w0[6], w0[7], w0[8]};
    double v2[3] = {w0[9], w0[10], w0[11]};
    double pos1_new[3], pos2_new[3], v1_new[3], v2_new[3];
    double axis[3];
    double angle;
    double rotation_matrix[3][3];

    normalize(orientation);

    // Compute axis of rotation (cross product of initial and target orientation)
    cross_product(initial_orientation, orientation, axis);

    // Compute angle of rotation (dot product) and create the rotation matrix
    angle = acos(dot_product(initial_orientation, orientation, 3));
    create_rotation_matrix(axis, angle, rotation_matrix);

    // Rotate positions and velocities of both masses
    rotate_vector(pos1, rotation_matrix, pos1_new);
    rotate_vector(pos2, rotation_matrix, pos2_new);
    rotate_vector(v1, rotation_matrix, v1_new);
    rotate_vector(v2, rotation_matrix, v2_new);

    // Update w0 with rotated positions and velocities
    w0[0] = pos1_new[0] + com_pos[0];
    w0[1] = pos1_new[1] + com_pos[1];
    w0[2] = pos1_new[2] + com_pos[2];
    w0[3] = pos2_new[0] + com_pos[0];
    w0[4] = pos2_new[1] + com_pos[1];
    w0[5] = pos2_new[2] + com_pos[2];
    w0[6] = v1_new[0];
    w0[7] = v1_new[1];
    w0[8] = v1_new[2];
    w0[9] = v2_new[0];
    w0[10] = v2_new[1];
    w0[11] = v2_new[2];
}


void binary_single_scattering_symmetric(double m_a, double m_b, double a, double e, double phi0,
                                        double d0, double v_rel, double b, double* orientation, double* w0)
/*Computes the inital values for a symmetric binary-single scattering setup. The arguments are:

m_a, m_b        masses of the binary members
a               semi-major axis of the binary
e               eccentricity of the binary
phi0            initial phase of the binary
d0              initial separation
v_rel           initial relative approach velocity   
b               scattering impact parameter           
orientation     binary orientation vector (if NULL the binary is oriented in its direction of motion)
w0              pointer to the array in which the initial positions and velocities of the bodies will be stored
                w0 = [pos_a, pos_b, pos_s, v_a, v_b, v_s]*/
{
    double v_x, v_y;
    double w0_binary[12];
    double binary_orientation[3];
    double binary_pos[3] = {0.0, 0.0, 0.0};

    // Create the binary initial values
    newtonian_binary(m_a, m_b, a, e, phi0, 3, w0_binary);

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

    // Add v_x and v_y velocity to the binary
    w0_binary[6] += v_x;
    w0_binary[7] -= v_y;
    w0_binary[9] += v_x;
    w0_binary[10] -= v_y;

    // Setting up the binary
    for(int i = 0; i < 6; i++)
        w0[i] = w0_binary[i];
    for(int i = 9; i < 15; i++)
        w0[i] = w0_binary[i-3];

    // Setting up the single body
    w0[6] = d0/2;
    w0[7] = 0.0;
    w0[8] = 0.0;
    w0[15] = -v_x;
    w0[16] = v_y;
    w0[17] = 0.0;

    return;
}


void binary_binary_scattering_symmetric(double m_a1, double m_b1, double m_a2, double m_b2, 
                                        double a1, double a2, double e1, double e2, double phi01, double phi02,
                                        double d0, double v_rel, double b, double* orientation1, double* orientation2, 
                                        double* w0)
/*Computes the inital values for a symmetric binary-binary scattering setup. The arguments are:

m_aX, m_bX      masses of the binary members and the scattering single body 
aX              semi-major axis of the binary
eX              eccentricity of the binary
phi0X           initial phase of the binary
d0              initial separation
v_rel           initial relative approach velocity   
b               scattering impact parameter           
orientationX    binary orientation vectors (if NULL the binary is oriented in its direction of motion)
w0              pointer to the array in which the initial positions and velocities of the bodies will be stored
                w0 = [pos_a1, pos_b1, pos_a2, pos_b2, v_a1, v_b1, v_a2, v_b2]*/
{
    double v_x, v_y;
    double w0_binary1[12];
    double w0_binary2[12];
    double binary1_orientation[3];
    double binary2_orientation[3];
    double binary1_pos[3] = {0.0, 0.0, 0.0};
    double binary2_pos[3] = {0.0, 0.0, 0.0};

    // Create the binary initial values
    newtonian_binary(m_a1, m_b1, a1, e1, phi01, 3, w0_binary1);
    newtonian_binary(m_a2, m_b2, a2, e2, phi02, 3, w0_binary2);

    // Compute the absolute values of v_x and v_y as well as the angle of the v-vector to the x-axis
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

    // Add v_x and v_y velocity to the binaries
    w0_binary1[6] += v_x;
    w0_binary1[7] -= v_y;
    w0_binary1[9] += v_x;
    w0_binary1[10] -= v_y;
    w0_binary2[6] -= v_x;
    w0_binary2[7] += v_y;
    w0_binary2[9] -= v_x;
    w0_binary2[10] += v_y;

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


void figure_eight_orbit(double width, struct ode_params* params, double* w0){
    double px, py, lambda;

    lambda = width/108.1;
    if (width > 10000 || width < 100)
        printf("Warning: width = %lf, figure eight orbit only accurate for 100 < width < 100000!\n", width);

    if(params->pn_terms[0] == 1 & params->pn_terms[1] == 0 & params->pn_terms[2] == 0 & params->pn_terms[3] == 0){
        px = -0.09324;
        py = -0.08647;
    }
    else if(params->pn_terms[0] == 1 & params->pn_terms[1] == 1 & params->pn_terms[2] == 0 & params->pn_terms[3] == 0){
        px = -sqrt(0.008693198284902323/lambda + 0.0007967218808497421/pow(lambda, 2) + 0.00013920898849859853/pow(lambda, 3) - 3.426433477622561e-06/pow(lambda, 4));
        py = -sqrt(0.007477506932144364/lambda + 0.0012750411566393282/pow(lambda, 2) + 0.00020194589670320117/pow(lambda, 3) - 5.316870339927638e-05/pow(lambda, 4));
    }
    else if(params->pn_terms[0] == 1 & params->pn_terms[1] == 1 & params->pn_terms[2] == 1 & params->pn_terms[3] == 0){
        px = -sqrt(0.008692910686038705/lambda + 0.0007977722653187864/pow(lambda, 2) + 6.351332174711012e-05/pow(lambda, 3) - 3.0103527312470005e-05/pow(lambda, 4));
        py = -sqrt(0.007477759360235814/lambda + 0.0012723359375445093/pow(lambda, 2) + 6.583447113705563e-05/pow(lambda, 3) - 5.3457362474338285e-06/pow(lambda, 4));
    }

    if (params->dim == 2){
        w0[0] = 97.0;
        w0[1] = -24.31;
        w0[2] = -97.0;
        w0[3] = 24.31;
        w0[4] = 0.0;
        w0[5] = 0.0;
        for (int i = 0; i < 6; i++)
            w0[i] *= lambda;
        w0[6] = -0.5*px;
        w0[7] = -0.5*py;
        w0[8] = -0.5*px;
        w0[9] = -0.5*py;
        w0[10] = px;
        w0[11] = py;
        if(params->pn_terms[0] == 1 & params->pn_terms[1] == 0 & params->pn_terms[2] == 0 & params->pn_terms[3] == 0){
            for (int i = 6; i < 12; i++)
                w0[i] *= 1/sqrt(lambda);
        }
    }
    else if (params->dim == 3){
        w0[0] = 97.0;
        w0[1] = -24.31;
        w0[2] = 0.0;
        w0[3] = -97.0;
        w0[4] = 24.31;
        w0[5] = 0.0;
        w0[6] = 0.0;
        w0[7] = 0.0;
        w0[8] = 0.0;
        for (int i = 0; i < 9; i++)
            w0[i] *= lambda;
        w0[9] = -0.5*px;
        w0[10] = -0.5*py;
        w0[11] = 0.0;
        w0[12] = -0.5*px;
        w0[13] = -0.5*py;
        w0[14] = 0.0;
        w0[15] = px;
        w0[16] = py;
        w0[17] = 0.0;
        if(params->pn_terms[0] == 1 & params->pn_terms[1] == 0 & params->pn_terms[2] == 0 & params->pn_terms[3] == 0){
            for (int i = 9; i < 18; i++)
                w0[i] *= 1/sqrt(lambda);
        }
    }
}