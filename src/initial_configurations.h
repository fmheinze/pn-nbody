#include "pn_eom.h"

#ifndef INITIAL_CONFIGURATIONS_H
#define INITIAL_CONFIGURATIONS_H


// Binary struct
struct binary_params {
    double a;       // Semi-major axis
    double b;       // Semi-minor axis
    double e;       // Eccentricity
    double r_a;     // Apoapsis distance
    double r_p;     // Periapsis distance
    double p;       // Semi-latus rectum
    double phi0;    // Initial phase
};


void ic_newtonian_binary(struct ode_params* ode_params, struct binary_params* binary_params,
    double* w0);

void ic_binary_single_scattering(struct ode_params* ode_params,
    struct binary_params* binary_params, double d0, double v0_rel, double b, double* orientation,
    double* w0);

void ic_binary_single_scattering_rel(double d0, double p0_rel, double b, double binary_phi0,
    double binary_r0, double binary_pt0, double binary_pr0, double* orientation, double* w0);

void ic_binary_binary_scattering(struct ode_params* ode_params, 
    struct binary_params* binary1_params, struct binary_params* binary2_params, double d0, 
    double v0_rel, double b, double* orientation_1, double* orientation_2, double* w0);

void ic_binary_binary_scattering_rel(double d0, double p0_rel, double b, double binary1_phi0, 
    double binary1_r0, double binary1_pt0, double binary1_pr0, double* orientation_1, 
    double binary_phi0_2, double binary2_r0, double binary2_pt0, double binary2_pr0, 
    double* orientation_2, double* w0);
    
void ic_figure_eight_orbit(struct ode_params* params, double width, double* w0);

#endif
