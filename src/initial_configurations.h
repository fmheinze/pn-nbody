#include "pn_eom.h"

#ifndef INITIAL_CONFIGURATIONS_H
#define INITIAL_CONFIGURATIONS_H

void newtonian_binary(struct ode_params* params, double* w0, double a, double e, double phi0);
void binary_single_scattering_symmetric(struct ode_params* params, double* w0, double a, double e, double phi0,
                                        double d0, double v_rel, double b, double* orientation);
void binary_single_scattering(struct ode_params* params, double* w0, double d0, double scatter_p0, double b, double phi0, 
                              double binary_r0, double binary_pt0, double binary_pr0, double* orientation, int invert_rotation);
void binary_binary_scattering_symmetric(struct ode_params* params, double* w0,
                                        double a1, double a2, double e1, double e2, double phi01, double phi02,
                                        double d0, double v_rel, double b, double* orientation1, double* orientation2);
void binary_binary_scattering(struct ode_params* params, double* w0, double d0, double scatter_p0, double b, 
                              double phi0_1, double phi0_2, double binary_r0_1, double binary_pt0_1, double binary_pr0_1,
                              double binary_r0_2, double binary_pt0_2, double binary_pr0_2, 
                              double* orientation_1, double* orientation_2, int invert_rotation_1, int invert_rotation_2);
void figure_eight_orbit(struct ode_params* params, double* w0, double width);

#endif