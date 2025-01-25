#include "pn_eom.h"

#ifndef INITIAL_CONFIGURATIONS_H
#define INITIAL_CONFIGURATIONS_H

void newtonian_binary(double m1, double m2, double a, double e, double phi0, int dim, double* w0);
void binary_single_scattering_symmetric(double m_a, double m_b, double a, double e, double phi0,
                                        double d0, double v_rel, double b, double* orientation, double* w0);
void binary_binary_scattering_symmetric(double m_a1, double m_b1, double m_a2, double m_b2, 
                                        double a1, double a2, double e1, double e2, double phi01, double phi02,
                                        double d0, double v_rel, double b, double* orientation1, double* orientation2, 
                                        double* w0);
void figure_eight_orbit(double width, struct ode_params* params, double* w0);

#endif