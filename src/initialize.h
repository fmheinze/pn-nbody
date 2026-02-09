#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "eom.h"

void initialize_parameters(void);
struct ode_params initialize_ode_params(void);
struct binary_params initialize_binary_params(int i);
double* initialize_state_vector(struct ode_params* params);
void initialize_newtonian_binary(struct ode_params* ode_params, double* w0);
void initialize_binary_single_scattering(struct ode_params* ode_params, double* w0);
void initialize_binary_single_scattering_rel(struct ode_params* ode_params, double* w0);
void initialize_binary_binary_scattering(struct ode_params* ode_params, double* w0);
void initialize_binary_binary_scattering_rel(struct ode_params* ode_params, double* w0);
void initialize_figure_eight(struct ode_params* ode_params, double* w0);


#endif
