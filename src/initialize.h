#include "pn_eom.h"

#ifndef OUTPUT_H
#define OUTPUT_H

void initialize_parameters();
double* initialize_state_vector(struct ode_params* params);
struct ode_params initialize_ode_params();


#endif