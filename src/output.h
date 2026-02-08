#include <stdio.h>
#include "pn_eom.h"

#ifndef OUTPUT_H
#define OUTPUT_H

void output_init(FILE** file_masses, FILE** file_pos, FILE** file_mom, FILE** file_energy, struct ode_params* params);
void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_energy, struct ode_params* params, double* w, double t);

#endif
