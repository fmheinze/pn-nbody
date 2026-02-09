#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>

struct ode_params;

void output_init(FILE** file_masses, FILE** file_pos, FILE** file_mom, FILE** file_energy,
    struct ode_params* ode_params);
void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_energy,
    struct ode_params* ode_params, double* w, double t);

#endif
