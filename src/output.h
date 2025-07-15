#include <stdio.h>
#include "pn_eom.h"

#ifndef OUTPUT_H
#define OUTPUT_H

void general_output_init(FILE** file, struct ode_params* params);
void general_output_write(FILE* file, struct ode_params* params, double* w, double t);
void eccentricity_output_init(FILE** file);
void eccentricity_output_write(FILE* file, struct ode_params* params, double* w, double t);
void semi_major_axes_output_init(FILE** file);
void semi_major_axes_output_write(FILE* file, struct ode_params* params, double* w, double t);

#endif