#ifndef ORBTIAL_ELEMENTS_H
#define ORBTIAL_ELEMENTS_H

void total_angular_momentum_com(int* bodies, int num_bodies, double* w, struct ode_params* params, double* L);
double eccentricity(int* bodies, double* w, struct ode_params* params);
void semi_major_axes(int* bodies, double* w, struct ode_params* params, double* a, double* a1_com, double* a2_com);

#endif