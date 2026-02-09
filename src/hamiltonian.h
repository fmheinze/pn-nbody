#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <complex.h>

struct ode_params;

typedef complex double (*c_hamiltonian)(complex double*, struct ode_params*, int);

#if HAVE_CUBA
void compute_dUTT4_dx(double*w, struct ode_params* ode_params, double *dUdx);
#endif

double H0PN(double* w, struct ode_params* ode_params);
double H1PN(double* w, struct ode_params* ode_params);
double H2PN(double* w, struct ode_params* ode_params, int utt4_flag);
complex double H2PN_base_complex(complex double* w, struct ode_params* ode_params, int p_flag);

#endif
