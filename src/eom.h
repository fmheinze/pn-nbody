#ifndef EOM_H
#define EOM_H

#include "hamiltonian.h"

struct ode_params {
    int num_dim;
    int num_bodies;
    int* pn_terms;
    double* masses;

    int use_impulse_method;
    int include_utt4;
    int utt4_mineval;
    int utt4_maxeval;
    double utt4_epsrel;
    double utt4_epsabs;
};

void rhs_pn_nbody(double t, double* w, struct ode_params* ode_params, double* dwdt);
void update_eom_hamiltonian_cs(double *w, c_hamiltonian H, double h, struct ode_params* ode_params,
    double *dwdt);

#endif
