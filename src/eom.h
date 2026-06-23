#ifndef EOM_H
#define EOM_H

#include "hamiltonian.h"

struct ode_params {
    // General setup
    int num_dim;
    int num_bodies_initial;
    int* pn_terms;
    double* masses;

    // EOM parameters (frequent access, parameter database queries would be expensive)
    int use_impulse_method;
    int include_utt4;
    int utt4_mineval;
    int utt4_maxeval;
    double utt4_epsrel;
    double utt4_epsabs;

    // Merger parameters
    int merge_activate;
    char* remnant_prescription;
    double merge_factor;

    // Merger history
    int num_active;
    int *active;
    long long *body_id;
    long long next_body_id;
    int *generation;
};

void rhs_pn_nbody(double t, double* w, struct ode_params* ode_params, double* dwdt);
void update_eom_hamiltonian_cs(double *w, c_hamiltonian H, double h, struct ode_params* ode_params,
    double *dwdt);

#endif
