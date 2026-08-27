#ifndef EOM_H
#define EOM_H

#include "hamiltonian.h"


struct PairCache;

struct ode_params {
    // General setup
    int num_dim;
    int num_bodies_initial;
    int* pn_terms;
    double* masses;

    // EOM parameters (frequent access, parameter database queries would be expensive)
    int use_impulse_method;
    int include_utt4;
    double utt4_ln_integral_epsrel;
    double utt4_ln_integral_epsabs;
    int utt4_ln_integral_min_order;
    int utt4_ln_integral_max_order;
    int utt4_ln_integral_adaptive;
    int utt4_ln_integral_max_depth;
    int utt4_ln_integral_parallel;

    // Merger parameters
    int merge_activate;
    char* remnant_prescription;
    double merge_factor;

    // Merger history
    int num_active;
    int *active;
    struct PairCache *pair_cache;  // Persistent, owned RHS/energy workspace
    long long *body_id;
    long long next_body_id;
    int *generation;
};

void compute_dUTT4_dx(double* w, struct ode_params* ode_params, double *dUdx);
void rhs_pn_nbody(double t, double* w, struct ode_params* ode_params, double* dwdt);
void update_eom_hamiltonian_cs(double *w, c_hamiltonian H, double h, struct ode_params* ode_params,
    double *dwdt);

#endif
