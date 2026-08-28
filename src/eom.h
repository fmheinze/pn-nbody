#ifndef EOM_H
#define EOM_H

#include "hamiltonian.h"


struct PairCache;
struct UTT4Cache;
struct DynamicsCache;
struct StateKey;

struct ode_params {
    // General setup
    int num_dim;
    int num_bodies_initial;
    int* pn_terms;
    double* masses;

    // EOM parameters (frequent access, parameter database queries would be expensive)
    int use_impulse_method;
    int include_utt4;
    double utt4_epsrel;
    double utt4_epsabs;
    int utt4_min_order;
    int utt4_max_order;
    int utt4_adaptive;
    int utt4_max_depth;
    int utt4_parallel;
    int utt4_verify_interval;

    // Merger parameters
    int merge_activate;
    char* remnant_prescription;
    double merge_factor;

    // Merger history
    int num_active;
    int *active;
    struct PairCache *pair_cache;  // Persistent, owned RHS/energy workspace
    struct UTT4Cache *utt4_cache;  // Persistent, owned UTT4 evaluation cache
    struct DynamicsCache *dynamics_cache;  // Exact-state derived velocity/RHS cache
    struct StateKey *state_key;    // Shared key every cache above is validated against
    long long *body_id;
    long long next_body_id;
    int *generation;
};

void compute_dUTT4_dx(double* w, struct ode_params* ode_params, double *dUdx);
void compute_coordinate_velocities(double* w, struct ode_params* ode_params, double* velocities);
void rhs_pn_nbody(double t, double* w, struct ode_params* ode_params, double* dwdt);
void update_eom_hamiltonian_cs(double *w, c_hamiltonian H, double h, struct ode_params* ode_params,
    double *dwdt);

#endif
