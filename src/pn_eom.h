#ifndef PN_EOM_H
#define PN_EOM_H

struct ode_params {
    int num_dim;
    int num_bodies;
    int use_impulse_method;
    int* pn_terms;
    double* masses;
};

void rhs_pn_twobody(double t, double* w, struct ode_params* params, double* dwdt);
void rhs_pn_nbody(double t, double* w, struct ode_params* params, double* dwdt);

#endif