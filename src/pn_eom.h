#ifndef PN_EOM_H
#define PN_EOM_H

struct ode_params{
    int dim;
    int num_bodies;
    int* pn_terms;
    double* masses;
};

void rhs_pn_twobody(double t, double* w, struct ode_params* params, double* dwdt);
void rhs_pn_threebody(double t, double* w, struct ode_params* params, double* dwdt);

#endif