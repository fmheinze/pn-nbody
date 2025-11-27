#ifndef PN_EOM_HAMILTONIANS_H
#define PN_EOM_HAMILTONIANS_H

double H0PN(double* w, struct ode_params* params);
double H1PN(double* w, struct ode_params* params);
double H2PN_threebody(double* w, struct ode_params* params, int p_flag);
double H2PN_nbody(double* w, struct ode_params* params, int p_flag);
void update_eom_hamiltonian(double *w, double *dwdt, double (*hamiltonian)(double*, struct ode_params*, int p_flag), double h, struct ode_params* params);

#endif