#ifndef PN_EOM_HAMILTONIANS_H
#define PN_EOM_HAMILTONIANS_H

double H0PN(double* w, struct ode_params* params);
double H1PN(double* w, struct ode_params* params);
double H2PN(double* w, struct ode_params* params);
double H2PN_threebody(double* w, struct ode_params* params);
void update_eom_hamiltonian(double *w, double *dwdt, double (*hamiltonian)(double*, struct ode_params*), double h, struct ode_params* params);

#endif