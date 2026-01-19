#ifndef PN_EOM_HAMILTONIANS_H
#define PN_EOM_HAMILTONIANS_H

static inline int role_of_body(int p, int body);
double ln_integral_sum(double* w, struct ode_params* params);
double UTT4_without_ln_integral(double* w, struct ode_params* params);
complex double UTT4_without_ln_integral_complex(complex double* w, struct ode_params* params);
void compute_dUTT4_dx(double*w, struct ode_params* params, double *dUdx);

double H0PN(double* w, struct ode_params* params);
complex double H0PN_complex(complex double* w, struct ode_params* params, int p_flag);
double H1PN(double* w, struct ode_params* params);
complex double H1PN_complex(complex double* w, struct ode_params* params, int p_flag);
double H2PN_threebody(double* w, struct ode_params* params, int p_flag);
double H2PN_nbody(double* w, struct ode_params* params, int p_flag, int utt4_flag);
complex double H2PN_nbody_base_complex(complex double* w, struct ode_params* params, int p_flag);

void update_eom_hamiltonian_fd(double *w, double *dwdt, double (*hamiltonian)(double*, struct ode_params*, int), double h, struct ode_params* params);
void update_eom_hamiltonian_cs(double *w, double *dwdt, complex double (*hamiltonian)(complex double*, struct ode_params*, int), double h, struct ode_params* params);

#endif