#include "pn_eom.h"

#ifndef INTEGRATION_H
#define INTEGRATION_H

void cash_karp_update(double* y, double* y_new, double* y_err, int y_size, double x, double dx, 
                      void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params, double* k1);
void cash_karp_driver(double* y, int y_size, double* x, double *dx, double dx_max, double rel_error, 
              void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params);
void ode_integrator(double* w, void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params);

typedef void (*utt4_grad_func)(double* w, struct ode_params* params, double* dUdx);
static void impulse_apply_kick(double* w, int num_bodies, int num_dim, double h, const double* dUdx);
static void impulse_advance_middle(double* w, int w_size, double t_start, double h, int n, const char* middle_method,
                                   double rel_error, void (*rhs_mid)(double, double*, struct ode_params*, double*),
                                   struct ode_params* params);
void impulse_integrator(double* w, void (*rhs_mid)(double, double*, struct ode_params*, double*),
                        utt4_grad_func grad_utt4, struct ode_params* params);

#endif