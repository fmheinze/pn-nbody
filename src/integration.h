#include "pn_eom.h"

#ifndef INTEGRATION_H
#define INTEGRATION_H

void cash_karp_update(double* y, double* y_new, double* y_err, int y_size, double x, double dx, 
                      void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params, double* k1);
void ode_step(double* y, int y_size, double* x, double *dx, double dx_max, double rel_error, 
              void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params);
void ode_integrator(double* w, void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params);

#endif