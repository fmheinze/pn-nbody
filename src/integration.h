#ifndef INTEGRATION_H
#define INTEGRATION_H

#define ODE_WS_MAXSIZE 9

struct ode_ws {
    int y_size;
    double *buf;
};

struct ode_params;

typedef enum { M_CASH_KARP, M_RK4, M_IMPLICIT_MIDPOINT, M_UNKNOWN } ode_method;

typedef void (*utt4_grad_func)(double* w, struct ode_params* ode_params, double* dUdx);
typedef void (*ode_rhs)(double t, double* w, struct ode_params* ode_params, double* dwdt);

void ode_integrator(double* w, ode_rhs rhs, struct ode_params* ode_params);
void ode_integrator_impulse(double* w, ode_rhs rhs_mid, utt4_grad_func grad_utt4, 
    struct ode_params* ode_params);


#endif
