#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "integration.h"
#include "pn_eom.h"
#include "output.h"
#include "parameters.h"


void rk4_update(double* y, double* y_new, int y_size, double x, double dx,
                void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params) {
    //printf("%lf\n", dx);
    double dx_half = 0.5 * dx;
    double x_dx_half = x + dx_half;

    double *y_temp, *k1, *k_temp1, *k_temp2;
    allocate_vector(&y_temp, y_size);
    allocate_vector(&k1, y_size);
    allocate_vector(&k_temp1, y_size);
    allocate_vector(&k_temp2, y_size);

    // Compute k1/dx
    ode_rhs(x, y, params, k1);

    // Compute k2/dx
    for (int i = 0; i < y_size; i++)
        y_temp[i] = y[i] + dx_half * k1[i];
    (*ode_rhs)(x_dx_half, y_temp, params, k_temp1);

    // Compute k3/dx
    for (int i = 0; i < y_size; i++)
        y_temp[i] = y[i] + dx_half * k_temp1[i];
    (*ode_rhs)(x_dx_half, y_temp, params, k_temp2);

    // Compute k4/dx
    for (int i = 0; i < y_size; i++) {
        y_temp[i] = y[i] + dx * k_temp2[i];
        k_temp1[i] += k_temp2[i]; // Combine k2/dx and k3/dx in the same variable as they have the same coefficient later
    }
    (*ode_rhs)(x + dx, y_temp, params, k_temp2);

    // Update y
    for (int i = 0; i < y_size; i++)
        y_new[i] = y[i] + dx/6 * (k1[i] + 2 * k_temp1[i] + k_temp2[i]);

    free_vector(y_temp);
    free_vector(k1);
    free_vector(k_temp1);
    free_vector(k_temp2);
    return;
}


void cash_karp_update(double* y, double* y_new, double* y_err, int y_size, double x, double dx, 
                      void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params, double* k1)
/* Updates the value or the array of values y by one step dx according to a first-order system of ordinary differential 
equations and updates the error estimates y_err (differences of the 5th and 4th order Runge-Kutta step).

y           pointer to the value or array of values y
y_new       pointer to the value or array of values that is going to be filled with the updated values
y_err       pointer to the value or array of values that is going to be filled with the error estimates
y_size      number of values in the array y
x           value of the independent variable of the function y(x)
dx          stepsize
ode_rhs     pointer to the function which takes x, y and outputs y'(x) (in the last argument)
params      additional parameters for ode_rhs that modify the dynamics of the system (can be set to NULL)
k1          pointer to the value or the array of values of k1/dx = ode_rhs(x, y, k1) = y'(x) for better performance 
            (for multiple iterations in adaptive stepsize control k1 doesn't need to be computed multiple times) */
{
    // Coefficients for the Cash-Karp method (5th order Runge-Kutta result, the 4th order result is not computed)
    static const double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
                        b21 = 0.2,
                        b31 = 0.075, b32 = 0.225, 
                        b41 = 0.3, b42 = -0.9, b43 = 1.2,
                        b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0, b54 = 35.0/27.0,
                        b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
                        c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    
    // Differences in the c-coefficients for computing the error estimate
                        dc1 = -277.0/64512.0, dc3 = 6925.0/370944.0, dc4 = -6925.0/202752.0, dc5 = -277.0/14336.0, dc6 = 277.0/7084.0;
    
    double *k2, *k3, *k4, *k5, *k6, *y_temp;
    allocate_vector(&k2, y_size);
    allocate_vector(&k3, y_size);
    allocate_vector(&k4, y_size);
    allocate_vector(&k5, y_size);
    allocate_vector(&k6, y_size);
    allocate_vector(&y_temp, y_size);

    // Compute k2/dx
    for (int i = 0; i < y_size; i++)
        y_temp[i] = y[i] + b21 * dx * k1[i];
    (*ode_rhs)(x + a2 * dx, y_temp, params, k2);

    // Compute k3/dx
    for (int i = 0; i < y_size; i++)
        y_temp[i] = y[i] + dx * (b31 * k1[i] + b32 * k2[i]);
    (*ode_rhs)(x + a3 * dx, y_temp, params, k3);

    // Compute k4/dx
    for (int i = 0; i < y_size; i++)
        y_temp[i] = y[i] + dx * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
    (*ode_rhs)(x + a4 * dx, y_temp, params, k4);

    // Compute k5/dx
    for (int i = 0; i < y_size; i++)
        y_temp[i] = y[i] + dx * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
    (*ode_rhs)(x + a5 * dx, y_temp, params, k5);

    // Compute k6/dx
    for (int i = 0; i < y_size; i++)
        y_temp[i] = y[i] + dx * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
    (*ode_rhs)(x + a6 * dx, y_temp, params, k6);

    // Update y (c2 and c5 are zero for RK5)
    for (int i = 0; i < y_size; i++)
        y_new[i] = y[i] + dx * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i]);
    
    // Compute the error estimates (c2 is also zero for RK4)
    for (int i = 0; i < y_size; i++)
        y_err[i] = dx * (dc1 * k1[i] + dc3 * k3[i] + dc4 * k4[i] + dc5 * k5[i] + dc6 * k6[i]);

    free_vector(k2);
    free_vector(k3);
    free_vector(k4);
    free_vector(k5);
    free_vector(k6);
    free_vector(y_temp);
    return;
}


void cash_karp_driver(double* y, int y_size, double* x, double *dx, double dx_max, double rel_error, 
              void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params)
/* Updates the value or the array of values y by one adaptive step using the Cash-Karp update method,
ensuring a specified error threshold per step. The next x-value and the stepsize also get updated.

y           pointer to the value or array of values y that is going to be updated (updated values overwrite old ones)
y_size      number of values in the array y
x           value of the independent variable of the function y(x)
dx          stepsize (might be modified and updated by the adaptive stepsize procedure to ensure specified accuracy)
dx_max      maximum step size
rel_error   relative error threshold per step size
ode_rhs     pointer to the function which takes x, y and outputs y'(x) (in the last argument) 
params      additional parameters for ode_rhs that modify the dynamics of the system (can be set to NULL) */
{
    double *k1, *y_err, *y_new;
    allocate_vector(&k1, y_size);
    allocate_vector(&y_err, y_size);
    allocate_vector(&y_new, y_size);

    // Compute k1/dx (the same for each time step, only needs to be computed once)
    ode_rhs(*x, y, params, k1);

    // Compute steps as long as the maximum relative error is larger than the specified relative error threshold
    while (1) {
        // Perform a step with the Cash-Karp method and get the updated values y and the error estimates y_err
        cash_karp_update(y, y_new, y_err, y_size, *x, *dx, ode_rhs, params, k1);

        // Compute max relative error
        double max_rel_error = 0.0;
        for (int i = 0; i < y_size; i++) {
            double denom = fabs(y[i]) + *dx * fabs(k1[i]);
            if (denom != 0)
                max_rel_error = fmax(max_rel_error, fabs(y_err[i]) / denom);
        }

        // Compute new step size
        double dx_new = *dx * pow(max_rel_error / rel_error, -0.2);
        dx_new = fmin(dx_new, dx_max);  // Ensure step size does not exceed max limit

        if (max_rel_error < rel_error) {
            *x += *dx;  // Accept step
            *dx = dx_new;
            break;
        }

        // Otherwise, reject step and retry with smaller step size
        *dx = dx_new;
        if (*dx == dx_max) break;
    }

    // Update final new y values 
    for(int i = 0; i < y_size; i++)
        y[i] = y_new[i];

    free_vector(k1);
    free_vector(y_err);
    free_vector(y_new);
    return;
}


void ode_integrator(double* w, void (*ode_rhs)(double, double*, struct ode_params*, double*), struct ode_params* params)
/* Integrates the ODE system w'(t) = ode_rhs(t, w) from t_start to t_end, using the ode_step function which employs
an adaptive stepsize control ensuring a specified maximum error per step. The result of each step is written to a file.

w           pointer to the value or array of values w that is going to be updated (updated values overwrite old ones)
ode_rhs     pointer to the function which takes t, w and outputs w'(x) (in the last argument)
ode_params  additional parameters for ode_rhs that modify the dynamics of the system */
{
    double t_end = get_parameter_double("t_end");
    double dt = get_parameter_double("dt");
    double dt_save = get_parameter_double("dt_save");
    double rel_error = get_parameter_double("rel_error");
    char *method = get_parameter_string("ode_integrator");

    double w_size = 2 * params->num_dim * params->num_bodies;
    double t_current = 0.0;
    double t_save = 0.0;

    FILE *file_masses, *file_pos, *file_mom, *file_energy;

    output_init(&file_masses, &file_pos, &file_mom, &file_energy, params);
    output_write_timestep(file_pos, file_mom, file_energy, params, w, t_current);

    // Carry out integration steps until x_end is reached and write values into new line
    while (t_current < t_end){
        if(strcmp(method, "cash-karp") == 0)
            cash_karp_driver(w, w_size, &t_current, &dt, dt_save, rel_error, ode_rhs, params);
        if(strcmp(method, "rk4") == 0) {
            rk4_update(w, w, w_size, t_current, dt, ode_rhs, params);
            t_current += dt;
        }
        if(t_current >= t_save + dt_save){
            output_write_timestep(file_pos, file_mom, file_energy, params, w, t_current);
            print_progress_bar(t_current/t_end*100);
            t_save += dt_save;
        }
    }
    fclose(file_masses);
    fclose(file_pos);
    fclose(file_mom);
    fclose(file_energy);
    return;
}