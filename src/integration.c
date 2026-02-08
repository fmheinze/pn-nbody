/**
 * @file integration.c
 * @brief Routines for the numerical integration of ODEs
 *
 * Routines for the numerical integration of ODEs, including single-step methods, method drivers,
 * high-level control flow for method selection and writing output, as well as helper functions.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "integration.h"
#include "pn_eom.h"
#include "output.h"
#include "parameters.h"


// ODE workspace functions to avoid allocating and freeing memory at each ODE integration step

static void ode_ws_init(struct ode_ws* ws, int y_size) {
    if (ws->buf) free_vector(ws->buf);
    ws->y_size = y_size;
    ws->buf = NULL;
    allocate_vector(&ws->buf, ODE_WS_MAXSIZE * y_size);
}

static void ode_ws_free(struct ode_ws* ws) {
    free_vector(ws->buf);
    ws->buf = NULL;
    ws->y_size = 0;
}

static inline double* ws_vec(struct ode_ws* ws, int idx) {
    return ws->buf + (size_t)idx * (size_t)ws->y_size;
}

static inline void ws_check(const struct ode_ws* ws, int y_size) {
    if (!ws || !ws->buf || ws->y_size != y_size)
        errorexit("ode_ws not initialized or size mismatch");
}


// Checks the validity of the user-specified parameters for the numerical ODE integration
static void check_integration_parameter_validity()
{
    double t_end = get_parameter_double("t_end");
    double dt = get_parameter_double("dt");
    int impulse_method = get_parameter_int("impulse_method");
    if(t_end < 0.0) errorexit("Please specify a valid t_end (must be t_end >= 0)");
    if(dt <= 0.0) errorexit("Please specify a valid dt (must be dt > 0)");
    if(impulse_method != 0 && impulse_method != 1) 
        errorexit("Please set impulse_method to 0 (off) or 1 (on)");

    // Cash-Karp method
    if(strcmp(get_parameter_string("ode_integrator"), "cash-karp") == 0) {
        double rtol = get_parameter_double("rtol");
        if(rtol <= 0.0) 
            errorexit("Please specify a valid rtol of the Cash-Karp method (must be rtol > 0)");
    }

    // Implicit-midpoint method
    if(strcmp(get_parameter_string("ode_integrator"), "implicit-midpoint") == 0) {
        double tol = get_parameter_double("tol");
        int max_iter = get_parameter_int("max_iter");
        if(tol <= 0.0)
            errorexit("Please specify a valid tol (must be tol > 0)");
        if(max_iter <= 0) errorexit("Please specify a valid max_iter (must be max_iter > 0)");
    }

    // Impulse method
    if(get_parameter_int("impulse_method") == 1) {
        int impulse_method_n = get_parameter_int("impulse_method_n");
        if(impulse_method_n <= 0) 
            errorexit("Please specify a valid impulse_method_n (must be impulse_method_n > 0)");
    }
}


/**
 * @brief Performs a fourth-order Runge-Kutta timestep.
 *
 * Performs a fourth-order Runge-Kutta timestep, updating the state w according to the ODE
 * w'(t) = ode_rhs(t, w).
 * 
 * @param[in,out]   w           Input state (gets updated)
 * @param[in]       w_size      Number of values in state w
 * @param[in]       t           Current time
 * @param[in]       dt          Timestep
 * @param[in]       ode_rhs     Pointer to the function that specifies the ODE right-hand side
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[in]       ws          ODE workspace for storing intermediate results
 */
static void rk4_update(double* w, int w_size, double t, double dt, ode_rhs ode_rhs, 
    struct ode_params* ode_params, struct ode_ws* ws) 
{
    // ODE workspace setup
    ws_check(ws, w_size);
    double *k1      = ws_vec(ws, 0);
    double *w_temp  = ws_vec(ws, 1);
    double *k_temp1 = ws_vec(ws, 2);
    double *k_temp2 = ws_vec(ws, 3);

    double dt_half = 0.5 * dt;
    double t_dt_half = t + dt_half;

    // Compute k1/dt
    ode_rhs(t, w, ode_params, k1);

    // Compute k2/dt
    for (int i = 0; i < w_size; i++)
        w_temp[i] = w[i] + dt_half * k1[i];
    ode_rhs(t_dt_half, w_temp, ode_params, k_temp1);

    // Compute k3/dt
    for (int i = 0; i < w_size; i++)
        w_temp[i] = w[i] + dt_half * k_temp1[i];
    ode_rhs(t_dt_half, w_temp, ode_params, k_temp2);

    // Compute k4/dt
    for (int i = 0; i < w_size; i++) {
        w_temp[i] = w[i] + dt * k_temp2[i];
        // Combine k2/dt and k3/dt in the same variable as they have the same coefficient later
        k_temp1[i] += k_temp2[i];
    }
    ode_rhs(t + dt, w_temp, ode_params, k_temp2);

    // Update y
    for (int i = 0; i < w_size; i++)
        w[i] = w[i] + dt/6 * (k1[i] + 2 * k_temp1[i] + k_temp2[i]);
}


/**
 * @brief Performs a timestep with the implicit midpoint method.
 *
 * Performs a timestep with the implicit midpoint method, updating the state w according to the
 * ODE w'(t) = ode_rhs(t, w). The implicit midpoint method is second order and both symmetric
 * and symplectic. The implicit equation is solved using fixed-point iterations with a convergence
 * criterion that is based on the magnitude of the iteration update. The result of the fixed-point
 * iteration is returned (0 converged, 1 failed).
 * 
 * @param[in,out]   w           State (gets updated)
 * @param[in]       w_size      Number of values in state w
 * @param[in]       t           Current time
 * @param[in]       dt          Timestep
 * @param[in]       ode_rhs     Pointer to the function that specifies the ODE right-hand side
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[in]       tol         Tolerance for the convergence criterion
 * @param[in]       max_iter    Maximum number of iterations in the fixed-point iteration
 * @param[in]       ws          ODE workspace for storing intermediate results
 * @return Result for the convergence of the fixed-point iteration (0 converged, 1 failed)
 * 
 * TODO: Implement a residual-based convergence criterion
 */
static int implicit_midpoint_update(double* w, int w_size, double t, double dt, ode_rhs ode_rhs, 
    struct ode_params* ode_params, double tol, int max_iter, struct ode_ws* ws)
{
    // ODE workspace setup
    ws_check(ws, w_size);
    double* f0      = ws_vec(ws, 0);
    double* f_mid   = ws_vec(ws, 1);
    double* w_mid   = ws_vec(ws, 2);
    double* w_guess = ws_vec(ws, 3);

    // Predictor: explicit Euler
    ode_rhs(t, w, ode_params, f0);
    for (int i = 0; i < w_size; ++i)
        w_guess[i] = w[i] + dt * f0[i];

    const double t_mid = t + 0.5 * dt;
    int failure = 1;

    // Fixed-point iteration on the implicit equation
    for (int it = 0; it < max_iter; ++it) {
        for (int i = 0; i < w_size; ++i)
            w_mid[i] = 0.5 * (w[i] + w_guess[i]);

        ode_rhs(t_mid, w_mid, ode_params, f_mid);

        double max_delta = 0.0;
        double max_scale = 1.0;
        for (int i = 0; i < w_size; ++i) {
            double y_next = w[i] + dt * f_mid[i];

            // Update max_delta (max component-wise change) and max_scale (max component magnitude)
            double delta  = fabs(y_next - w_guess[i]);
            if (delta > max_delta) max_delta = delta;
            double scale = fabs(y_next);
            if (scale > max_scale) max_scale = scale;

            // Update current iteration
            w_guess[i] = y_next;
        }

        // Convergence criterion: stop when update is small compared to the size of the solution
        // If w is very small this becomes: || w_k+1 - w_k || <= tol (absolute criterion)
        // If w is large this becomes: || w_k+1 - w_k || / || w_k+1 || <= tol (relative criterion)
        if (max_delta <= tol * (1.0 + max_scale)) {
            failure = 0;
            break;
        }
    }

    // Write converged result in the final output and return whether the method has converged
    for (int i = 0; i < w_size; ++i)
        w[i] = w_guess[i];

    return failure;
}


/**
 * @brief Performs a timestep with the Cash-Karp method (embedded 5th-order Runge-Kutta method).
 *
 * Performs a timestep with the Cash-Karp method (embedded 5th-order Runge-Kutta method) from a 
 * state w to a new state w_new, according to the ODE w'(t) = ode_rhs(t, w). The result is computed
 * using a fifth-order Runge-Kutta method and the error is estimated with the difference to the
 * result that uses a fourth-order Runge-Kutta method. For this, the RK4 result is not explicitly
 * computed but the differences in the coefficients are used to compute the error estimate.
 * 
 * @param[in]   w           Input state
 * @param[out]  w_new       Updated state
 * @param[out]  w_err       Estimated error
 * @param[in]   w_size      Number of values in state w
 * @param[in]   t           Current time
 * @param[in]   dt          Timestep
 * @param[in]   ode_rhs     Pointer to the function that specifies the ODE right-hand side
 * @param[in]   ode_params  Parameter struct containing general information about the system
 * @param[in]   k1          Pointer to k1/dt (doesn't need to be computed multiple times)
 * @param[in]   ws          ODE workspace for storing intermediate results
 */
static void cash_karp_update(double* w, double* w_new, double* w_err, int w_size, double t, 
    double dt, ode_rhs ode_rhs, struct ode_params* ode_params, double* k1, struct ode_ws* ws)
{
    // ODE workspace setup
    ws_check(ws, w_size);
    double *k2     = ws_vec(ws, 3);
    double *k3     = ws_vec(ws, 4);
    double *k4     = ws_vec(ws, 5);
    double *k5     = ws_vec(ws, 6);
    double *k6     = ws_vec(ws, 7);
    double *w_temp = ws_vec(ws, 8);

    // Coefficients for the Cash-Karp method
    static const double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
                        b21 = 0.2,
                        b31 = 0.075, b32 = 0.225, 
                        b41 = 0.3, b42 = -0.9, b43 = 1.2,
                        b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0, b54 = 35.0/27.0,
                        b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0, 
                        b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
                        c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    
    // Differences in the c-coefficients for computing the error estimate
                        dc1 = -277.0/64512.0, dc3 = 6925.0/370944.0, dc4 = -6925.0/202752.0,
                        dc5 = -277.0/14336.0, dc6 = 277.0/7084.0;

    // Compute k2/dt
    for (int i = 0; i < w_size; i++)
        w_temp[i] = w[i] + b21 * dt * k1[i];
    ode_rhs(t + a2 * dt, w_temp, ode_params, k2);

    // Compute k3/dt
    for (int i = 0; i < w_size; i++)
        w_temp[i] = w[i] + dt * (b31 * k1[i] + b32 * k2[i]);
    ode_rhs(t + a3 * dt, w_temp, ode_params, k3);

    // Compute k4/dt
    for (int i = 0; i < w_size; i++)
        w_temp[i] = w[i] + dt * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
    ode_rhs(t + a4 * dt, w_temp, ode_params, k4);

    // Compute k5/dt
    for (int i = 0; i < w_size; i++)
        w_temp[i] = w[i] + dt * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
    ode_rhs(t + a5 * dt, w_temp, ode_params, k5);

    // Compute k6/dt
    for (int i = 0; i < w_size; i++)
        w_temp[i] = w[i] + dt * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
    ode_rhs(t + a6 * dt, w_temp, ode_params, k6);

    // Update y (c2 and c5 are zero for RK5)
    for (int i = 0; i < w_size; i++)
        w_new[i] = w[i] + dt * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i]);
    
    // Compute the error estimates (c2 is also zero for RK4)
    for (int i = 0; i < w_size; i++)
        w_err[i] = dt * (dc1 * k1[i] + dc3 * k3[i] + dc4 * k4[i] + dc5 * k5[i] + dc6 * k6[i]);
}


/**
 * @brief Driver for the Cash-Karp method. Performs a timestep, ensuring a specified maximum error.
 *
 * Driver for the Cash-Karp method (embedded 5th-order Runge-Kutta method). It performs a timestep,
 * updating the state w according to the ODE w'(t) = ode_rhs(t, w). The timestepping is adaptive -
 * the timestep gets reduced until the estimated error from the Cash-Karp method is below a 
 * specified error tolerance. In case the estimated error is below the error tolerance, the
 * timestep gets increased for the next call. The design and tuning constants are largely inspired
 * by the book "Numerical Recipes in C".
 * 
 * @param[in,out]   w           State (gets updated)
 * @param[in]       w_size      Number of values in state w
 * @param[in,out]   t           Time (gets updated)
 * @param[in,out]   dt          Timestep (gets updated)
 * @param[in]       dt_max      Maximum timestep
 * @param[in]       rtol        Error tolerance
 * @param[in]       ode_rhs     Pointer to the function that specifies the ODE right-hand side
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[in]       ws          ODE workspace for storing intermediate results
 * @return Outcome of the Cash-Karp step (0 success, 1 failed)
 */
static int cash_karp_driver(double* w, int w_size, double* t, double *dt, double dt_max, 
    double rtol, ode_rhs ode_rhs, struct ode_params* ode_params, struct ode_ws *ws)
{
    if (!(rtol > 0.0)) return 1;

    // Tuning constants
    const double SAFETY = 0.9;
    const double PGROW  = -0.2;
    const double PSHRNK = -0.25;
    const double ERRCON = 1.89e-4;   // (5/SAFETY)^(1/PGROW) for max. factor 5 timestep increase
    const double TINY   = 1e-30;
    const int MAX_REJECT = 100;

    // ODE workspace setup
    ws_check(ws, w_size);
    double *k1     = ws_vec(ws, 0);
    double *w_err  = ws_vec(ws, 1);
    double *w_temp = ws_vec(ws, 2);

    // Compute k1/dx (the same for each time step, only needs to be computed once)
    ode_rhs(*t, w, ode_params, k1);

    // Compute step as long as the estimated error is larger than the specified error threshold
    double dt_temp;
    for (int it = 0; it < MAX_REJECT; it++) {

        // Perform a Cash-Karp step and get the updated values w and the error estimates w_err
        cash_karp_update(w, w_temp, w_err, w_size, *t, *dt, ode_rhs, ode_params, k1, ws);

        // Compute scaled error estimate
        double max_error = 0.0;
        for (int i = 0; i < w_size; i++) {
            double denom = fabs(w[i]) + fabs(*dt * k1[i]) + TINY;
            max_error = fmax(max_error, fabs(w_err[i]) / denom);
        }

        max_error /= rtol;

        // If estimated error below rtol, accept the step
        if (max_error <= 1.0) {
            for (int i = 0; i < w_size; i++)
                w[i] = w_temp[i];

            *t += *dt;

            // Increase next step (maximum factor 5 increase)
            if (max_error > ERRCON)
                dt_temp = SAFETY * (*dt) * pow(max_error, PGROW);
            else
                dt_temp = 5.0 * (*dt);
            *dt = fmin(dt_temp, dt_max);
            return 0;
        }   

        // Otherwise, reject step and retry with smaller step size (maximum factor 10 decrease)
        dt_temp = SAFETY * (*dt) * pow(max_error, PSHRNK);
        *dt = fmax(dt_temp, 0.1 * (*dt));
        
        // Check if the step size is too small (underflow)
        if (*t + *dt == *t) return 1;
    }
    return 1;
}


/**
 * @brief Numerically integrates the ODE w'(t) = ode_rhs(t, w) and writes output to files.
 *
 * Numerically integrates the ODE w'(t) = ode_rhs(t, w). The final time, timestep, error tolerances
 * and integration method are read from the user-specified parameters in the parameter file. If the
 * user-specified integration method cannot be found it uses a 4th-order Runge-Kutta method. Output
 * is written to files at user-specified intervals.
 * 
 * @param[in,out]   w           State vector (on input initial state, on output final state)
 * @param[in]       ode_rhs     Pointer to the function that specifies the ODE right-hand side
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * 
 * TODO: Add interpolator for writing outputs at exactly the right times
 */
void ode_integrator(double* w, ode_rhs ode_rhs, struct ode_params* ode_params)
{
    // Determine integration method 
    ode_method m = M_UNKNOWN;
    char *method = get_parameter_string("ode_integrator");
    if (strcmp(method, "cash-karp") == 0) m = M_CASH_KARP;
    else if (strcmp(method, "rk4") == 0) m = M_RK4;
    else if (strcmp(method, "implicit-midpoint") == 0) m = M_IMPLICIT_MIDPOINT;
    if (m == M_UNKNOWN) 
        fprintf(stderr, "Warning: The specified method %s is not implemented; using RK4.\n", 
            method);

    // Load the specified parameters from the parameter database
    check_integration_parameter_validity();
    double rtol = 0.0, tol = 0.0, dt_max = 0.0;
    int max_iter = 0;
    double t_end = get_parameter_double("t_end");
    double dt = get_parameter_double("dt");
    double dt_save = get_parameter_double("dt_save");
    if (!is_set_double(dt_save)) dt_save = dt;
    if (m == M_CASH_KARP) {
        rtol = get_parameter_double("rtol");
        dt_max = get_parameter_double("dt_max");
        if (!is_set_double(dt_max)) dt_max = dt;
    }
    if (m == M_IMPLICIT_MIDPOINT) { 
        tol = get_parameter_double("tol"); 
        max_iter = get_parameter_int("max_iter"); 
    }

    int w_size = 2 * ode_params->num_dim * ode_params->num_bodies;
    int failure = 0;

    double t_current = 0.0;
    long long save_idx = 1;
    double next_save = dt_save;
    const double eps_time = 1e-12 * fmax(1.0, fabs(dt_save));

    // ODE workspace setup
    struct ode_ws ws = {0};
    ode_ws_init(&ws, w_size);

    // Initialize output files
    FILE *file_mass, *file_pos, *file_mom, *file_energy;
    output_init(&file_mass, &file_pos, &file_mom, &file_energy, ode_params);
    output_write_timestep(file_pos, file_mom, file_energy, ode_params, w, t_current);

    // Iterate until the final time
    while (t_current < t_end) {
    
        // Set step size (so we don't step past t_end)
        if (t_current + dt > t_end) dt = t_end - t_current;

        // Use specified ODE integration method
        switch (m) {
            case M_CASH_KARP:
                // The Cash-Karp driver updates t_current and dt internally
                failure = cash_karp_driver(w, w_size, &t_current, &dt, dt_max, rtol, ode_rhs,
                    ode_params, &ws);

                // Use RK4 if the Cash-Karp method is stuck
                if (failure) {
                    progress_bar_break_line();
                    fprintf(stderr, "Warning: The Cash-Karp method is stuck at " 
                        "t = %.7g, using RK4 for this timestep instead.\n", t_current);
                    if (t_current + dt == t_current)
                        errorexit("Stepsize underflow in ode_integrator");
                    rk4_update(w, w_size, t_current, dt, ode_rhs, ode_params, &ws);
                    t_current += dt;
                }
                break;

            case M_RK4:
                rk4_update(w, w_size, t_current, dt, ode_rhs, ode_params, &ws);
                t_current += dt;
                break;

            case M_IMPLICIT_MIDPOINT:
                failure = implicit_midpoint_update(w, w_size, t_current, dt, ode_rhs,
                    ode_params, tol, max_iter, &ws);
            
                // Use RK4 if the fixed-point iteration did not converge
                if (failure) {
                    progress_bar_break_line();
                    fprintf(stderr, "Warning: The implicit midpoint method did not converge at "
                        "t = %.7g, using RK4 for this timestep instead.\n", t_current);
                    rk4_update(w, w_size, t_current, dt, ode_rhs, ode_params, &ws);
                }
                t_current += dt;
                break;

            default:
                rk4_update(w, w_size, t_current, dt, ode_rhs, ode_params, &ws);
                t_current += dt;
                break;
        }

        // Write output if the current time is an output time
        if (t_current + eps_time >= next_save) {
            output_write_timestep(file_pos, file_mom, file_energy, ode_params, w, t_current);
            print_progress_bar((int)(100.0 * t_current / t_end));

            // Advance output schedule robustly
            while (t_current + eps_time >= next_save) {
                save_idx++;
                next_save = save_idx * dt_save;
            }
        }
    }

    // Cleanup
    ode_ws_free(&ws);
    fclose(file_mass);
    fclose(file_pos);
    fclose(file_mom);
    fclose(file_energy);
}


// ------------------------------------------------------------------------------------------------
// Impulse method functions (multiple time stepping)
// ------------------------------------------------------------------------------------------------


// Apply a kick generated by a potential U(x): p <- p - h * dU/dx, x remains unchanged
static void impulse_apply_kick(double* w, int num_bodies, int num_dim, double h, 
    const double* dUdx)
{
    int array_half = num_bodies * num_dim;
    for (int i = 0; i < array_half; ++i)
        w[array_half + i] -= h * dUdx[i];
}


/**
 * @brief Advances the state w for timestep h with middle method of the impulse method
 *
 * Advances the state w for timestep h with middle method of the impulse method. The middle method
 * can be any of the regular ODE integration methods. In case the method is an adaptive method, the
 * integrator integrates adaptively until t advances by h. For non-adaptive methods the integrator
 * performs n fixed substeps of size h/n.
 * 
 * @param[in,out]   w           State vector (on input initial state, on output final state)
 * @param[in]       w_size      Number of values in state w
 * @param[in]       t_start     Start time
 * @param[in]       h           Outer timestep of the impulse method
 * @param[in]       n           Number of substeps for the middle method (if non-adaptive method)
 * @param[in]       m           Numerical ODE integration method used as the middle method
 * @param[in]       err_tol     Error tolerance for numerical ODE integration method
 * @param[in]       max_iter    Maximum number of iterations for numerical ODE integration method
 * @param[in]       rhs_mid     Pointer to the ODE right-hand side corresponding to H0 without UTT4
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[in]       ws          ODE workspace for storing intermediate results
 */
static void impulse_advance_middle(double* w, int w_size, double t_start, double h, int n, 
    ode_method m, double err_tol, int max_iter, ode_rhs rhs_mid, struct ode_params* ode_params, 
    struct ode_ws* ws)
{
    ws_check(ws, w_size);
    int failure = 0;
    double dt = h/n;
    double t = t_start;
    double t_end = t_start + h;

    // Use specified ODE integration method
    switch (m) {
        case M_RK4:
            for (int k = 0; k < n; ++k) {
                rk4_update(w, w_size, t, dt, rhs_mid, ode_params, ws);
                t += dt;
            }
            return;

        case M_CASH_KARP:
            while (t < t_end) {
                double remaining = t_end - t;
                if (remaining <= 0.0) break;

                // Ensure we don't step past the end of the middle interval.
                if (dt > remaining) dt = remaining;

                // Limit max step to remaining as well (keeps the driver well-behaved).
                double dt_max = remaining;

                failure = cash_karp_driver(w, w_size, &t, &dt, dt_max, err_tol, rhs_mid, 
                    ode_params, ws);

                // Use RK4 if the Cash-Karp method fails
                if (failure) {
                    progress_bar_break_line();
                    fprintf(stderr, "Warning: The Cash-Karp failed at " 
                        "t = %.7g, using RK4 for this timestep instead.\n", t);
                    dt = h/n;
                    if (dt > remaining) dt = remaining;
                    if (t + dt == t) {
                        errorexit("Stepsize underflow in impulse_advance_middle");
                    }
                    rk4_update(w, w_size, t, dt, rhs_mid, ode_params, ws);
                    t += dt; 
                }
            }
            return;

        case M_IMPLICIT_MIDPOINT:
            for (int k = 0; k < n; ++k) {
                failure = implicit_midpoint_update(w, w_size, t, dt, rhs_mid, ode_params, err_tol,
                    max_iter, ws);
                
                // Use RK4 if the fixed-point iteration did not converge
                if (failure) {
                    progress_bar_break_line();
                    fprintf(stderr, "Warning: The implicit midpoint method did not converge at "
                        "t = %.7g, using RK4 for this timestep instead.\n", t);
                    rk4_update(w, w_size, t, dt, rhs_mid, ode_params, ws);
                }
                t += dt;
            }
            return;

        default:
            for (int k = 0; k < n; ++k) {
                rk4_update(w, w_size, t, dt, rhs_mid, ode_params, ws);
                t += dt;
            }
            return;
    }
}


/**
 * @brief Integrates the ODE w'(t) = ode_rhs(t, w) with the impulse method and writes the output.
 *
 * Numerically integrates the ODE w'(t) = ode_rhs(t, w) using the impulse method:
 * \Psi_h = \Phi^(TT4)_{h/2} o (\Phi^(0)_{h/n})^n o \Phi^(TT4)_{h/2},
 * where n is the number of substeps for the method associated with H0, h is the outer timestep of
 * the full composition method \Psi_h, h/n is the inner (effective) timestep, and \Phi^(0)_{h/n} 
 * and \Phi^(TT4)_{h/2} are numerical integrators consistent with the Hamiltonians H0 and UTT4, 
 * respectively (see Heinze, Schäfer and Brügmann 2026).
 * The evaluation of the gradient of UTT4 is very expensive and therefore \Phi^(TT4)_{h/2} is
 * performed only twice per outer timestep (grad_utt4 only needs to be computed once per outer 
 * timestep). The final time, timestep, error tolerances and integration method are read from the 
 * user-specified parameters in the parameter file. If the user-specified integration method cannot
 * be found it uses a 4th-order Runge-Kutta method. Output is written to files at user-specified 
 * intervals.
 * 
 * @param[in,out]   w           State vector (on input initial state, on output final state)
 * @param[in]       rhs_mid     Pointer to the ODE right-hand side corresponding to H0 without UTT4
 * @param[in]       grad_utt4   Pointer to the function that computes the gradient of UTT4
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * 
 * TODO: Add interpolator for writing outputs at exactly the right times
 */
void ode_integrator_impulse(double* w, ode_rhs rhs_mid, utt4_grad_func grad_utt4, 
    struct ode_params* params)
{
    // Determine integration method 
    ode_method m = M_UNKNOWN;
    char *middle_method = get_parameter_string("ode_integrator");
    if (strcmp(middle_method, "cash-karp") == 0) m = M_CASH_KARP;
    else if (strcmp(middle_method, "rk4") == 0) m = M_RK4;
    else if (strcmp(middle_method, "implicit-midpoint") == 0) m = M_IMPLICIT_MIDPOINT;
    if (m == M_UNKNOWN)
        fprintf(stderr, "Warning: The specified method %s is not implemented; using RK4.\n", 
            middle_method);
    
    // Load the specified parameters from the parameter database
    check_integration_parameter_validity();
    double t_end     = get_parameter_double("t_end");
    double dt_full   = get_parameter_double("dt");
    double dt_save   = get_parameter_double("dt_save");
    if (!is_set_double(dt_save)) dt_save = dt_full;
    int n = get_parameter_int("impulse_method_n");
    double err_tol = 0.0;
    int max_iter = 0;
    if (m == M_CASH_KARP) 
        err_tol = get_parameter_double("rtol");
    if (m == M_IMPLICIT_MIDPOINT) { 
        err_tol = get_parameter_double("tol"); 
        max_iter = get_parameter_int("max_iter"); 
    }

    int num_bodies = params->num_bodies;
    int num_dim    = params->num_dim;
    int array_half = num_bodies * num_dim;
    int w_size     = 2 * array_half;

    double t_current = 0.0;
    long long save_idx = 1;
    double next_save = dt_save;
    const double eps_time = 1e-12 * fmax(1.0, fabs(dt_save));

    // ODE workspace setup
    struct ode_ws ws = {0};
    ode_ws_init(&ws, w_size);

    // Initialize output files
    FILE *file_mass, *file_pos, *file_mom, *file_energy;
    output_init(&file_mass, &file_pos, &file_mom, &file_energy, params);
    output_write_timestep(file_pos, file_mom, file_energy, params, w, t_current);

    // Cache gradient to reuse between steps:
    // After finishing a step, positions do not change before the next step's first half-kick,
    // so the end gradient is also the next step's start gradient.
    double* dUdx = NULL;
    allocate_vector(&dUdx, array_half);
    int grad_valid = 0;

    // Iterate until the final time
    while (t_current < t_end) {
        double h = dt_full;
        if (t_current + h > t_end) h = t_end - t_current;
        if (h <= 0.0) break;

        // First half kick (uses cached gradient if available)
        if (!grad_valid) {
            grad_utt4(w, params, dUdx);
            grad_valid = 1;
        }
        impulse_apply_kick(w, num_bodies, num_dim, 0.5 * h, dUdx);

        // Middle evolution for timestep h
        impulse_advance_middle(w, w_size, t_current, h, n, m, err_tol, max_iter,
            rhs_mid, params, &ws);
        t_current += h;

        // Compute gradient at end positions, do second half kick, keep for reuse
        grad_utt4(w, params, dUdx);
        impulse_apply_kick(w, num_bodies, num_dim, 0.5 * h, dUdx);
        grad_valid = 1;

        // Write output if the current time is an output time
        if (t_current + eps_time >= next_save) {
            output_write_timestep(file_pos, file_mom, file_energy, params, w, t_current);
            print_progress_bar((int)(100.0 * t_current / t_end));

            // Advance output schedule robustly
            while (t_current + eps_time >= next_save) {
                save_idx++;
                next_save = save_idx * dt_save;
            }
        }
    }

    // Cleanup
    fclose(file_mass);
    fclose(file_pos);
    fclose(file_mom);
    fclose(file_energy);
    free_vector(dUdx);
    ode_ws_free(&ws);
}
