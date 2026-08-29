/**
 * @file pn_binary.c
 * @brief Relativistic PN binary utilities.
 */

#include "pn_binary.h"

#include <float.h>
#include <math.h>
#include <stddef.h>


/**
 * @brief Reduced conservative ADM Hamiltonian of an isolated nonspinning binary.
 *
 * The dimensionless variables are x = r/M, pr_hat = p_r/mu and j = J/(mu M), where M is the total
 * mass and mu the reduced mass. The returned value is H_rel/mu and always includes the Newtonian
 * term.
 */
double pn_binary_reduced_hamiltonian(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn)
{
    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;
    const double inv_x3 = inv_x2 * inv_x;

    const double pr2 = pr_hat * pr_hat;
    const double j2 = j * j;
    const double p2 = pr2 + j2 * inv_x2;

    const double p4 = p2 * p2;
    const double p6 = p4 * p2;

    double h = 0.5 * p2 - inv_x;

    if (use_1pn) {
        h += 0.125 * (3.0 * nu - 1.0) * p4
           - 0.5 * inv_x * ((3.0 + nu) * p2 + nu * pr2)
           + 0.5 * inv_x2;
    }

    if (use_2pn) {
        h += (1.0 - 5.0 * nu + 5.0 * nu * nu) * p6 / 16.0

           + inv_x / 8.0
             * ((5.0 - 20.0 * nu - 3.0 * nu * nu) * p4
                - 2.0 * nu * nu * pr2 * p2
                - 3.0 * nu * nu * pr2 * pr2)

           + 0.5 * inv_x2 * ((5.0 + 8.0 * nu) * p2 + 3.0 * nu * pr2)

           - 0.25 * (1.0 + 3.0 * nu) * inv_x3;
    }

    return h;
}


/** @brief Analytic derivative of the reduced Hamiltonian with respect to x = r/M. */
double pn_binary_reduced_hamiltonian_dx(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn)
{
    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;
    const double inv_x3 = inv_x2 * inv_x;
    const double pr2 = pr_hat * pr_hat;
    const double j2 = j * j;
    const double p2 = pr2 + j2 * inv_x2;

    /* Differentiate with respect to u = 1/x first, then use du/dx = -u^2. */
    double dh_du = j2 * inv_x - 1.0;

    if (use_1pn) {
        dh_du += 0.5 * (3.0 * nu - 1.0) * j2 * inv_x * p2
               - 0.5 * ((3.0 + nu) * p2 + nu * pr2)
               - (3.0 + nu) * j2 * inv_x2
               + inv_x;
    }

    if (use_2pn) {
        const double nu2 = nu * nu;
        const double p4 = p2 * p2;
        const double p4_coefficient = 5.0 - 20.0 * nu - 3.0 * nu2;
        const double radial_polynomial =
            p4_coefficient * p4 - 2.0 * nu2 * pr2 * p2 - 3.0 * nu2 * pr2 * pr2;
        const double potential_polynomial = (5.0 + 8.0 * nu) * p2 + 3.0 * nu * pr2;

        dh_du += 3.0 * (1.0 - 5.0 * nu + 5.0 * nu2)
                     * j2 * inv_x * p4 / 8.0
               + radial_polynomial / 8.0
               + 0.5 * j2 * inv_x2 * (p4_coefficient * p2 - nu2 * pr2)
               + inv_x * potential_polynomial
               + (5.0 + 8.0 * nu) * j2 * inv_x3
               - 0.75 * (1.0 + 3.0 * nu) * inv_x2;
    }

    return -inv_x2 * dh_du;
}


/** @brief Analytic derivative of the reduced Hamiltonian with respect to p_r/mu. */
double pn_binary_reduced_hamiltonian_dpr(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn)
{
    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;
    const double pr2 = pr_hat * pr_hat;
    const double p2 = pr2 + j * j * inv_x2;
    double coefficient = 1.0;

    if (use_1pn)
        coefficient += 0.5 * (3.0 * nu - 1.0) * p2
                     - (3.0 + 2.0 * nu) * inv_x;

    if (use_2pn) {
        const double nu2 = nu * nu;
        const double p4_coefficient = 5.0 - 20.0 * nu - 3.0 * nu2;
        coefficient += 3.0 * (1.0 - 5.0 * nu + 5.0 * nu2) * p2 * p2 / 8.0
                     + 0.5 * inv_x
                         * ((p4_coefficient - nu2) * p2 - 4.0 * nu2 * pr2)
                     + (5.0 + 11.0 * nu) * inv_x2;
    }

    return pr_hat * coefficient;
}


/**
 * @brief Analytic derivative of the reduced Hamiltonian with respect to reduced angular momentum.
 *
 * For H = mu*h(x, pr_hat, j), Hamilton's equation gives the physical azimuthal frequency
 * Omega = dH/dJ = (1/M)*dh/dj. Returning dh/dj keeps this helper in the dimensionless conventions
 * used throughout this file.
 */
double pn_binary_reduced_hamiltonian_dj(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn)
{
    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;
    const double pr2 = pr_hat * pr_hat;
    const double p2 = pr2 + j*j*inv_x2;

    double dh_dp2 = 0.5;

    if (use_1pn) {
        dh_dp2 += 0.25*(3.0*nu - 1.0)*p2
                  - 0.5*inv_x*(3.0 + nu);
    }

    if (use_2pn) {
        const double nu2 = nu*nu;
        dh_dp2 += 3.0*(1.0 - 5.0*nu + 5.0*nu2)*p2*p2/16.0
                  + inv_x/8.0
                    *(2.0*(5.0 - 20.0*nu - 3.0*nu2)*p2 - 2.0*nu2*pr2)
                  + 0.5*inv_x2*(5.0 + 8.0*nu);
    }

    return 2.0*j*inv_x2*dh_dp2;
}


// Evaluate the reduced Hamiltonian at a radial turning point (pr_hat = 0).
double pn_binary_turning_hamiltonian(double x, double j, double nu,
    int use_1pn, int use_2pn)
{
    return pn_binary_reduced_hamiltonian(x, 0.0, j, nu, use_1pn, use_2pn);
}


typedef double (*ScalarResidual)(double variable, const void *context);


static int same_nonzero_sign(double a, double b)
{
    return (a < 0.0 && b < 0.0) || (a > 0.0 && b > 0.0);
}


// Shared safeguarded bisection used by every isolated-binary scalar solve.
static int bisect_scalar_root(ScalarResidual residual, const void *context,
    double x_a, double f_a, double x_b, double f_b, double *root)
{
    if (root == NULL || isnan(f_a) || isnan(f_b) || same_nonzero_sign(f_a, f_b))
        return 0;
    if (f_a == 0.0) {
        *root = x_a;
        return 1;
    }
    if (f_b == 0.0) {
        *root = x_b;
        return 1;
    }

    for (int iteration = 0; iteration < 192; ++iteration) {
        const double x_mid = x_a + 0.5 * (x_b - x_a);
        if (x_mid == x_a || x_mid == x_b)
            break;

        const double f_mid = residual(x_mid, context);
        if (isnan(f_mid))
            return 0;
        if (same_nonzero_sign(f_a, f_mid)) {
            x_a = x_mid;
            f_a = f_mid;
        }
        else {
            x_b = x_mid;
            f_b = f_mid;
        }
    }

    (void)f_b;
    *root = x_a + 0.5 * (x_b - x_a);
    return isfinite(*root);
}


// Bracket a nonnegative root by repeatedly doubling the upper endpoint.
static int solve_nonnegative_root(ScalarResidual residual, const void *context,
    double initial_upper, int max_expansions, double *root)
{
    double lower = 0.0;
    double f_lower = residual(lower, context);
    double upper = initial_upper;
    double f_upper = residual(upper, context);

    if (isnan(f_lower) || isnan(f_upper))
        return 0;

    for (int expansion = 0;
         same_nonzero_sign(f_lower, f_upper) && expansion < max_expansions;
         ++expansion) {
        if (!(upper <= 0.5 * DBL_MAX))
            return 0;
        upper *= 2.0;
        f_upper = residual(upper, context);
        if (isnan(f_upper))
            return 0;
    }

    return bisect_scalar_root(residual, context,
        lower, f_lower, upper, f_upper, root);
}


typedef struct AngularMomentumProblem {
    double x_inner;
    double x_outer;
    double nu;
    int use_1pn;
    int use_2pn;
} AngularMomentumProblem;


// Radial derivative of the turning-point Hamiltonian at fixed angular momentum.
static double pn_binary_dh_dx_turning_numeric(double x, double j, double nu,
    int use_1pn, int use_2pn)
{
    double dx = 1e-6 * fabs(x);

    if (dx < 1e-8)
        dx = 1e-8;
    if (x - dx <= 0.0)
        dx = 0.5 * x;

    const double hp = pn_binary_turning_hamiltonian(x + dx, j, nu, use_1pn, use_2pn);
    const double hm = pn_binary_turning_hamiltonian(x - dx, j, nu, use_1pn, use_2pn);
    return (hp - hm) / (2.0 * dx);
}


static double circular_angular_momentum_residual(double j, const void *context)
{
    const AngularMomentumProblem *problem = context;
    return pn_binary_dh_dx_turning_numeric(problem->x_inner, j, problem->nu,
        problem->use_1pn, problem->use_2pn);
}


static double eccentric_angular_momentum_residual(double j, const void *context)
{
    const AngularMomentumProblem *problem = context;
    return pn_binary_turning_hamiltonian(problem->x_inner, j, problem->nu,
               problem->use_1pn, problem->use_2pn)
         - pn_binary_turning_hamiltonian(problem->x_outer, j, problem->nu,
               problem->use_1pn, problem->use_2pn);
}


/**
 * @brief Solve dH(x, p_r=0, j)/dx = 0 for the circular reduced angular momentum.
 * 
 * Returns NAN if a nonnegative root cannot be bracketed.
 */
double pn_binary_solve_circular_j(double x, double nu, int use_1pn, int use_2pn)
{
    if (!isfinite(x) || !(x > 0.0) || !isfinite(nu) || !(nu > 0.0) || nu > 0.25)
        return NAN;

    const AngularMomentumProblem problem = {
        .x_inner = x,
        .x_outer = x,
        .nu = nu,
        .use_1pn = use_1pn,
        .use_2pn = use_2pn
    };
    double j;
    if (!solve_nonnegative_root(circular_angular_momentum_residual, &problem,
            2.0 * sqrt(x) + 1.0, 80, &j))
        return NAN;
    return j;
}


/**
 * @brief Solve H(xp, 0, j) = H(xa, 0, j) for the eccentric reduced angular momentum.
 * 
 * Returns NAN if a nonnegative root cannot be bracketed.
 */
double pn_binary_solve_eccentric_j(double xp, double xa, double nu,
    int use_1pn, int use_2pn)
{
    if (!isfinite(xp) || !(xp > 0.0) || !isfinite(xa) || !(xa > xp)
        || !isfinite(nu) || !(nu > 0.0) || nu > 0.25)
        return NAN;

    const AngularMomentumProblem problem = {
        .x_inner = xp,
        .x_outer = xa,
        .nu = nu,
        .use_1pn = use_1pn,
        .use_2pn = use_2pn
    };
    double j;
    if (!solve_nonnegative_root(eccentric_angular_momentum_residual, &problem,
            2.0 * sqrt(0.5 * (xp + xa)) + 1.0, 80, &j))
        return NAN;
    return j;
}


typedef struct RadialMomentumProblem {
    double x;
    double j;
    double energy;
    double nu;
    int use_1pn;
    int use_2pn;
} RadialMomentumProblem;


static double radial_momentum_squared_residual(double pr_hat_squared, const void *context)
{
    const RadialMomentumProblem *problem = context;
    return pn_binary_reduced_hamiltonian(problem->x, sqrt(pr_hat_squared), problem->j,
        problem->nu, problem->use_1pn, problem->use_2pn) - problem->energy;
}


/**
 * @brief Solve H(x, p_r, j) = reduced_energy for the nonnegative absolute radial momentum.
 * 
 * Returns NAN if no real root can be bracketed.
 */
double pn_binary_solve_pr_hat_abs(double x, double j, double reduced_energy, double nu,
    int use_1pn, int use_2pn)
{
    if (!isfinite(x) || !(x > 0.0) || !isfinite(j) || !(j >= 0.0)
        || !isfinite(reduced_energy) || !isfinite(nu) || !(nu > 0.0) || nu > 0.25)
        return NAN;

    const RadialMomentumProblem problem = {
        .x = x,
        .j = j,
        .energy = reduced_energy,
        .nu = nu,
        .use_1pn = use_1pn,
        .use_2pn = use_2pn
    };
    const double residual_at_zero = radial_momentum_squared_residual(0.0, &problem);
    if (fabs(residual_at_zero) < 1e-13)
        return 0.0;
    if (residual_at_zero > 0.0)
        return residual_at_zero < 1e-10 ? 0.0 : NAN;

    double pr_hat_squared;
    if (!solve_nonnegative_root(radial_momentum_squared_residual, &problem,
            1e-12, 120, &pr_hat_squared))
        return NAN;
    return sqrt(pr_hat_squared);
}


// Return A in dH/dpr_hat = A pr_hat + O(pr_hat^3) at pr_hat = 0.
double pn_binary_dh_dpr_coeff_at_zero(double x, double j, double nu,
    int use_1pn, int use_2pn)
{
    if (!isfinite(x) || !(x > 0.0) || !isfinite(j) || !(j >= 0.0)
        || !isfinite(nu) || !(nu > 0.0) || nu > 0.25)
        return NAN;

    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;
    const double j2 = j * j;
    const double tangential_momentum_squared = j2 * inv_x2;
    double coefficient = 1.0;

    if (use_1pn) {
        coefficient += 0.5 * (3.0 * nu - 1.0) * tangential_momentum_squared
                     - (3.0 + 2.0 * nu) * inv_x;
    }

    if (use_2pn) {
        const double p6_coefficient =
            (1.0 - 5.0 * nu + 5.0 * nu * nu) / 16.0;
        const double p4_coefficient = 5.0 - 20.0 * nu - 3.0 * nu * nu;

        coefficient += 6.0 * p6_coefficient
                         * tangential_momentum_squared * tangential_momentum_squared
                     + inv_x / 8.0
                         * (4.0 * p4_coefficient * tangential_momentum_squared
                            - 4.0 * nu * nu * tangential_momentum_squared)
                     + inv_x2 * (5.0 + 11.0 * nu);
    }

    return coefficient;
}


typedef struct TurningPointProblem {
    double energy;
    double j;
    double nu;
    int use_1pn;
    int use_2pn;
} TurningPointProblem;


static double turning_residual(double x, const void *context)
{
    const TurningPointProblem *problem = context;
    return pn_binary_turning_hamiltonian(x, problem->j, problem->nu,
        problem->use_1pn, problem->use_2pn) - problem->energy;
}


static double turning_residual_scale(double x, const TurningPointProblem *problem)
{
    const double inv_x = 1.0 / x;
    return fmax(fabs(problem->energy),
        fmax(inv_x, 0.5 * problem->j * problem->j * inv_x * inv_x));
}


static int find_enclosing_turning_root(double current_x, int direction,
    const TurningPointProblem *problem, double *root)
{
    double x_inside = current_x;
    double f_inside = turning_residual(x_inside, problem);
    const double initial_tolerance =
        1024.0 * DBL_EPSILON * turning_residual_scale(current_x, problem);

    if (fabs(f_inside) <= initial_tolerance)
        f_inside = 0.0;
    if (f_inside > 0.0)
        return 0;

    double relative_step = 1e-8;
    for (int iteration = 0; iteration < 512; ++iteration) {
        double x_outside;
        if (direction < 0)
            x_outside = x_inside / (1.0 + relative_step);
        else
            x_outside = x_inside * (1.0 + relative_step);

        if (!(x_outside > 0.0) || !isfinite(x_outside) || x_outside == x_inside)
            return 0;

        double f_outside = turning_residual(x_outside, problem);
        if (!isfinite(f_outside)) {
            /* Overflow at very small radius is a usable positive barrier only. */
            if (direction < 0 && signbit(f_outside))
                return 0;
        }

        const double tolerance =
            1024.0 * DBL_EPSILON * turning_residual_scale(x_outside, problem);
        if (isfinite(f_outside) && fabs(f_outside) <= tolerance)
            f_outside = 0.0;

        if (f_outside > 0.0 || (isinf(f_outside) && !signbit(f_outside))) {
            return bisect_scalar_root(turning_residual, problem,
                x_inside, f_inside, x_outside, f_outside, root)
                && *root > 0.0;
        }

        x_inside = x_outside;
        f_inside = f_outside;
        relative_step = fmin(0.25, 1.6 * relative_step);
    }

    return 0;
}


/**
 * @brief Solve for the two positive radial turning points enclosing current_x.
 *
 * The inputs use the same dimensionless conventions as pn_binary_reduced_hamiltonian(). 
 * Only bound states with two finite positive turning points succeed. On success, 
 * x_pericenter <= current_x <= x_apocenter.
 *
 * @return 1 on success and 0 when no physical enclosing pair can be bracketed.
 */
int pn_binary_solve_turning_points(double reduced_energy, double j, double nu,
    int use_1pn, int use_2pn, double current_x,
    double *x_pericenter, double *x_apocenter)
{
    if (x_pericenter == NULL || x_apocenter == NULL)
        return 0;

    *x_pericenter = NAN;
    *x_apocenter = NAN;

    if (!isfinite(reduced_energy) || !(reduced_energy < 0.0)
        || !isfinite(j) || !(j > 0.0)
        || !isfinite(nu) || !(nu > 0.0) || nu > 0.25
        || !isfinite(current_x) || !(current_x > 0.0))
        return 0;

    const TurningPointProblem problem = {
        .energy = reduced_energy,
        .j = j,
        .nu = nu,
        .use_1pn = use_1pn,
        .use_2pn = use_2pn
    };

    const double current_residual = turning_residual(current_x, &problem);
    const double current_scale = turning_residual_scale(current_x, &problem);
    if (fabs(current_residual) <= 1024.0 * DBL_EPSILON * current_scale) {
        double derivative_step = 1e-5 * current_x;
        if (derivative_step < 1e-10)
            derivative_step = 1e-10;
        if (derivative_step < current_x) {
            const double derivative =
                (turning_residual(current_x + derivative_step, &problem)
                 - turning_residual(current_x - derivative_step, &problem))
                / (2.0 * derivative_step);
            if (isfinite(derivative)
                && fabs(derivative) * current_x <= 1e-8 * current_scale) {
                *x_pericenter = current_x;
                *x_apocenter = current_x;
                return 1;
            }
        }
    }

    double inner_root;
    double outer_root;
    if (!find_enclosing_turning_root(current_x, -1, &problem, &inner_root)
        || !find_enclosing_turning_root(current_x, +1, &problem, &outer_root))
        return 0;

    if (inner_root > outer_root) {
        const double temporary = inner_root;
        inner_root = outer_root;
        outer_root = temporary;
    }

    const double enclosure_tolerance =
        64.0 * DBL_EPSILON * fmax(current_x, outer_root);
    if (current_x < inner_root - enclosure_tolerance
        || current_x > outer_root + enclosure_tolerance)
        return 0;

    *x_pericenter = inner_root;
    *x_apocenter = outer_root;
    return 1;
}
