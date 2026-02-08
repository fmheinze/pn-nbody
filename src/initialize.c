/**
 * @file initialize.c
 * @brief Functions for initializing parameters and initial condition presets
 *
 * Functions for the initialization of the parameter database, the ODE parameters, as well as the
 * state vector based on specific initial condition presets. The parameters are checked for 
 * consistency and validity and missing parameters are computed from user-specified ones.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "initialize.h"
#include "eom.h"
#include "parameters.h"
#include "utils.h"
#include "initial_configurations.h"

#define NUM_PN_TERMS 4


/**
 * @brief Initializes needed parameters in the parameter database.
 * 
 * Initializes needed parameters in the parameter database, together with default values 
 * (-1 usually means not set, 1 on, 0 off), as well as a short description that contains which 
 * values are valid.
 */
void initialize_parameters() 
{
    // --------------------------------------------------------------------------------------------
    // General parameters
    // --------------------------------------------------------------------------------------------
    add_parameter("num_dim", "3", "number of dimensions [2, 3]");
    add_parameter("num_bodies", "0", "number of bodies [1, 2, 3, ...]");
    add_parameter("0pn_terms", "1", "whether to include 0PN (Newtonian) terms [0, 1]");
    add_parameter("1pn_terms", "0", "whether to include 1PN terms [0, 1]");
    add_parameter("2pn_terms", "0", "whether to include 2PN terms [0, 1]");
    add_parameter("2.5pn_terms", "0", "whether to include 2.5PN terms [0, 1]");

    // --------------------------------------------------------------------------------------------
    // Numerical parameters
    // --------------------------------------------------------------------------------------------
    add_parameter("t_end", "0.0", "total duration of the simulation [>= 0]");
    add_parameter("dt", "0.1", "time step / initial time step for adaptive ODE integrators [> 0]");
    add_parameter("ode_integrator", "rk4", "which ODE integrator to use [rk4, cash-karp, ...]");
    add_parameter("impulse_method", "0", "whether to use the impulse method");

    // Cash-Karp method
    if(strcmp(get_parameter_string("ode_integrator"), "cash-karp") == 0) {
        add_parameter("rtol", "-1", "target relative error tolerance [> 0]");
        add_parameter("dt_max", "-1", "target relative error tolerance [> 0]");
    }

    // Implicit-midpoint
    if(strcmp(get_parameter_string("ode_integrator"), "implicit-midpoint") == 0) {
        add_parameter("tol", "-1", "target error tolerance [> 0]");
        add_parameter("max_iter", "50", "maximum number of fixed-point interations [> 0]");
    }

    // Impulse method
    if(get_parameter_int("impulse_method") == 1) {
        add_parameter("impulse_method_n", "1", "number of substeps for the impulse method");
    }

    // --------------------------------------------------------------------------------------------
    // Output parameters
    // --------------------------------------------------------------------------------------------
    add_parameter("dt_save", "-1", "times at which quantities are written to a file [>= 0]");

    // --------------------------------------------------------------------------------------------
    // Initial configuration presets
    // --------------------------------------------------------------------------------------------
    add_parameter("ic_preset", "-1", "specify initial condition preset");

    // Newtonian binary
    if(strcmp(get_parameter_string("ic_preset"), "newtonian_binary") == 0) {
        add_parameter("binary_a", "-1", "semi-major axis [> 0]");
        add_parameter("binary_b", "-1", "semi-minor axis [> 0]");
        add_parameter("binary_e", "-1", "eccentricity [>= 0]");
        add_parameter("binary_ra", "-1", "apoapsis [> 0]");
        add_parameter("binary_rp", "-1", "periapsis [> 0]");
        add_parameter("binary_p", "-1", "semi-parameter [> 0]");
        add_parameter("binary_phi0", "0.0", "initial phase");
    }

    // Newtonian binary-single scattering
    if(strcmp(get_parameter_string("ic_preset"), "binary_single_scattering") == 0) {
        add_parameter("binary_a", "-1", "semi-major axis [> 0]");
        add_parameter("binary_b", "-1", "semi-minor axis [> 0]");
        add_parameter("binary_e", "-1", "eccentricity [>= 0]");
        add_parameter("binary_ra", "-1", "apoapsis [> 0]");
        add_parameter("binary_rp", "-1", "periapsis [> 0]");
        add_parameter("binary_p", "-1", "semi-parameter [> 0]");
        add_parameter("binary_phi0", "0.0", "initial phase");
        add_parameter("d0", "-1", "initial distance of the binary and the single [> 0]");
        add_parameter("v0_rel", "-1", "initial relative approach velocity [>= 0]");
        add_parameter("b", "-1", "scattering impact parameter [>= 0]");
        add_parameter("orientation", "-1", "components of the orientation of the binary");
    }

    // Relativistic binary-single scattering
    if(strcmp(get_parameter_string("ic_preset"), "binary_single_scattering_rel") == 0) {
        add_parameter("d0", "-1", "initial distance of the binary and the single [> 0]");
        add_parameter("p0_rel", "-1", "initial relative approach momentum [>= 0]");
        add_parameter("b", "-1", "scattering impact parameter [>= 0]");
        add_parameter("orientation", "-1", "components of the orientation of the binary");
        add_parameter("binary_r0", "-1", "initial binary radius [> 0]");
        add_parameter("binary_pt0", "-1", "tangential initial momentum [>= 0]");
        add_parameter("binary_pr0", "-1", "radial initial momentum [>= 0]");
        add_parameter("binary_phi0", "0.0", "initial phase of the binary [>= 0]");
    }

    // Newtonian binary-binary scattering
    if(strcmp(get_parameter_string("ic_preset"), "binary_binary_scattering") == 0) {
        add_parameter("d0", "-1", "initial distance of the binary and the single [> 0]");
        add_parameter("v0_rel", "-1", "initial relative approach velocity [>= 0]");
        add_parameter("b", "-1", "scattering impact parameter [>= 0]");
        add_parameter("binary1_a", "-1", "semi-major axis of binary 1 [> 0]");
        add_parameter("binary1_b", "-1", "semi-minor axis of binary 1 [> 0]");
        add_parameter("binary1_e", "-1", "eccentricity of binary 1 [>= 0]");
        add_parameter("binary1_ra", "-1", "apoapsis of binary 1 [> 0]");
        add_parameter("binary1_rp", "-1", "periapsis of binary 1 [> 0]");
        add_parameter("binary1_p", "-1", "semi-parameter of binary 1 [> 0]");
        add_parameter("binary1_phi0", "0.0", "initial phase of binary 1");
        add_parameter("orientation_1", "-1", "components of the orientation of binary 1");
        add_parameter("binary2_a", "-1", "semi-major axis of binary 2 [> 0]");
        add_parameter("binary2_b", "-1", "semi-minor axis of binary 2 [> 0]");
        add_parameter("binary2_e", "-1", "eccentricity of binary 2 [>= 0]");
        add_parameter("binary2_ra", "-1", "apoapsis of binary 2 [> 0]");
        add_parameter("binary2_rp", "-1", "periapsis of binary 2 [> 0]");
        add_parameter("binary2_p", "-1", "semi-parameter of binary 2 [> 0]");
        add_parameter("binary2_phi0", "0.0", "initial phase of binary 2");
        add_parameter("orientation_2", "-1", "components of the orientation of binary 2");
    }

    // Relativistic binary-binary scattering
    if(strcmp(get_parameter_string("ic_preset"), "binary_binary_scattering_rel") == 0) {
        add_parameter("d0", "-1", "initial distance of the binary and the single [> 0]");
        add_parameter("p0_rel", "-1", "initial relative approach momentum [>= 0]");
        add_parameter("b", "-1", "scattering impact parameter [>= 0]");
        add_parameter("orientation_1", "-1", "components of the orientation of binary 1");
        add_parameter("binary1_r0", "-1", "initial radius of binary 1 [> 0]");
        add_parameter("binary1_pt0", "-1", "tangential initial momentum for binary 1 [>= 0]");
        add_parameter("binary1_pr0", "-1", "radial initial momentum for binary 1 [>= 0]");
        add_parameter("binary1_phi0", "0.0", "initial phase for binary 1 [>= 0]");
        add_parameter("orientation_2", "-1", "components of the orientation of binary 2");
        add_parameter("binary2_r0", "-1", "initial radius of binary 2 [> 0]");
        add_parameter("binary2_pt0", "-1", "tangential initial momentum for binary 2 [>= 0]");
        add_parameter("binary2_pr0", "-1", "radial initial momentum for binary 2 [>= 0]");
        add_parameter("binary2_phi0", "0.0", "initial phase for binary 2 [>= 0]");
    }

    // Figure-eight orbit
    if(strcmp(get_parameter_string("ic_preset"), "figure_eight") == 0) {
        add_parameter("figure_eight_width", "1000.0", "width of the figure-eight orbit [> 0]");
    }

    // --------------------------------------------------------------------------------------------
    // Parameters associated with individual bodies
    // --------------------------------------------------------------------------------------------
    int num_bodies = get_parameter_int("num_bodies"); 
    for(int i = 1; i < num_bodies+1; i++) {
        add_parameter_i("mass", i, "0.0", "mass of body i [>= 0]");
        add_parameter_i("pos", i, "0.0 0.0 0.0", "coordinates of the initial position of body i");
        add_parameter_i("p", i, "0.0 0.0 0.0", "components of the initial momentum of body i");
    }
}


/**
 * @brief Construct and return an ode_params struct from the parameter database.
 *
 * Loads all global simulation/ODE configuration needed by the integrator and the ODE
 * right-hand side (RHS), and caches them in an ode_params struct to avoid repeated
 * parameter-database lookups during RHS evaluations. The returned struct contains heap-allocated 
 * arrays (params.masses and params.pn_terms) and the caller is responsible for freeing these.
 *
 * @return Fully initialized ode_params with cached configuration and allocated arrays.
 */
struct ode_params initialize_ode_params()
{
    struct ode_params params;

    // General parameters
    params.num_dim = get_parameter_int("num_dim");
    params.num_bodies = get_parameter_int("num_bodies");
    params.use_impulse_method = get_parameter_int("impulse_method");

    // Masses
    allocate_vector(&params.masses, params.num_bodies);
    for (int i = 0; i < params.num_bodies; i++)
        params.masses[i] = get_parameter_double_i("mass", i+1);

    // PN terms
    params.pn_terms = (int *)malloc(NUM_PN_TERMS * sizeof(double));
    if (params.pn_terms == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    params.pn_terms[0] = get_parameter_int("0pn_terms");
    params.pn_terms[1] = get_parameter_int("1pn_terms");
    params.pn_terms[2] = get_parameter_int("2pn_terms");
    params.pn_terms[3] = get_parameter_int("2.5pn_terms");

    // Check validity
    if (params.num_dim != 2 && params.num_dim != 3) 
        errorexit("Please specify a valid num_dim (must be 2 or 3)");
    if (params.num_bodies <= 0) 
        errorexit("Please specify a valid num_bodies (must be num_bodies > 0)");
    if (params.use_impulse_method != 0 && params.use_impulse_method != 1) 
        errorexit("Please set impulse_method to 0 (off) or 1 (on)");
    for (int i = 0; i < params.num_bodies; i++) {
        if (params.masses[i] < 0) 
            errorexit("Please specify valid masses (must be mass >= 0)");
    }
    for (int i = 0; i < NUM_PN_TERMS; i++) {
        if (params.pn_terms[i] != 0 && params.pn_terms[i] != 1) 
            errorexit("Please set pn_terms to 0 (off) or 1 (on)");
    }

    return params;
}


/**
 * @brief Load, complete, and validate orbital parameters for a binary.
 *
 * Reads the orbital parameter set for binary i from the parameter database and returns a
 * fully-initialized binary_params struct. Unspecified values must be set to -1 in the database.
 * The function loads all supported orbital parameters, validates any user-specified values,
 * infers missing parameters from the specified ones, and verifies that the resulting set is 
 * self-consistent. At least two orbital parameters (in addition to phi0) must be specified to 
 * determine a unique orbit, otherwise the function aborts with an error.
 * 
 * @param[in]   i   Binary index (use 0 for the default/unindexed "binary_*" parameter set)
 * @return A fully initialized binary_params struct
 */
struct binary_params initialize_binary_params(int i)
{
    struct binary_params params;

    // Load specified parameters from the parameter database (-1 if not specified)
    params.a    = get_binary_parameter_double_i("a", i);
    params.b    = get_binary_parameter_double_i("b", i);
    params.e    = get_binary_parameter_double_i("e", i);
    params.r_a  = get_binary_parameter_double_i("ra", i);
    params.r_p  = get_binary_parameter_double_i("rp", i);
    params.p    = get_binary_parameter_double_i("p", i);
    params.phi0 = get_binary_parameter_double_i("phi0", i);

    // Check the user-specified values
    if (is_set_double(params.a) && params.a <= 0.0) 
        errorexit("Please specify a valid a (must be > 0)");
    if (is_set_double(params.e) && (params.e < 0.0 || params.e >= 1.0)) 
        errorexit("Please specify a valid e (must be 0 <= e < 1 for an elliptical orbit)");
    if (is_set_double(params.b) && params.b <= 0.0) 
        errorexit("Please specify a valid b (must be > 0)");
    if (is_set_double(params.a) && is_set_double(params.b) && params.a < params.b) 
        errorexit("Invalid a and b (must be a >= b)");
    if (is_set_double(params.r_p) && params.r_p <= 0.0) 
        errorexit("Please specify a valid r_p (must be > 0)");
    if (is_set_double(params.r_a) && params.r_a <= 0.0)
        errorexit("Please specify a valid r_a (must be > 0)");
    if (is_set_double(params.r_a) && is_set_double(params.r_p) && params.r_a < params.r_p) 
        errorexit("Invalid r_a and r_p (must be r_a >= r_p)");
    if (is_set_double(params.p) && params.p <= 0.0) 
        errorexit("Please specify a valid p (must be > 0)");

    // Iteratively infer missing values from whatever is known
    int changed = 1;
    for (int iter = 0; iter < 20 && changed; ++iter) {
        changed = 0;

        // --- From (a, e) ---
        if (is_set_double(params.a) && params.a > 0.0 && 
            is_set_double(params.e) && params.e >= 0.0 && params.e < 1.0)
        {
            double one_minus_e2 = clamp0(1.0 - params.e*params.e);
            changed |= set_if_unset_double(&params.b,   params.a * sqrt(one_minus_e2));
            changed |= set_if_unset_double(&params.r_p, params.a * (1.0 - params.e));
            changed |= set_if_unset_double(&params.r_a, params.a * (1.0 + params.e));
            changed |= set_if_unset_double(&params.p,   params.a * one_minus_e2);
        }

        // --- From (a, b) ---
        if (is_set_double(params.a) && params.a > 0.0 && 
            is_set_double(params.b) && params.b > 0.0 && params.b <= params.a)
        {
            double e2 = clamp0(1.0 - (params.b*params.b)/(params.a*params.a));
            changed |= set_if_unset_double(&params.e, sqrt(e2));
            changed |= set_if_unset_double(&params.p, (params.b*params.b)/params.a);
        }

        // --- From (a, p) ---
        if (is_set_double(params.a) && params.a > 0.0 && 
            is_set_double(params.p) && params.p > 0.0 && params.p <= params.a)
        {
            double e2 = clamp0(1.0 - params.p/params.a);
            changed |= set_if_unset_double(&params.e, sqrt(e2));
            changed |= set_if_unset_double(&params.b, sqrt(params.a*params.p));
            if (is_set_double(params.e) && params.e < 1.0) {
                double e_now = params.e;
                changed |= set_if_unset_double(&params.r_p, params.p/(1.0 + e_now));
                changed |= set_if_unset_double(&params.r_a, params.p/(1.0 - e_now));
            }
        }

        // --- From (r_a, r_p) ---
        if (is_set_double(params.r_a) && params.r_a > 0.0 && is_set_double(params.r_p) && 
            params.r_p > 0.0 && params.r_a >= params.r_p)
        {
            double aa = 0.5 * (params.r_a + params.r_p);
            double ee = (params.r_a - params.r_p) / (params.r_a + params.r_p);
            changed |= set_if_unset_double(&params.a, aa);
            changed |= set_if_unset_double(&params.e, ee);
            changed |= set_if_unset_double(&params.p, (2.0 * params.r_a * params.r_p) / 
                (params.r_a + params.r_p));
            changed |= set_if_unset_double(&params.b, sqrt(params.r_a * params.r_p));
        }

        // --- From (r_p, e) ---
        if (is_set_double(params.r_p) && params.r_p > 0.0 && is_set_double(params.e) && 
            params.e >= 0.0 && params.e < 1.0)
        {
            double aa = params.r_p / (1.0 - params.e);
            changed |= set_if_unset_double(&params.a, aa);
            changed |= set_if_unset_double(&params.r_a, aa * (1.0 + params.e));
            changed |= set_if_unset_double(&params.p, params.r_p * (1.0 + params.e));
        }

        // --- From (r_a, e) ---
        if (is_set_double(params.r_a) && params.r_a > 0.0 && is_set_double(params.e) &&
            params.e >= 0.0 && params.e < 1.0)
        {
            double aa = params.r_a / (1.0 + params.e);
            changed |= set_if_unset_double(&params.a, aa);
            changed |= set_if_unset_double(&params.r_p, aa * (1.0 - params.e));
            changed |= set_if_unset_double(&params.p, params.r_a * (1.0 - params.e));
        }

        // --- From (p, e) ---
        if (is_set_double(params.p) && params.p > 0.0 && is_set_double(params.e) &&
            params.e >= 0.0 && params.e < 1.0)
        {
            double denom = clamp0(1.0 - params.e*params.e);
            if (denom > 0.0) {
                double aa = params.p / denom;
                changed |= set_if_unset_double(&params.a, aa);
                changed |= set_if_unset_double(&params.b, params.p / sqrt(denom));
                changed |= set_if_unset_double(&params.r_p, params.p/(1.0 + params.e));
                changed |= set_if_unset_double(&params.r_a, params.p/(1.0 - params.e));
            }
        }

        // --- From (p, r_p) ---
        if (is_set_double(params.p) && params.p > 0.0 && is_set_double(params.r_p) &&
            params.r_p > 0.0)
        {
            double ee = params.p/params.r_p - 1.0;
            if (ee >= 0.0 && ee < 1.0) {
                changed |= set_if_unset_double(&params.e, ee);
            }
        }

        // --- From (p, r_a) ---
        if (is_set_double(params.p) && params.p > 0.0 && is_set_double(params.r_a) &&
            params.r_a > 0.0)
        {
            double ee = 1.0 - params.p/params.r_a;
            if (ee >= 0.0 && ee < 1.0) {
                changed |= set_if_unset_double(&params.e, ee);
            }
        }

        // --- From (b, p) ---
        if (is_set_double(params.b) && params.b > 0.0 && is_set_double(params.p) && 
            params.p > 0.0)
        {
            double aa = (params.b*params.b)/params.p;
            double e2 = clamp0(1.0 - (params.p*params.p)/(params.b*params.b));
            changed |= set_if_unset_double(&params.a, aa);
            changed |= set_if_unset_double(&params.e, sqrt(e2));
        }

        // --- From (b, e) ---
        if (is_set_double(params.b) && params.b > 0.0 && is_set_double(params.e) &&
            params.e >= 0.0 && params.e < 1.0)
        {
            double denom = clamp0(1.0 - params.e*params.e);
            if (denom > 0.0) {
                double aa = params.b / sqrt(denom);
                changed |= set_if_unset_double(&params.a, aa);
                changed |= set_if_unset_double(&params.p, aa * denom);
            }
        }
    }

    // At this point, if the user gave enough info, we expect everything to be set
    if (!is_set_double(params.a) || !is_set_double(params.b) || !is_set_double(params.e) ||
        !is_set_double(params.r_a) || !is_set_double(params.r_p) || !is_set_double(params.p)) {
        errorexit("Insufficient information to determine the binary parameters "
            "(need a consistent pair like (a, e), (r_a, r_p), (p, e), etc.)");
    }

    // Consistency checks (useful if user specified more than two values)
    const double rel_eps = 1e-9;
    {
        double b_expected  = params.a * sqrt(clamp0(1.0 - params.e*params.e));
        double rp_expected = params.a * (1.0 - params.e);
        double ra_expected = params.a * (1.0 + params.e);
        double p_expected  = params.a * (1.0 - params.e*params.e);

        if (!almost_equal(params.b,   b_expected,  rel_eps)) 
            errorexit("Inconsistent parameters (b doesn't match a,e)");
        if (!almost_equal(params.r_p, rp_expected, rel_eps)) 
            errorexit("Inconsistent parameters (r_p doesn't match a,e)");
        if (!almost_equal(params.r_a, ra_expected, rel_eps)) 
            errorexit("Inconsistent parameters (r_a doesn't match a,e)");
        if (!almost_equal(params.p,   p_expected,  rel_eps)) 
            errorexit("Inconsistent parameters (p doesn't match a,e)");
    }

    // Print the final binary parameters
    if (i == 0) {
        printf("Full list of binary parameters:\n");
        printf("binary_a\t= %10.16e\n", params.a);
        printf("binary_b\t= %10.16e\n", params.b);
        printf("binary_e\t= %10.16e\n", params.e);
        printf("binary_ra\t= %10.16e\n", params.r_a);
        printf("binary_rp\t= %10.16e\n", params.r_p);
        printf("binary_p\t= %10.16e\n", params.p);
        printf("binary_phi0\t= %10.16e\n", params.phi0);
    } else {
        printf("Full list of parameters for binary %d:\n", i);
        printf("binary%d_a\t= %10.16e\n", i, params.a);
        printf("binary%d_b\t= %10.16e\n", i, params.b);
        printf("binary%d_e\t= %10.16e\n", i, params.e);
        printf("binary%d_ra\t= %10.16e\n", i, params.r_a);
        printf("binary%d_rp\t= %10.16e\n", i, params.r_p);
        printf("binary%d_p\t= %10.16e\n", i, params.p);
        printf("binary%d_phi0\t= %10.16e\n", i, params.phi0);
    }
    print_divider();

    return params;
}


/**
 * @brief Build and return the initial state vector w0.
 *
 * Allocates and initializes the state vector w0 = [x_1, ..., x_N, p_1, ..., p_N] with 
 * N = ode_params->num_bodies and each x_i, p_i having dimension D = ode_params->num_dim.
 * The initialization is done either from an initial-condition preset (ic_preset), or from 
 * user-specified per-body arrays ("pos" and "p") if no preset is selected.
 *
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @return Pointer to a newly allocated state vector of length 2 * N * D.
 *         The caller owns the memory and must free it with the appropriate deallocator.
 */
double* initialize_state_vector(struct ode_params* ode_params)
{
    int num_dim = ode_params->num_dim;
    int num_bodies = ode_params->num_bodies;

    // Allocate state vector
    double* w0;
    allocate_vector(&w0, 2 * num_dim * num_bodies);

    // Check whether an initial condition preset has been selected and initialize w0 accordinly
    const char *preset = get_parameter_string("ic_preset");
    if (strcmp(preset, "-1") != 0) {
        if (strcmp(preset, "newtonian_binary") == 0)
            initialize_newtonian_binary(ode_params, w0);

        else if (strcmp(preset, "binary_single_scattering") == 0)
            initialize_binary_single_scattering(ode_params, w0);

        else if (strcmp(preset, "binary_single_scattering_rel") == 0)
            initialize_binary_single_scattering_rel(ode_params, w0);

        else if (strcmp(preset, "binary_binary_scattering") == 0)
            initialize_binary_binary_scattering(ode_params, w0);
        
        else if (strcmp(preset, "binary_binary_scattering_rel") == 0)
            initialize_binary_binary_scattering_rel(ode_params, w0);

        else if (strcmp(preset, "figure_eight") == 0)
            initialize_figure_eight(ode_params, w0);

        else errorexit("Unknown ic_preset");

        // Print the computed initial positions and momenta
        printf("Computed initial positions and momenta from the initial condition preset:\n");
        print_state_vector(w0, num_bodies, num_dim);
        print_divider();
    }
    // If no initial condition preset has been selected, use specified positions and momenta
    else {
        for (int i = 0; i < num_bodies; i++) {
            double* pos = get_parameter_double_array_i("pos", i+1); 
            double* p = get_parameter_double_array_i("p", i+1); 
            for (int j = 0; j < num_dim; j++) {
                w0[i * num_dim + j] = pos[j];
                w0[num_dim * num_bodies + i * num_dim + j] = p[j];
            }
            free_vector(pos);
            free_vector(p);
        }  
    }
    return w0;
}


/**
 * @brief Initializes the state vector for a Newtonian binary
 * 
 * Initializes the state vector for a Newtonian binary and checks the consistency and validity of
 * the involved parameters.
 * 
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[out]  w0              Initialized state vector, w0 = [positions, momenta]
 */
void initialize_newtonian_binary(struct ode_params* ode_params, double* w0)
{
    printf("Setting up initial parameters for a Newtonian binary...\n");
    print_divider();
    struct binary_params binary_params = initialize_binary_params(0);
    ic_newtonian_binary(ode_params, &binary_params, w0);
}


/**
 * @brief Initializes the state vector for a Newtonian binary-single scattering
 * 
 * Initializes the state vector for a Newtonian binary-single scattering and checks the 
 * consistency and validity of the involved parameters.
 * 
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[out]  w0              Initialized state vector, w0 = [positions, momenta]
 */
void initialize_binary_single_scattering(struct ode_params* ode_params, double* w0)
{
    printf("Setting up initial parameters for a Newtonian binary-single scattering...\n");
    print_divider();

    // Load specified values
    struct binary_params binary_params = initialize_binary_params(0);
    double d0 = get_parameter_double("d0");
    double v0_rel = get_parameter_double("v0_rel");
    double b = get_parameter_double("b");
    double* orientation;
    if (strcmp(get_parameter_string("orientation"), "-1") == 0)
        orientation = NULL;
    else {
        allocate_vector(&orientation, 3);
        for (int i = 0; i < 3; i++)
            orientation[i] = get_parameter_double_array("orientation")[i];
    }

    // Check specified values for validity
    if (d0 < 0) errorexit("Please specify a valid d0 (must be d0 >= 0)");
    if (v0_rel < 0) errorexit("Please specify a valid v0_rel (must be v0_rel >= 0)");
    if (fabs(b) > d0) errorexit("Please specify a valid b (must be |b| <= d0)");

    ic_binary_single_scattering(ode_params, &binary_params, d0, v0_rel, b, orientation, w0);
    free_vector(orientation);
}


/**
 * @brief Initializes the state vector for a relativistic binary-single scattering
 * 
 * Initializes the state vector for a relativistic binary-single scattering and checks the 
 * consistency and validity of the involved parameters.
 * 
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[out]  w0              Initialized state vector, w0 = [positions, momenta]
 */
void initialize_binary_single_scattering_rel(struct ode_params* ode_params, double* w0)
{
    printf("Setting up initial parameters for a relativistic binary-single scattering...\n");
    print_divider();

    // Load specified values
    double d0 = get_parameter_double("d0");
    double p0_rel = get_parameter_double("p0_rel");
    double b = get_parameter_double("b");
    double binary_r0 = get_parameter_double("binary_r0");
    double binary_pt0 = get_parameter_double("binary_pt0");
    double binary_pr0 = get_parameter_double("binary_pr0");
    double binary_phi0 = get_parameter_double("binary_phi0");
    double* orientation;
    if (strcmp(get_parameter_string("orientation"), "-1") == 0)
        orientation = NULL;
    else {
        allocate_vector(&orientation, 3);
        for (int i = 0; i < 3; i++)
            orientation[i] = get_parameter_double_array("orientation")[i];
    }

    // Check specified values for validity
    if (d0 < 0) errorexit("Please specify a valid d0 (must be d0 >= 0)");
    if (p0_rel < 0) errorexit("Please specify a valid p0_rel (must be v0_rel >= 0)");
    if (fabs(b) > d0) errorexit("Please specify a valid b (must be |b| <= d0)");
    if (binary_r0 <= 0) errorexit("Please specify a valid binary_r0 (must be binary_r0 > 0)");
    if (binary_pt0 < 0) errorexit("Please specify a valid binary_pt0 (must be binary_pt0 >= 0)");
    if (binary_pr0 < 0) errorexit("Please specify a valid binary_pr0 (must be binary_pr0 >= 0)");
    if (binary_phi0 < 0) errorexit("Please specify a valid binary_phi0 (must be binary_phi0 >= 0)");
    if (ode_params->masses[0] != ode_params->masses[1])
        errorexit("Currently only equal-mass binaries supported in rel. binary-single scattering");

    ic_binary_single_scattering_rel(d0, p0_rel, b, binary_phi0, binary_r0, binary_pt0, 
        binary_pr0, orientation, w0);

    free_vector(orientation);
}


/**
 * @brief Initializes the state vector for a Newtonian binary-binary scattering
 * 
 * Initializes the state vector for a Newtonian binary-binary scattering and checks the 
 * consistency and validity of the involved parameters.
 * 
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[out]  w0              Initialized state vector, w0 = [positions, momenta]
 */
void initialize_binary_binary_scattering(struct ode_params* ode_params, double* w0)
{
    printf("Setting up initial parameters for a Newtonian binary-binary scattering...\n");
    print_divider();

    // Load specified values
    struct binary_params binary1_params = initialize_binary_params(1);
    struct binary_params binary2_params = initialize_binary_params(2);
    double d0 = get_parameter_double("d0");
    double v0_rel = get_parameter_double("v0_rel");
    double b = get_parameter_double("b");
    double* orientation_1;
    if (strcmp(get_parameter_string("orientation_1"), "-1") == 0)
        orientation_1 = NULL;
    else {
        allocate_vector(&orientation_1, 3);
        for (int i = 0; i < 3; i++)
            orientation_1[i] = get_parameter_double_array("orientation_1")[i];
    }
    double* orientation_2;
    if (strcmp(get_parameter_string("orientation_2"), "-1") == 0)
        orientation_2 = NULL;
    else {
        allocate_vector(&orientation_2, 3);
        for (int i = 0; i < 3; i++)
            orientation_2[i] = get_parameter_double_array("orientation_2")[i];
    }

    // Check specified values for validity
    if (d0 < 0) errorexit("Please specify a valid d0 (d0 >= 0)");
    if (v0_rel < 0) errorexit("Please specify a valid v0_rel (v0_rel >= 0)");
    if (fabs(b) > d0) errorexit("Please specify a valid b (|b| <= d0)");

    ic_binary_binary_scattering(ode_params, &binary1_params, &binary2_params, d0, v0_rel, b,
        orientation_1, orientation_2, w0);

    free_vector(orientation_1);
    free_vector(orientation_2);
}


/**
 * @brief Initializes the state vector for a relativistic binary-binary scattering
 * 
 * Initializes the state vector for a relativistic binary-binary scattering and checks the 
 * consistency and validity of the involved parameters.
 * 
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[out]  w0              Initialized state vector, w0 = [positions, momenta]
 */
void initialize_binary_binary_scattering_rel(struct ode_params* ode_params, double* w0)
{
    printf("Setting up initial parameters for a relativistic binary-binary scattering...\n");
    print_divider();

    // Load specified values
    double d0 = get_parameter_double("d0");
    double p0_rel = get_parameter_double("p0_rel");
    double b = get_parameter_double("b");
    double binary1_r0 = get_parameter_double("binary1_r0");
    double binary1_pt0 = get_parameter_double("binary1_pt0");
    double binary1_pr0 = get_parameter_double("binary1_pr0");
    double binary1_phi0 = get_parameter_double("binary1_phi0");
    double binary2_r0 = get_parameter_double("binary2_r0");
    double binary2_pt0 = get_parameter_double("binary2_pt0");
    double binary2_pr0 = get_parameter_double("binary2_pr0");
    double binary2_phi0 = get_parameter_double("binary2_phi0");
    double* orientation_1;
    if (strcmp(get_parameter_string("orientation_1"), "-1") == 0)
        orientation_1 = NULL;
    else {
        allocate_vector(&orientation_1, 3);
        for (int i = 0; i < 3; i++)
            orientation_1[i] = get_parameter_double_array("orientation_1")[i];
    }
    double* orientation_2;
    if (strcmp(get_parameter_string("orientation_2"), "-1") == 0)
        orientation_2 = NULL;
    else {
        allocate_vector(&orientation_2, 3);
        for (int i = 0; i < 3; i++)
            orientation_2[i] = get_parameter_double_array("orientation_2")[i];
    }

    // Check specified values for validity
    if (d0 < 0) errorexit("Please specify a valid d0 (must be d0 >= 0)");
    if (p0_rel < 0) errorexit("Please specify a valid p0_rel (must be v0_rel >= 0)");
    if (fabs(b) > d0) errorexit("Please specify a valid b (must be |b| <= d0)");
    if (binary1_r0 <= 0) errorexit("Please specify a valid binary1_r0 (must be binary1_r0 > 0)");
    if (binary1_pt0 < 0) errorexit("Please specify a valid binary1_pt0 (must be binary1_pt0 >= 0)");
    if (binary1_pr0 < 0) errorexit("Please specify a valid binary1_pr0 (must be binary1_pr0 >= 0)");
    if (binary1_phi0 < 0) errorexit("Please specify a valid binary1_phi0 (must be binary1_phi0 >= 0)");
    if (binary2_r0 <= 0) errorexit("Please specify a valid binary2_r0 (must be binary2_r0 > 0)");
    if (binary2_pt0 < 0) errorexit("Please specify a valid binary2_pt0 (must be binary2_pt0 >= 0)");
    if (binary2_pr0 < 0) errorexit("Please specify a valid binary2_pr0 (must be binary2_pr0 >= 0)");
    if (binary2_phi0 < 0) errorexit("Please specify a valid binary2_phi0 (must be binary2_phi0 >= 0)");
    if (ode_params->masses[0] != ode_params->masses[1] || 
        ode_params->masses[2] != ode_params->masses[3])
        errorexit("Currently only equal-mass binaries supported in rel. binary-binary scattering");

    ic_binary_binary_scattering_rel(d0, p0_rel, b, binary1_phi0, binary1_r0, binary1_pt0,
        binary1_pr0, orientation_1, binary2_phi0, binary2_r0, binary2_pt0, binary2_pr0, 
        orientation_2, w0);

    free_vector(orientation_1);
    free_vector(orientation_2);
}


/**
 * @brief Initializes the state vector for a figure-eight orbit
 * 
 * Initializes the state vector for a figure-eight orbit (both Newtonian and relativistic) and 
 * checks the consistency and validity of the involved parameters.
 * 
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[out]  w0              Initialized state vector, w0 = [positions, momenta]
 */
void initialize_figure_eight(struct ode_params* ode_params, double* w0)
{
    printf("Setting up initial parameters for a figure-eight orbit...\n");
    print_divider();

    // Check specified values for validity
    if (ode_params->masses[0] != ode_params->masses[1] || 
        ode_params->masses[0] != ode_params->masses[2] ||
        ode_params->masses[1] != ode_params->masses[2])
        errorexit("The figure-eight orbit requires all three masses to be equal!");

    double width = get_parameter_double("figure_eight_width");
    if ((width > 10000 || width < 100) && (ode_params->pn_terms[1] == 1 || 
        ode_params->pn_terms[2] == 1))
        printf("Warning: figure_eight_width = %lf, post-Newtonian figure-eight orbit "
            "only accurate for 100 < width < 100000!\n", width);

    ic_figure_eight_orbit(ode_params, width, w0);
}
