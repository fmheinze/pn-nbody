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
#include "pn_eom.h"
#include "parameters.h"
#include "utils.h"
#include "initial_configurations.h"

#define NUM_PN_TERMS 4


/**
 * @brief Initializes needed parameters in the parameter database.
 * 
 * Initializes needed parameters in the parameter database, together with default values 
 * (-1 usually means not set), as well as a short description that contains which values
 * are valid.
 */
void initialize_parameters() 
{
    // --------------------------------------------------------------------------------------------
    // General parameters
    // --------------------------------------------------------------------------------------------
    add_parameter("num_dim", "3", "number of dimensions [2, 3]");
    add_parameter("num_bodies", "0", "number of bodies [1, 2, 3, ...]");
    add_parameter("0pn_terms", "1", "whether to include 0PN (Newtonian) terms [0, 1]");
    add_parameter("1pn_terms", "1", "whether to include 1PN terms [0, 1]");
    add_parameter("2pn_terms", "1", "whether to include 2PN terms [0, 1]");
    add_parameter("2.5pn_terms", "1", "whether to include 2.5PN terms [0, 1]");

    // --------------------------------------------------------------------------------------------
    // Numerical parameters
    // --------------------------------------------------------------------------------------------
    add_parameter("t_end", "0.0", "total duration of the simulation [>= 0]");
    add_parameter("dt", "0.1", "time step / initial time step for adaptive ODE integrators [> 0]");
    add_parameter("ode_integrator", "rk4", "which ODE integrator to use [rk4, cash-karp, ...]");
    add_parameter("impulse_method", "0", "whether to use the impulse method");

    // Cash-Karp method specific
    if(strcmp(get_parameter_string("ode_integrator"), "cash-karp") == 0) {
        add_parameter("rel_error", "1e-6", "target relative error [> 0]");
    }

    // Implicit-midpoint method specific
    if(strcmp(get_parameter_string("ode_integrator"), "implicit-midpoint") == 0) {
        add_parameter("rel_error", "1e-6", "target relative error [> 0]");
    }

    // Impulse method specific
    if(get_parameter_int("impulse_method") == 1) {
        add_parameter("impulse_method_n", "1", "number of substeps for the impulse method");
    }

    // --------------------------------------------------------------------------------------------
    // Initial configuration presets
    // --------------------------------------------------------------------------------------------
    add_parameter("ic_preset", "-1", "specify initial condition preset");

    // Parameters for the Newtonian binary
    if(strcmp(get_parameter_string("ic_preset"), "newtonian_binary") == 0) {
        add_parameter("binary_a", "-1", "semi-major axis [> 0]");
        add_parameter("binary_b", "-1", "semi-minor axis [> 0]");
        add_parameter("binary_e", "-1", "eccentricity [>= 0]");
        add_parameter("binary_ra", "-1", "apoapsis [> 0]");
        add_parameter("binary_rp", "-1", "periapsis [> 0]");
        add_parameter("binary_p", "-1", "semi-parameter [> 0]");
        add_parameter("binary_phi0", "0.0", "initial phase");
    }

    // Parameters for the Newtonian binary-single scattering
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

    // Parameters for the relativistic binary-single scattering
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

    // Parameters for the Newtonian binary-single scattering
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

    // Parameters for the relativistic binary-single scattering
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

    // Figure-eight orbit specific
    if(strcmp(get_parameter_string("ic_preset"), "figure_eight") == 0) {
        add_parameter("figure_eight_width", "1000.0", "width of the figure-eight orbit [> 0]");
    }

    // --------------------------------------------------------------------------------------------
    // Parameters associated with individual bodies
    // --------------------------------------------------------------------------------------------
    int num_bodies = get_parameter_int("num_bodies"); 
    for(int i = 1; i < num_bodies+1; i++) {
        add_parameter_i("mass", i, "0.0", "mass of body i");
        add_parameter_i("pos", i, "0.0 0.0 0.0", "coordinates of the initial position of body i");
        add_parameter_i("p", i, "0.0 0.0 0.0", "components of the initial momentum of body i");
    }
}


/**
 * @brief Initializes the parameters of ode_params.
 * 
 * Initializes the parameters of ode_params, with the specified values from the parameter database.
 * These contain general information about the system and specifications for the ODE right-hand 
 * side. These parameters are often used at each ODE right-hand side evaluation and are therefore 
 * cached to avoid repeated costly database lookups.
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
    allocate_vector(&params.pn_terms, NUM_PN_TERMS);
    params.pn_terms[0] = get_parameter_int("0pn_terms");
    params.pn_terms[1] = get_parameter_int("1pn_terms");
    params.pn_terms[2] = get_parameter_int("2pn_terms");
    params.pn_terms[3] = get_parameter_int("2.5pn_terms");

    return params;
}


/**
 * @brief Initializes the state vector w0 based on specified values.
 * 
 * Initializes the state vector w0 based on either the individually specified values for the 
 * initial positions and momenta of each body, or on the specified parameters for an initial 
 * condition preset.
 * 
 * @param[in]   ode_params      Parameter struct containing general information about the system
 * @param[out]  w0              Initialized state vector, w0 = [positions, momenta]
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
    if (preset && preset[0] != '\0') {
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
 * @brief Initializes a binary_params struct based on specified parameters.
 * 
 * Initializes a binary_params struct based on the specified parameters of binary i in the 
 * parameter database. It loads the values from the database, checks them for validity, computes
 * missing values based on the specified ones (at least two binary parameters need to be specified
 * besides phi0), checks for consistency and finally prints all the binary parameters.
 * 
 * @param[in]    i              Index of binary (set to 0 if the system contains only one binary)
 * @param[out]   binary_params  Struct containing the binary parameters
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
        errorexit("Invalid a (must be > 0)");
    if (is_set_double(params.e) && (params.e < 0.0 || params.e >= 1.0)) 
        errorexit("Invalid e (must be 0 <= e < 1 for an elliptical orbit)");
    if (is_set_double(params.b) && params.b <= 0.0) 
        errorexit("Invalid b (must be > 0)");
    if (is_set_double(params.a) && is_set_double(params.b) && params.a < params.b) 
        errorexit("Invalid a and b (must be a >= b)");
    if (is_set_double(params.r_p) && params.r_p <= 0.0) 
        errorexit("Invalid r_p (must be > 0)");
    if (is_set_double(params.r_a) && params.r_a <= 0.0)
        errorexit("Invalid r_a (must be > 0)");
    if (is_set_double(params.r_a) && is_set_double(params.r_p) && params.r_a < params.r_p) 
        errorexit("Invalid r_a and r_p (must be r_a >= r_p)");
    if (is_set_double(params.p) && params.p <= 0.0) 
        errorexit("Invalid p (must be > 0)");

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
    if (d0 < 0) errorexit("Please specify a valid d0 (d0 >= 0)");
    if (v0_rel < 0) errorexit("Please specify a valid v0_rel (v0_rel >= 0)");
    if (fabs(b) > d0) errorexit("Please specify a valid b (|b| <= d0)");

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
    if (d0 < 0) errorexit("Please specify a valid d0 (d0 >= 0)");
    if (p0_rel < 0) errorexit("Please specify a valid p0_rel (v0_rel >= 0)");
    if (fabs(b) > d0) errorexit("Please specify a valid b (|b| <= d0)");
    if (binary_r0 <= 0) errorexit("Please specify a valid binary_r0 (binary_r0 > 0)");
    if (binary_pt0 < 0) errorexit("Please specify a valid binary_pt0 (binary_pt0 >= 0)");
    if (binary_pr0 < 0) errorexit("Please specify a valid binary_pr0 (binary_pr0 >= 0)");
    if (binary_phi0 < 0) errorexit("Please specify a valid binary_phi0 (binary_phi0 >= 0)");
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
    if (d0 < 0) errorexit("Please specify a valid d0 (d0 >= 0)");
    if (p0_rel < 0) errorexit("Please specify a valid p0_rel (v0_rel >= 0)");
    if (fabs(b) > d0) errorexit("Please specify a valid b (|b| <= d0)");
    if (binary1_r0 <= 0) errorexit("Please specify a valid binary1_r0 (binary1_r0 > 0)");
    if (binary1_pt0 < 0) errorexit("Please specify a valid binary1_pt0 (binary1_pt0 >= 0)");
    if (binary1_pr0 < 0) errorexit("Please specify a valid binary1_pr0 (binary1_pr0 >= 0)");
    if (binary1_phi0 < 0) errorexit("Please specify a valid binary1_phi0 (binary1_phi0 >= 0)");
    if (binary2_r0 <= 0) errorexit("Please specify a valid binary2_r0 (binary2_r0 > 0)");
    if (binary2_pt0 < 0) errorexit("Please specify a valid binary2_pt0 (binary2_pt0 >= 0)");
    if (binary2_pr0 < 0) errorexit("Please specify a valid binary2_pr0 (binary2_pr0 >= 0)");
    if (binary2_phi0 < 0) errorexit("Please specify a valid binary2_phi0 (binary2_phi0 >= 0)");
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