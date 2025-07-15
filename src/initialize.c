#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "initialize.h"
#include "pn_eom.h"
#include "parameters.h"
#include "utils.h"
#include "initial_configurations.h"

#define NUM_PN_TERMS 4


void initialize_parameters() {
    // General parameters
    add_parameter("num_dim", "3", "number of dimensions [2, 3]");
    add_parameter("num_bodies", "0", "number of bodies [1, 2, 3, ...]");
    add_parameter("0pn_terms", "1", "whether to include 0PN (Newtonian) terms [0, 1]");
    add_parameter("1pn_terms", "1", "whether to include 1PN terms [0, 1]");
    add_parameter("2pn_terms", "1", "whether to include 2PN terms [0, 1]");
    add_parameter("2.5pn_terms", "1", "whether to include 2.5PN terms [0, 1]");

    // Numerical parameters
    add_parameter("t_end", "0.0", "total duration of the simulation [>= 0]");
    add_parameter("dt", "0.1", "fixed time step, or initial time step for adaptive ODE integrators [> 0]");

    // Cash-Karp method specific
    add_parameter("rel_error", "1e-6", "target relative error [> 0]");

    // Initial configuration presets
    add_parameter("ic_preset", "-1", "specify initial condition preset");
    char *preset = get_parameter_string("ic_preset");
    // Newtonian binary
    if(strcmp(preset, "newtonian_binary") == 0) {
        add_parameter("a", "-1", "semimajor axis [> 0]");
        add_parameter("e", "0.0", "eccentricity [>= 0]");
        add_parameter("phi0", "0.0", "initial phase");
    }

    // Individual masses, positions and momenta for each body
    int num_bodies = get_parameter_int("num_bodies"); 
    for(int i = 1; i < num_bodies+1; i++) {
        add_parameter_i("mass", i, "0.0", "mass of body i");
        add_parameter_i("pos", i, "0.0 0.0 0.0", "x, y and z coordinates of the initial position of body i");
        add_parameter_i("p", i, "0.0 0.0 0.0", "x, y and z components of the initial momentum of body i");
    }

    // Output parameters
    add_parameter("dt_save", "0.1", "time after which outputs are saved [>= dt]");
}


double* initialize_state_vector(struct ode_params* params) {
    // Add part that checks for initial condition presets
    int num_dim = params->num_dim;
    int num_bodies = params->num_bodies;

    double* w;
    allocate_vector(&w, 2 * num_dim * num_bodies);
    // Check whether an initial condition preset has been selected
    if(!get_parameter_int("ic_preset")){
        char *preset = get_parameter_string("ic_preset");
        if(strcmp(preset, "newtonian_binary") == 0) {
            double a = get_parameter_double("a");
            double e = get_parameter_double("e");
            double phi0 = get_parameter_double("phi0");
            newtonian_binary(params, w, a, e, phi0);
        }
    }
    // If no initial condition preset has been selected, use specified positions and momenta
    else {
        for(int i = 0; i < num_bodies; i++){
            double* pos = get_parameter_double_array_i("pos", i+1); 
            double* p = get_parameter_double_array_i("p", i+1); 
            for(int j = 0; j < num_dim; j++){
                w[i * num_dim + j] = pos[j];
                w[num_dim * num_bodies + i * num_dim + j] = p[j];
            }
        free_vector(pos);
        free_vector(p);
        }  
    }

    return w;
}

/* Cache frequently used parameters in ode_params to avoid repeated database lookups. */
struct ode_params initialize_ode_params() {
    struct ode_params params;
    // General parameters
    params.num_dim = get_parameter_int("num_dim");
    params.num_bodies = get_parameter_int("num_bodies");

    // Masses
    allocate_vector(&params.masses, params.num_bodies);
    for (int i = 0; i < params.num_bodies; i++)
        params.masses[i] = get_parameter_double_i("mass", i+1);

    // PN terms
    params.pn_terms = (int *)malloc(NUM_PN_TERMS * sizeof(int));
    if (params.pn_terms == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    params.pn_terms[0] = get_parameter_int("0pn_terms");
    params.pn_terms[1] = get_parameter_int("1pn_terms");
    params.pn_terms[2] = get_parameter_int("2pn_terms");
    params.pn_terms[3] = get_parameter_int("2.5pn_terms");

    return params;
}