#include <stdlib.h>
#include <stdio.h>
#include "integration.h"
#include "utils.h"
#include "pn_eom.h"
#include "initial_configurations.h"


int main(){
    double* w0;
    size_t w_size;
    double t_end, dt0, dt_save, rel_error;
    struct ode_params params;
    params.pn_terms = calloc(6, sizeof(int));

    // Set number of dimensions and number of bodies
    params.dim = 3;
    params.num_bodies = 3;

    w_size = 2 * params.dim * params.num_bodies;
    allocate_vector(&w0, w_size);
    allocate_vector(&params.masses, params.num_bodies);

    // Simulation parameters
    t_end = 1e5;
    dt0 = 0.1;
    dt_save = 10;
    rel_error = 1e-8;
    params.pn_terms[0] = 1;
    params.pn_terms[1] = 1;
    params.pn_terms[2] = 1;
    params.pn_terms[3] = 0;

    // Body masses and initial positions/velocities
    params.masses[0] = 1.0;
    params.masses[1] = 1.0;
    params.masses[2] = 1.0;
    
    //double a1 = 100.0;
    //double e1 = 0.0;
    //double phi01 = 1.2;
    //double d0 = 3000;
    //double v_rel = 0.1;
    //double b = 400;
    //double orientation[3] = {0, 0, 1};

    // Compute initial positions and velocities
    // binary_single_scattering_symmetric(params.masses[0], params.masses[1], params.masses[2], 
    //                                   a1, e1, phi01, d0, v_rel, b, NULL, w0);

    figure_eight_orbit(108.1, &params, w0);

    printf("Starting simulation...\n");

    // Run simulation
    ode_integrator(w0, w_size, 0.0, t_end, dt0, dt_save, rel_error, rhs_pn_threebody, &params);

    free_vector(w0);
    free_vector(params.masses);
    free(params.pn_terms);
    params.pn_terms = NULL;
    return 0;
}