#include <stdlib.h>
#include <stdio.h>
#include "integration.h"
#include "utils.h"
#include "pn_eom.h"
#include "initial_configurations.h"
#include "transformations.h"
#include "orbital_elements.h"


int main(){
    double* w0;
    double t_end, dt0, dt_save, rel_error;
    struct ode_params params;
    params.pn_terms = calloc(6, sizeof(int));

    // Set number of dimensions and number of bodies
    params.dim = 3;
    params.num_bodies = 2;

    allocate_vector(&w0, 2 * params.dim * params.num_bodies);
    allocate_vector(&params.masses, params.num_bodies);

    // Simulation parameters
    t_end =2.9e7;
    dt0 = 0.1;
    dt_save = 1000;
    rel_error = 1e-12;
    params.pn_terms[0] = 1;
    params.pn_terms[1] = 1;
    params.pn_terms[2] = 0;
    params.pn_terms[3] = 1;

    // Body masses and initial positions/velocities
    params.masses[0] = 1.0;
    params.masses[1] = 1.0;
    params.masses[2] = 0.001;
    
    double a1 = 300.0;
    double e1 = 0.5;
    double phi01 = 3.1;
    double d0 = 10000;
    double v_rel = 0.001;
    double b = 3000;
    double orientation[3] = {0, 0, 1};

    // Compute initial positions and velocities
    //binary_single_scattering_symmetric(&params, w0, a1, e1, phi01, d0, v_rel, b, orientation);
    newtonian_binary(&params, w0, a1, e1, phi01);

    //figure_eight_orbit(&params, w0, 108.3);

    printf("Starting simulation...\n"); 

    // Run simulation
    ode_integrator(w0, t_end, dt0, dt_save, rel_error, rhs_pn_threebody, &params);

    free_vector(w0);
    free_vector(params.masses);
    free(params.pn_terms);
    params.pn_terms = NULL;
    return 0;
}