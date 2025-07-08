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
    params.num_bodies = 3;

    allocate_vector(&w0, 2 * params.dim * params.num_bodies);
    allocate_vector(&params.masses, params.num_bodies);

    // Simulation parameters
    t_end = 1000;
    dt0 = 0.1;
    dt_save = 1;
    rel_error = 1e-6;
    params.pn_terms[0] = 1;
    params.pn_terms[1] = 1;
    params.pn_terms[2] = 1;
    params.pn_terms[3] = 1;

    // Body masses and initial positions/velocities
    params.masses[0] = 0.5;
    params.masses[1] = 0.5;
    params.masses[2] = 0.5;
    
    double a1 = 60.0;
    double e1 = 0.0;
    double phi01 = 5.2;
    double d0 = 200;
    double v_rel = 0.05;
    double b = 150;
    double orientation[3] = {0, 0, 1};

    // Compute initial positions and velocities
    //binary_single_scattering_symmetric(&params, w0, a1, e1, phi01, d0, v_rel, b, NULL);
    //newtonian_binary(&params, w0, a1, e1, phi01);

    //figure_eight_orbit(&params, w0, 108.3);

    d0 = 150.0;
    //double scatter_p0 = 0.11456439/2;
    double scatter_p0 = 0.05;
    b = 52.0;
    double phi0 = 1.7;
    double binary_r0 = 5.9859;
    double binary_pt0 = 0.085168;
    double binary_pr0 = 0.000388;

    //binary_single_scattering(&params, w0, d0, scatter_p0, b, phi0, binary_r0, binary_pt0, binary_pr0, orientation, 1);
    //binary_binary_scattering(&params, w0, d0, scatter_p0, b, 0.0, 4.2, binary_r0, binary_pt0, binary_pr0,
    //                          binary_r0, binary_pt0, binary_pr0, orientation, orientation, 0, 0);

    w0[0] = -47.007050;
    w0[1] = 5.183941;
    w0[2] = 0.0;
    w0[3] = -52.992950;
    w0[4] = -5.183941;
    w0[5] = 0.0;
    w0[6] = 50.0;
    w0[7] = 0.0;
    w0[8] = 0.0;
    w0[9] = -0.020966;
    w0[10] = 0.020481;
    w0[11] = 0.0;
    w0[12] = 0.126937;
    w0[13] = -0.064015;
    w0[14] = 0.0;
    w0[15] = -0.052985;
    w0[16] = 0.021767;
    w0[17] = 0.0;



    //binary_single_scattering_symmetric(&params, w0, 150, e1, phi01, 10000, 0.05, 300, orientation);

    printf("Starting simulation...\n"); 

    for (int i = 0; i < 2 * params.dim * params.num_bodies; i++)
        printf("w0[%d] = %lf\n", i, w0[i]);

    // Run simulation
    ode_integrator(w0, t_end, dt0, dt_save, rel_error, rhs_pn_threebody, &params);

    free_vector(w0);
    free_vector(params.masses);
    free(params.pn_terms);
    params.pn_terms = NULL;
    return 0;
}