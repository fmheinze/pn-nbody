#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pn_eom.h"
#include "utils.h"


void transform_to_com(int* bodies, int num_bodies, struct ode_params* params, double* w, double* w_transform) {
    int array_half = params->dim * params->num_bodies; 
    double M = 0.0;
    double* R_COM, *V_COM;
    
    // Allocate and initialize position and velocity of the COM
    allocate_vector(&R_COM, params->dim);
    allocate_vector(&V_COM, params->dim);
    for (int i = 0; i < params->dim; i++) {
        R_COM[i] = 0.0;
        V_COM[i] = 0.0;
    }

    // Compute total mass
    for (int a = 0; a < num_bodies; a++)
        M += params->masses[bodies[a]];
    
    // Compute position and velocity of the COM
    for (int a = 0; a < num_bodies; a++) {
        for (int i = 0; i < params->dim; i++) {
            R_COM[i] += params->masses[bodies[a]] * w[params->dim * bodies[a] + i];
            V_COM[i] += w[array_half + params->dim * bodies[a] + i];
        }
    }

    // Normalize COM position and velocity
    for (int i = 0; i < params->dim; i++) {
        R_COM[i] /= M;
        V_COM[i] /= M;
    }

    // Transform all positions and velocities to the COM frame
    for (int a = 0; a < params->num_bodies; a++) {
        for (int i = 0; i < params->dim; i++) {
            w_transform[params->dim * a + i] = w[params->dim * a + i] - R_COM[i];
            w_transform[array_half + params->dim * a + i] = w[array_half + params->dim * a + i] - V_COM[i];
        }
    }

    free_vector(R_COM);
    free_vector(V_COM);
}