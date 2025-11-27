#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pn_eom.h"
#include "utils.h"
#include "transformations.h"
#include "pn_eom_hamiltonians.h"


void total_angular_momentum_com(int* bodies, int num_bodies, double* w, struct ode_params* params, double* L) {
    int a, index_r, index_p;
    double r_x, r_y, r_z, p_x, p_y, p_z;
    int array_half = params->num_dim * params->num_bodies;

    // Initialize angular momentum to zero
    for (int i = 0; i < params->num_dim; i++)
        L[i] = 0.0; 

    // Transform coordinates to center of mass
    double* w_COM;
    allocate_vector(&w_COM, 2 * array_half);
    transform_to_com(bodies, num_bodies, params, w, w_COM);

    // Compute angular momentum
    for (a = 0; a < num_bodies; a++) {
        index_r = params->num_dim * bodies[a];
        index_p = array_half + params->num_dim * bodies[a];

        r_x = w_COM[index_r];
        r_y = w_COM[index_r + 1];
        r_z = w_COM[index_r + 2];

        p_x = w_COM[index_p];
        p_y = w_COM[index_p + 1];
        p_z = w_COM[index_p + 2];

        if (params->num_dim == 3) {
            L[0] += r_y * p_z - r_z * p_y;
            L[1] += r_z * p_x - r_x * p_z;
            L[2] += r_x * p_y - r_y * p_x;
        } else if (params->num_dim == 2) {
            L[0] += r_x * p_y - r_y * p_x;
        }
    }
    free_vector(w_COM);
}


double eccentricity(int* bodies, double* w, struct ode_params* params) {
    int array_half = params->num_dim * params->num_bodies;
    double L2 = 0.0; 
    // Compute squared magnitude of the total angular momentum
    if (params->num_dim == 2) {
        total_angular_momentum_com(bodies, 2, w, params, &L2);
        L2 = L2 * L2;
    }
    if (params->num_dim == 3) {
        double* L;
        allocate_vector(&L, params->num_dim);
        total_angular_momentum_com(bodies, 2, w, params, L);
        L2 = dot_product(L, L, params->num_dim);
        free_vector(L);
    }

    double* w_bodies;
    struct ode_params params_bodies;
    params_bodies.pn_terms = calloc(6, sizeof(int));
    params_bodies.num_dim = params->num_dim;
    params_bodies.num_bodies = 2;
    allocate_vector(&w_bodies, 2 * params_bodies.num_dim * 2);
    allocate_vector(&params_bodies.masses, 2);
    params_bodies.pn_terms[0] = params->pn_terms[0];
    params_bodies.pn_terms[1] = params->pn_terms[1];
    params_bodies.pn_terms[2] = params->pn_terms[2];
    params_bodies.pn_terms[3] = params->pn_terms[3];
    params_bodies.masses[0] = params->masses[bodies[0]];
    params_bodies.masses[1] = params->masses[bodies[1]];
    for (int a = 0; a < 2; a++) {
        for (int i = 0; i < params->num_dim; i++) {
            w_bodies[a * params->num_dim + i] = w[bodies[a] * params->num_dim + i];
            w_bodies[2 * params->num_dim + a * params->num_dim + i] = w[array_half + bodies[a] * params->num_dim + i];
        }
    }
    
    // Compute conservative part of the Hamiltonian
    double H_c = 0.0;
    if (params->pn_terms[0] == 1)
        H_c += H0PN(w_bodies, &params_bodies);
    if (params->pn_terms[1] == 1)
        H_c += H1PN(w_bodies, &params_bodies);
    if (params->pn_terms[2] == 1)
        H_c += H2PN_threebody(w_bodies, &params_bodies, 0);

    // Compute and return eccentricity
    double M = params->masses[bodies[0]] + params->masses[bodies[1]];
    double mu = params->masses[bodies[0]] * params->masses[bodies[1]] / M;

    free_vector(w_bodies);
    free_vector(params_bodies.masses);
    free(params_bodies.pn_terms);
    params_bodies.pn_terms = NULL;

    // Compute and return eccentricity
    return sqrt(fabs(1 + 2 * L2 * H_c / (pow(mu, 3) * pow(M, 2))));
}


void semi_major_axes(int* bodies, double* w, struct ode_params* params, double* a, double* a1_com, double* a2_com) {
    double M = params->masses[bodies[0]] + params->masses[bodies[1]];

    // Compute conservative part of the Hamiltonian
    double H_c = 0.0;
    if (params->pn_terms[0] == 1)
        H_c += H0PN(w, params);
    if (params->pn_terms[1] == 1)
        H_c += H1PN(w, params);
    if (params->pn_terms[2] == 1)
        H_c += H2PN_threebody(w, params, 0);
    
    *a = - params->masses[bodies[0]] * params->masses[bodies[1]] / (2 * H_c);
    *a1_com = params->masses[bodies[1]] / M * *a;
    *a2_com = params->masses[bodies[0]] / M * *a;
}


double total_energy_conservative(double* w, struct ode_params* params) {
    // Compute conservative part of the Hamiltonian
    double H_c = 0.0;
    if (params->pn_terms[0] == 1)
        H_c += H0PN(w, params);
    if (params->pn_terms[1] == 1)
        H_c += H1PN(w, params);
    if (params->pn_terms[2] == 1)
        H_c += H2PN_threebody(w, params, 0);
    return H_c;
}