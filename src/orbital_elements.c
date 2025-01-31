#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pn_eom.h"
#include "utils.h"
#include "transformations.h"


void total_angular_momentum_com(int* bodies, int num_bodies, double* w, struct ode_params* params, double* L) {
    int array_half = params->dim * params->num_bodies; 
    double* w_COM, *r_COM, *p_COM, *L_body;

    allocate_vector(&r_COM, params->dim);
    allocate_vector(&p_COM, params->dim);
    allocate_vector(&L_body, params->dim);
    allocate_vector(&w_COM, 2 * array_half);
    for (int i = 0; i < params->dim; i++)
        L[i] = 0.0;

    transform_to_com(bodies, num_bodies, params, w, w_COM);

    for (int a = 0; a < num_bodies; a++) {
        for (int i = 0; i < params->dim; i++) {
            r_COM[i] = w_COM[params->dim * bodies[a] + i];
            p_COM[i] = params->masses[bodies[a]] * w_COM[array_half + params->dim * bodies[a] + i];
        }
        cross_product(r_COM, p_COM, L_body);
        for (int i = 0; i < params->dim; i++)
            L[i] += L_body[i];
    }

    free_vector(r_COM);
    free_vector(p_COM);
    free_vector(w_COM);
    free_vector(L_body);
}


double eccentricity(int body1, int body2, struct ode_params* params) {

}