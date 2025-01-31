#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "pn_eom.h"
#include "pn_eom_hamiltonians.h"


void rhs_pn_twobody(double t, double* w, struct ode_params* params, double* dwdt)
/* Right-hand side of the first order ODE system describing the evolution of two bodies under the influence of gravity,
including Post-Newtonian correction terms up to 2.5PN order. More than two bodies can be added but nonlinear N-body 
effects will not be accounted for.

t           time
w           pointer to the state vector w = [x1, y1, z1, ..., xn, yn, zn, vx1, vy1, vz1, ..., vxn, vyn, vzn]
params      contains the number of dimensions, the number of bodies, the bodies masses 
            as well the PN terms that should be included in the computation of the accelerations of the bodies
dwdt        pointer to the array of values that will be updated with the right-hand side of the ODE*/
{
    double  a_temp, dist2, 
            mimj, mi2, mj2, mi3, mj3,
            v_rel2, vi2, vj2, vi4, vj4, vi_dot_vj, vi_dot_vj2,
            pos_dot_vi, pos_dot_vj, pos_dot_vi2, pos_dot_vj2, pos_dot_vi3, pos_dot_vj3, pos_dot_vi4, pos_dot_vj4, pos_dot_v_rel;
    double* pos_rel, *pos_dist3, *pos_dist4, *pos_dist5, *v_rel, *v_dist2, *v_dist3, *v_dist4;
    int array_half, vi_coord_index, vj_coord_index; 

    allocate_vector(&pos_rel, params->dim);
    allocate_vector(&pos_dist3, params->dim);
    allocate_vector(&pos_dist4, params->dim);
    allocate_vector(&pos_dist5, params->dim);
    allocate_vector(&v_rel, params->dim);
    allocate_vector(&v_dist2, params->dim);
    allocate_vector(&v_dist3, params->dim);
    allocate_vector(&v_dist4, params->dim);
    array_half = params->dim * params->num_bodies;
    
    for (int i = 0; i < array_half; i++){
        // Initialize accelerations to zero
        dwdt[array_half + i] = 0.0;

        // Set velocities dx/dt = v 
        dwdt[i] = w[array_half + i];
    }

    // Compute accelerations dv/dt = a (looping over all bodies i, j with i < j)
    for (int i = 0; i < params->num_bodies; i++){
        for (int j = i + 1; j < params->num_bodies; j++){

            // Compute relative position of the bodies i and j and their distance
            dist2 = 0;
            for (int coord = 0; coord < params->dim; coord++){
                pos_rel[coord] = w[params->dim * i + coord] - w[params->dim * j + coord];
                dist2 += pow(pos_rel[coord], 2);
            }
            for(int coord = 0; coord < params->dim; coord++)
                pos_dist3[coord] = pos_rel[coord] / pow(dist2, 1.5);


            // Compute the accelerations of the bodies i and j using Newton's 3rd law

            // Add 0PN (Newtonian) terms
            if(params->pn_terms[0]){
                for (int coord = 0; coord < params->dim; coord++){
                    dwdt[array_half + params->dim * i + coord] -= params->masses[j] * pos_dist3[coord];
                    dwdt[array_half + params->dim * j + coord] += params->masses[i] * pos_dist3[coord];
                }
            }

            // Frequently used quantities in the >=1PN terms
            if(params->pn_terms[1] | params->pn_terms[2] | params->pn_terms[3] | params->pn_terms[4] | params->pn_terms[5]){
                vi2 = vj2 = vi_dot_vj = pos_dot_vi = pos_dot_vj = 0;
                for(int coord = 0; coord < params->dim; coord++){
                    vi_coord_index = array_half + params->dim * i + coord;
                    vj_coord_index = array_half + params->dim * j + coord;
                    v_rel[coord] = w[vi_coord_index] - w[vj_coord_index];
                    v_dist2[coord] = v_rel[coord] / dist2;
                    pos_dist4[coord] = pos_rel[coord] / pow(dist2, 2);
                    pos_dot_vi += pos_rel[coord] / dist2 * w[vi_coord_index];
                    pos_dot_vj += pos_rel[coord] / dist2 * w[vj_coord_index];
                    vi_dot_vj += w[vi_coord_index] * w[vj_coord_index];
                    vi2 += pow(w[vi_coord_index], 2);
                    vj2 += pow(w[vj_coord_index], 2);
                }
                pos_dot_vi2 = pow(pos_dot_vi, 2);
                pos_dot_vj2 = pow(pos_dot_vj, 2);
                mimj = params->masses[i] * params->masses[j];
                mi2 = pow(params->masses[i], 2);
                mj2 = pow(params->masses[j], 2);
            }

            // Add 1PN terms
            if(params->pn_terms[1]){
                for (int coord = 0; coord < params->dim; coord++){
                    vi_coord_index = array_half + params->dim * i + coord;
                    vj_coord_index = array_half + params->dim * j + coord;

                    a_temp = 5.0 * mimj * pos_dist4[coord];
                    dwdt[vi_coord_index] += a_temp;
                    dwdt[vj_coord_index] -= a_temp;

                    a_temp = 4.0 * pos_dist4[coord];
                    dwdt[vi_coord_index] += mj2 * a_temp;
                    dwdt[vj_coord_index] -= mi2 * a_temp;

                    a_temp = 1.5 * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vj2 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vi2 * a_temp;

                    dwdt[vi_coord_index] -= params->masses[j] * vi2 * pos_dist3[coord];
                    dwdt[vj_coord_index] += params->masses[i] * vj2 * pos_dist3[coord];

                    a_temp = 4.0 * pos_dist3[coord] * vi_dot_vj;
                    dwdt[vi_coord_index] += params->masses[j] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * a_temp;

                    a_temp = -2.0 * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * vj2 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * vi2 * a_temp;

                    a_temp = 4.0 * v_dist2[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vi * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vj * a_temp;

                    a_temp = -3.0 * v_dist2[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vj * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vi * a_temp;
                }
            }

            // Frequently used quantities in the >=2PN terms
            if(params->pn_terms[2] | params->pn_terms[3] | params->pn_terms[4] | params->pn_terms[5]){
                for(int coord = 0; coord < params->dim; coord++){
                    v_dist3[coord] = v_rel[coord] / pow(dist2, 1.5);
                    pos_dist5[coord] = pos_rel[coord] / pow(dist2, 2.5);
                }
                pos_dot_vi3 = pow(pos_dot_vi, 3);
                pos_dot_vj3 = pow(pos_dot_vj, 3);
                pos_dot_vi4 = pow(pos_dot_vi, 4);
                pos_dot_vj4 = pow(pos_dot_vj, 4);
                vi_dot_vj2 = pow(vi_dot_vj, 2);
                mi3 = pow(params->masses[i], 3);
                mj3 = pow(params->masses[j], 3);
                vi4 = pow(vi2, 2);
                vj4 = pow(vj2, 2);
            }

            // Add 2PN terms
            if(params->pn_terms[2]){
                for (int coord = 0; coord < params->dim; coord++){
                    vi_coord_index = array_half + params->dim * i + coord;
                    vj_coord_index = array_half + params->dim * j + coord;

                    a_temp = -14.25 * mimj * pos_dist5[coord];
                    dwdt[vi_coord_index] += params->masses[i] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[j] * a_temp;

                    a_temp = -34.5 * mimj * pos_dist5[coord];
                    dwdt[vi_coord_index] += params->masses[j] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * a_temp;

                    a_temp = -9.0 * pos_dist5[coord];
                    dwdt[vi_coord_index] += mj3 * a_temp;
                    dwdt[vj_coord_index] -= mi3 * a_temp;

                    a_temp = -1.875 * pos_dot_vj4 * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * a_temp;

                    a_temp = 1.5 * pos_dot_vj2 * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * vi2 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * vj2 * a_temp;

                    a_temp = -6.0 * pos_dot_vj2 * pos_dist3[coord] * vi_dot_vj;
                    dwdt[vi_coord_index] += params->masses[j] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * a_temp;

                    a_temp = -2.0 * vi_dot_vj2 * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * a_temp;

                    a_temp = 4.5 * pos_dot_vj2 * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * vj2 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * vi2 * a_temp;

                    a_temp = 4.0 * vi_dot_vj * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * vj2 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * vi2 * a_temp;

                    a_temp = -2.0 * pos_dist3[coord];
                    dwdt[vi_coord_index] += params->masses[j] * vj4 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * vi4 * a_temp;

                    a_temp = 19.5 * mimj * pos_dist4[coord];
                    dwdt[vi_coord_index] += pos_dot_vi2 * a_temp;
                    dwdt[vj_coord_index] -= pos_dot_vj2 * a_temp;

                    a_temp = -39.0 * mimj * pos_dot_vi * pos_dot_vj * pos_dist4[coord];
                    dwdt[vi_coord_index] += a_temp;
                    dwdt[vj_coord_index] -= a_temp;

                    a_temp = 8.5 * mimj * pos_dist4[coord];
                    dwdt[vi_coord_index] += pos_dot_vj2 * a_temp;
                    dwdt[vj_coord_index] -= pos_dot_vi2 * a_temp;

                    a_temp = -3.75 * mimj * pos_dist4[coord];
                    dwdt[vi_coord_index] += vi2 * a_temp;
                    dwdt[vj_coord_index] -= vj2 * a_temp;

                    a_temp = -2.5 * mimj * vi_dot_vj * pos_dist4[coord];
                    dwdt[vi_coord_index] += a_temp;
                    dwdt[vj_coord_index] -= a_temp;

                    a_temp = 2.0 * pos_dist4[coord];
                    dwdt[vi_coord_index] += mj2 * pos_dot_vi2 * a_temp;
                    dwdt[vj_coord_index] -= mi2 * pos_dot_vj2 * a_temp;

                    a_temp = -4.0 * pos_dot_vi * pos_dot_vj * pos_dist4[coord];
                    dwdt[vi_coord_index] += mj2 * a_temp;
                    dwdt[vj_coord_index] -= mi2 * a_temp;

                    a_temp = -6.0 * pos_dist4[coord];
                    dwdt[vi_coord_index] += mj2 * pos_dot_vj2 * a_temp;
                    dwdt[vj_coord_index] -= mi2 * pos_dot_vi2 * a_temp;

                    a_temp = -8.0 * vi_dot_vj * pos_dist4[coord];
                    dwdt[vi_coord_index] += mj2 * a_temp;
                    dwdt[vj_coord_index] -= mi2 * a_temp;

                    a_temp = 4.0 * pos_dist4[coord];
                    dwdt[vi_coord_index] += mj2 * vj2 * a_temp;
                    dwdt[vj_coord_index] -= mi2 * vi2 * a_temp;

                    a_temp = -2.0 * v_dist3[coord];
                    dwdt[vi_coord_index] += mj2 * pos_dot_vi * a_temp;
                    dwdt[vj_coord_index] -= mi2 * pos_dot_vj * a_temp;

                    a_temp = -2.0 * v_dist3[coord];
                    dwdt[vi_coord_index] += mj2 * pos_dot_vj * a_temp;
                    dwdt[vj_coord_index] -= mi2 * pos_dot_vi * a_temp;

                    a_temp = -15.75 * mimj * v_dist3[coord];
                    dwdt[vi_coord_index] += pos_dot_vi * a_temp;
                    dwdt[vj_coord_index] -= pos_dot_vj * a_temp;

                    a_temp = 13.75 * mimj * v_dist3[coord];
                    dwdt[vi_coord_index] += pos_dot_vj * a_temp;
                    dwdt[vj_coord_index] -= pos_dot_vi * a_temp;

                    a_temp = -6.0 * pos_dot_vi * pos_dot_vj * v_dist2[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vj * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vi * a_temp;

                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vj * vi2 * v_dist2[coord];
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vi * vj2 * v_dist2[coord];

                    a_temp = -4.0 * vi_dot_vj * v_dist2[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vi * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vj * a_temp;

                    a_temp = 4.0 * vi_dot_vj * v_dist2[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vj * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vi * a_temp;

                    a_temp = 4.0 * v_dist2[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vi * vj2 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vj * vi2 * a_temp;

                    a_temp = -5.0 * v_dist2[coord];
                    dwdt[vi_coord_index] += params->masses[j] * pos_dot_vj * vj2 * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * pos_dot_vi * vi2 * a_temp;
                }
            }

            // Frequently used quantities in the >=2.5PN terms
            if(params->pn_terms[3] | params->pn_terms[4] | params->pn_terms[5]){
                pos_dot_v_rel = v_rel2 = 0;
                for(int coord = 0; coord < params->dim; coord++){
                    pos_dot_v_rel += pos_rel[coord] / dist2 * v_rel[coord];
                    v_rel2 += pow(v_rel[coord], 2);
                }
            }

            // Add 2.5PN terms
            if(params->pn_terms[3]){
                for (int coord = 0; coord < params->dim; coord++){
                    vi_coord_index = array_half + params->dim * i + coord;
                    vj_coord_index = array_half + params->dim * j + coord;

                    a_temp = 208.0/15.0 * mimj * pos_dot_v_rel * pos_dist5[coord];
                    dwdt[vi_coord_index] += params->masses[j] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * a_temp;

                    a_temp = -4.8 * mimj * pos_dot_v_rel * pos_dist5[coord];
                    dwdt[vi_coord_index] += params->masses[i] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[j] * a_temp;

                    a_temp = 2.4 * mimj * pos_dot_v_rel * v_rel2 * pos_dist4[coord];
                    dwdt[vi_coord_index] += a_temp;
                    dwdt[vj_coord_index] -= a_temp;

                    a_temp = 1.6 * mimj * v_dist4[coord];
                    dwdt[vi_coord_index] += params->masses[i] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[j] * a_temp;

                    a_temp = -6.4 * mimj * v_dist4[coord];
                    dwdt[vi_coord_index] += params->masses[j] * a_temp;
                    dwdt[vj_coord_index] -= params->masses[i] * a_temp;

                    a_temp = -0.8 * mimj * v_rel2 * v_dist3[coord];
                    dwdt[vi_coord_index] += a_temp;
                    dwdt[vj_coord_index] -= a_temp;
                }
            }
        }
    }
    free_vector(pos_rel);
    free_vector(pos_dist3);
    free_vector(pos_dist4);
    free_vector(pos_dist5);
    free_vector(v_rel);
    free_vector(v_dist2);
    free_vector(v_dist3);
    free_vector(v_dist4);
}


void rhs_pn_threebody(double t, double* w, struct ode_params* params, double* dwdt)
/* Right-hand side of the first order ODE system describing the evolution of three bodies under the influence of gravity,
including Post-Newtonian correction terms up to 2.5PN order.

t           time
w           pointer to the state vector w = [x1, y1, z1, ..., x3, y3, z3, vx1, vy1, vz1, ..., vx3, vy3, vz3]
params      contains the number of dimensions, the number of bodies, the bodies masses 
            as well the PN terms that should be included in the computation of the accelerations of the bodies
dwdt        pointer to the array of values that will be updated with the right-hand side of the ODE*/
{
    int array_half = 3 * params->dim; 

    // Masses
    double* m;
    allocate_vector(&m, params->num_bodies);
    for (int a = 0; a < params->num_bodies; a++)
        m[a] = params->masses[a];
    
    // Momenta
    double** p;
    allocate_array(&p, params->num_bodies, params->dim);
    for (int a = 0; a < params->num_bodies; a++)
        for (int i = 0; i < params->dim; i++)
            p[a][i] = m[a] * w[array_half + a * params->dim + i];
    
    // Relative positions and distances
    double*** x_rel, **r, ***n;
    allocate_3d_array(&x_rel, params->num_bodies, params->num_bodies, params->dim);
    allocate_3d_array(&n, params->num_bodies, params->num_bodies, params->dim);
    allocate_array(&r, params->num_bodies, params->num_bodies);
    for (int a = 0; a < params->num_bodies; a++) {
        for (int b = a; b < params->num_bodies; b++) {
            for (int i = 0; i < params->dim; i++){
                x_rel[a][b][i] = w[a * params->dim + i] - w[b * params->dim + i];
                x_rel[b][a][i] = -x_rel[a][b][i];
            } 
            r[a][b] = norm(x_rel[a][b], params->dim);
            r[b][a] = r[a][b];
            for (int i = 0; i < params->dim; i++){
                if (a == b){
                    n[a][b][i] = 0.0;
                    n[b][a][i] = 0.0;
                } else {
                    n[a][b][i] = x_rel[a][b][i] / r[a][b];
                    n[b][a][i] = -n[a][b][i];
                }
            }
        }
    }

    // Time derivatives
    double**p_dot, **x_dot;
    allocate_array(&p_dot, params->num_bodies, params->dim);
    allocate_array(&x_dot, params->num_bodies, params->dim);
    for (int a = 0; a < params->num_bodies; a++) {
        for (int i = 0; i < params->dim; i++) {
            p_dot[a][i] = 0.0;
            x_dot[a][i] = 0.0;
        }
    }
            
    // Initialize all velocities and accelerations to zero
    for (int i = 0; i < 6 * params->dim; i++)   
        dwdt[i] = 0.0;

    // Add 0PN (Newtonian) terms
    if (params->pn_terms[0]) {
        // Velocities
        for (int a = 0; a < params->num_bodies; a++)
            for (int i = 0; i < params->dim; i++) {
                dwdt[a * params->dim + i] += w[array_half + a * params->dim + i];
                x_dot[a][i] += w[array_half + a * params->dim + i];
            }
    
        //Accelerations
        for (int a = 0; a < params->num_bodies; a++) {
            for (int b = a+1; b < params->num_bodies; b++) {
                for (int i = 0; i < params->dim; i++) {
                    dwdt[array_half + a * params->dim + i] += -m[b] * 1/pow(r[a][b], 2) * n[a][b][i];
                    dwdt[array_half + b * params->dim + i] += -m[a] * 1/pow(r[a][b], 2) * n[b][a][i];
                    p_dot[a][i] += -m[b] * 1/pow(r[a][b], 2) * n[a][b][i];
                    p_dot[b][i] += -m[a] * 1/pow(r[a][b], 2) * n[b][a][i];
                }
            }
        }
    }

    // Add 1PN terms
    if (params->pn_terms[1]) {
        // Velocities
        for (int a = 0; a < params->num_bodies; a++) {
            for (int i = 0; i < params->dim; i++) {
                x_dot[a][i] += -0.5 * dot_product(p[a], p[a], params->dim) / pow(m[a], 3) * p[a][i];
                dwdt[a * params->dim + i] += -0.5 * dot_product(p[a], p[a], params->dim) / pow(m[a], 3) * p[a][i];
                for (int b = 0; b < params->num_bodies; b++) {
                    if (b != a) {
                        x_dot[a][i] += -0.5 * 1/r[a][b] * (6 * m[b]/m[a] * p[a][i] - 7 * p[b][i] - dot_product(n[a][b], p[b], params->dim) * n[a][b][i]);
                        dwdt[a * params->dim + i] += -0.5 * 1/r[a][b] * (6 * m[b]/m[a] * p[a][i] - 7 * p[b][i] - dot_product(n[a][b], p[b], params->dim) * n[a][b][i]);
                    }
                }
            }
        }

        // Accelerations
        for (int a = 0; a < params->num_bodies; a++) {
            for (int i = 0; i < params->dim; i++) {
                for (int b = 0; b < params->num_bodies; b++) {
                    if (b != a) {
                        p_dot[a][i] += -0.5 / pow(r[a][b], 2) * (3 * m[b]/m[a] * dot_product(p[a], p[a], params->dim) + 3 * m[a]/m[b] * dot_product(p[b], p[b], params->dim) 
                                               - 7 * dot_product(p[a], p[b], params->dim) - 3 * dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][b], p[b], params->dim)) * n[a][b][i];
                        dwdt[array_half + a * params->dim + i] += -0.5 / pow(r[a][b], 2) * (3 * m[b]/m[a] * dot_product(p[a], p[a], params->dim) + 3 * m[a]/m[b] * dot_product(p[b], p[b], params->dim) 
                                               - 7 * dot_product(p[a], p[b], params->dim) - 3 * dot_product(n[a][b], p[a], params->dim) * dot_product(n[a][b], p[b], params->dim)) * n[a][b][i];
                        p_dot[a][i] += -0.5 / pow(r[a][b], 2) * (dot_product(n[a][b], p[b], params->dim) * p[a][i] + dot_product(n[a][b], p[a], params->dim) * p[b][i]);
                        dwdt[array_half + a * params->dim + i] += -0.5 / pow(r[a][b], 2) * (dot_product(n[a][b], p[b], params->dim) * p[a][i] + dot_product(n[a][b], p[a], params->dim) * p[b][i]);
                    }
                    for (int c = 0; c < params->num_bodies; c++) {
                        if (b != a && c != a) {
                            p_dot[a][i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[a][c]) * n[a][b][i];
                            dwdt[array_half + a * params->dim + i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[a][c]) * n[a][b][i];
                        }
                        if (b != a && c != b) {
                            p_dot[a][i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[b][c]) * n[a][b][i];
                            dwdt[array_half + a * params->dim + i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[b][c]) * n[a][b][i];
                        }
                    }
                }
            }
        }
    }

    // Add 2PN terms (using finite differencing on the hamiltonian)
    if (params->pn_terms[2]) {
        update_eom_hamiltonian(w, dwdt, H2PN_threebody, 1e-4, params);
    }

    // Add 2.5PN terms
    if (params->pn_terms[3]) {
        double** chi_dot, ****dp_chi, ****dx_chi;
        double*** x_rel_dot;

        allocate_array(&chi_dot, params->dim, params->dim);
        allocate_3d_array(&x_rel_dot, params->num_bodies, params->num_bodies, params->dim);
        allocate_4d_array(&dp_chi, params->num_bodies, params->dim, params->dim, params->dim);
        allocate_4d_array(&dx_chi, params->num_bodies, params->dim, params->dim, params->dim);

        for (int a = 0; a < params->num_bodies; a++) {
            for (int b = 0; b < params->num_bodies; b++) {
                for (int i = 0; i < params->dim; i++) {
                    x_rel_dot[a][b][i] = x_dot[a][i] - x_dot[b][i];
                }
            }
        }

        // Initialize chi arrays to zero
        for (int i = 0; i < params->dim; i++) {
            for (int j = 0; j < params->dim; j++) {
                chi_dot[i][j] = 0.0;
            }
        }
        for (int a = 0; a < params->dim; a++) {
            for (int i = 0; i < params->dim; i++) {
                for (int j = 0; j < params->dim; j++) {
                    for (int k = 0; k < params->dim; k++) {
                        dp_chi[a][i][j][k] = 0.0;
                        dx_chi[a][i][j][k] = 0.0;
                    }
                }
            }
        }    
        // Compute chi values
        for (int i = 0; i < params->dim; i++) {
            for (int j = 0; j < params->dim; j++) {
                for (int a = 0; a < params->num_bodies; a++) {
                    chi_dot[i][j] += 2.0 / m[a] * (2 * dot_product(p_dot[a], p[a], params->dim) * kronecker_delta(i, j) - 3 * (p_dot[a][i] * p[a][j] + p[a][i] * p_dot[a][j]));
                }
                for (int a = 0; a < params->num_bodies; a++) {
                    for (int b = 0; b < params->num_bodies; b++) {
                        if (b != a) {
                            chi_dot[i][j] += m[a] * m[b] / pow(r[a][b], 2) * (3 * (x_rel_dot[a][b][i] * n[a][b][j] + n[a][b][i] * x_rel_dot[a][b][j]) +
                                                                              dot_product(n[a][b], x_rel_dot[a][b], params->dim) * (kronecker_delta(i, j) - 9 * n[a][b][i] * n[a][b][j]));
                        }
                    }
                }
            }
        }

        for (int c = 0; c < params->num_bodies; c++) {
            for (int i = 0; i < params->dim; i++) {
                for (int j = 0; j < params->dim; j++) {
                    for (int k = 0; k < params->dim; k++) {
                        dp_chi[c][i][j][k] += 2.0 / m[c] * (2 * p[c][k] * kronecker_delta(i, j) - 3 * (p[c][j] * kronecker_delta(i, k) + p[c][i] * kronecker_delta(j, k)));
                    }
                }
            }
        }
        for (int c = 0; c < params->num_bodies; c++) {
            for (int i = 0; i < params->dim; i++) {
                for (int j = 0; j < params->dim; j++) {
                    for (int k = 0; k < params->dim; k++) {
                        for (int a = 0; a < params->num_bodies; a++) {
                            for (int b = 0; b < params->num_bodies; b++) {
                                if (b != a) {
                                    dx_chi[c][i][j][k] += m[a] * m[b] / pow(r[a][b], 2) * (kronecker_delta(a, c) - kronecker_delta(b, c)) * 
                                    (3 * (kronecker_delta(i, k) * n[a][b][j] + kronecker_delta(j, k) * n[a][b][i]) - 9 * n[a][b][k] * n[a][b][i] * n[a][b][j] + kronecker_delta(i, j) * n[a][b][k]);
                                }
                            }
                        }
                    }    
                }
            }
        }

        for (int a = 0; a < params->num_bodies; a++) {
            for (int k = 0; k < params->dim; k++) {
                for (int i = 0; i < params->dim; i++) {
                    for (int j = 0; j < params->dim; j++) {
                        dwdt[a * params->dim + k] += 1/45.0 * chi_dot[i][j] * dp_chi[a][i][j][k];
                        dwdt[array_half + a * params->dim + k] += -1/45.0 * chi_dot[i][j] * dx_chi[a][i][j][k];
                    }
                }
            }
        }

        free_array(chi_dot, params->dim);
        free_3d_array(x_rel_dot, params->num_bodies, params->num_bodies);
        free_4d_array(dp_chi, params->num_bodies, params->dim, params->dim);
        free_4d_array(dx_chi, params->num_bodies, params->dim, params->dim);
    }

    free_vector(m);
    free_array(p, params->num_bodies);
    free_array(r, params->num_bodies);
    free_array(p_dot, params->num_bodies);
    free_array(x_dot, params->num_bodies);
    free_3d_array(x_rel, params->num_bodies, params->num_bodies);
    free_3d_array(n, params->num_bodies, params->num_bodies);
}