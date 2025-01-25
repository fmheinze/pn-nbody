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
    double *p1, *p2, *p3, *x12, *x13, *x23, *n12, *n13, *n23;
    double p1_norm2, p2_norm2, p3_norm2, r12_inv, r13_inv, r23_inv, m1, m2, m3;
    double p1_p2, p1_p3, p2_p3, p1_n12, p1_n13, p1_n23, p2_n12, p2_n13, p2_n23, p3_n12, p3_n13, p3_n23;
    double r12_inv2, r13_inv2, r23_inv2, m12, m22, m32;
    int array_half; 
    array_half = 3 * params->dim;

    allocate_vector(&p1, params->dim);
    allocate_vector(&p2, params->dim);
    allocate_vector(&p3, params->dim);
    allocate_vector(&x12, params->dim);
    allocate_vector(&x13, params->dim);
    allocate_vector(&x23, params->dim);
    allocate_vector(&n12, params->dim);
    allocate_vector(&n13, params->dim);
    allocate_vector(&n23, params->dim);

    // Compute momenta
    for (int i = 0; i < params->dim; i++){
        p1[i] = params->masses[0] * w[array_half + i];
        p2[i] = params->masses[1] * w[array_half + params->dim + i];
        p3[i] = params->masses[2] * w[array_half + 2 * params->dim + i];
    }

    // Compute other frequently used quantities
    m1 = params->masses[0];
    m2 = params->masses[1];
    m3 = params->masses[2];
    for (int i = 0; i < params->dim; i++){
        x12[i] = w[i] - w[params->dim + i];
        x13[i] = w[i] - w[2 * params->dim + i];
        x23[i] = w[params->dim + i] - w[2 * params->dim + i];
    }
    p1_norm2 = dot_product(p1, p1, params->dim);
    p2_norm2 = dot_product(p2, p2, params->dim);
    p3_norm2 = dot_product(p3, p3, params->dim);
    r12_inv2 = 1/dot_product(x12, x12, params->dim);
    r13_inv2 = 1/dot_product(x13, x13, params->dim);
    r23_inv2 = 1/dot_product(x23, x23, params->dim);
    r12_inv = sqrt(r12_inv2);
    r13_inv = sqrt(r13_inv2);
    r23_inv = sqrt(r23_inv2);
    m12 = pow(m1, 2);
    m22 = pow(m2, 2);
    m32 = pow(m3, 2);
    p1_p2 = dot_product(p1, p2, params->dim);
    p1_p3 = dot_product(p1, p3, params->dim);
    p2_p3 = dot_product(p2, p3, params->dim);
    for (int i = 0; i < params->dim; i++){
        n12[i] = x12[i] * r12_inv;
        n13[i] = x13[i] * r13_inv;
        n23[i] = x23[i] * r23_inv;
    }
    p1_n12 = dot_product(p1, n12, params->dim);
    p1_n13 = dot_product(p1, n13, params->dim);
    p1_n23 = dot_product(p1, n23, params->dim);
    p2_n12 = dot_product(p2, n12, params->dim);
    p2_n13 = dot_product(p2, n13, params->dim);
    p2_n23 = dot_product(p2, n23, params->dim);
    p3_n12 = dot_product(p3, n12, params->dim);
    p3_n13 = dot_product(p3, n13, params->dim);
    p3_n23 = dot_product(p3, n23, params->dim);

    // Initialize all velocities and accelerations to zero
    for (int i = 0; i < 6 * params->dim; i++)   
        dwdt[i] = 0.0;

    // Add 0PN (Newtonian) terms
    if(params->pn_terms[0]){
        for (int i = 0; i < params->dim; i++){
            // Body 1
            dwdt[i] += w[array_half + i];
            dwdt[array_half + i] += - m2 * r12_inv2 * n12[i] - m3 * r13_inv2 * n13[i];
            // Body 2
            dwdt[params->dim + i] += w[array_half + params->dim + i];
            dwdt[array_half + params->dim + i] += m1 * r12_inv2 * n12[i] - m3 * r23_inv2 * n23[i];
            // Body 3
            dwdt[2 * params->dim + i] += w[array_half + 2 * params->dim + i];
            dwdt[array_half + 2 * params->dim + i] += m1 * r13_inv2 * n13[i] + m2 * r23_inv2 * n23[i];
        }
    }

    // Add 1PN terms
    if(params->pn_terms[1]){
        for (int i = 0; i < params->dim; i++){
            // Body 1
            dwdt[i] += - 0.5 * p1_norm2/pow(m1, 3) * p1[i] 
                       - 0.5 * r12_inv * (6 * m2/m1 * p1[i] - 7 * p2[i] - p2_n12 * n12[i])
                       - 0.5 * r13_inv * (6 * m3/m1 * p1[i] - 7 * p3[i] - p3_n13 * n13[i]);
            dwdt[array_half + i] += - n12[i] * m1*m2*r12_inv2 * (- m3*(r23_inv+r13_inv) - (m1+m2)*r12_inv 
                                             + 0.5*(3*p1_norm2/m12 + 3*p2_norm2/m22 - 7*p1_p2/(m1*m2) - 3*p1_n12*p2_n12/(m1*m2)))
                                    - n13[i] * m1*m3*r13_inv2 * (- m2*(r23_inv+r12_inv) - (m1+m3)*r13_inv 
                                             + 0.5*(3*p1_norm2/m12 + 3*p3_norm2/m32 - 7*p1_p3/(m1*m3) - 3*p1_n13*p3_n13/(m1*m3)))
                                    - 0.5*r12_inv2 * (p2_n12*p1[i] + p1_n12*p2[i]) - 0.5*r13_inv2 * (p1_n13*p3[i] + p3_n13*p1[i]);
            dwdt[array_half + i] /= m1;
            // Body 2
            dwdt[params->dim + i] += - 0.5 * p2_norm2/pow(m2, 3) * p2[i] 
                                     - 0.5 * r23_inv * (6 * m3/m2 * p2[i] - 7 * p3[i] - p3_n23 * n23[i])
                                     - 0.5 * r12_inv * (6 * m1/m2 * p2[i] - 7 * p1[i] - p1_n12 * n12[i]);
            dwdt[array_half + params->dim + i] += - n23[i] * m2*m3*r23_inv2 * (- m1*(r13_inv+r12_inv) - (m2+m3)*r23_inv 
                                                           + 0.5*(3*p2_norm2/m22 + 3*p3_norm2/m32 - 7*p2_p3/(m2*m3) - 3*p2_n23*p3_n23/(m2*m3)))
                                                  + n12[i] * m1*m2*r12_inv2 * (- m3*(r13_inv+r23_inv) - (m1+m2)*r12_inv 
                                                           + 0.5*(3*p2_norm2/m22 + 3*p1_norm2/m12 - 7*p1_p2/(m1*m2) - 3*p2_n12*p1_n12/(m1*m2)))
                                                  - 0.5*r23_inv2 * (p3_n23*p2[i] + p2_n23*p3[i]) + 0.5*r12_inv2 * (p2_n12*p1[i] + p1_n12*p2[i]);
            dwdt[array_half + params->dim + i] /= m2;                                   
            // Body 3
            dwdt[2 * params->dim + i] += - 0.5 * p3_norm2/pow(m3, 3) * p3[i] 
                                         - 0.5 * r13_inv * (6 * m1/m3 * p3[i] - 7 * p1[i] - p1_n13 * n13[i])
                                         - 0.5 * r23_inv * (6 * m2/m3 * p3[i] - 7 * p2[i] - p2_n23 * n23[i]);
            dwdt[array_half + 2 * params->dim + i] += n13[i] * m1*m3*r13_inv2 * (- m2*(r12_inv+r23_inv) - (m1+m3)*r13_inv 
                                                               + 0.5*(3*p3_norm2/m32 + 3*p1_norm2/m12 - 7*p1_p3/(m1*m3) - 3*p3_n13*p1_n13/(m1*m3)))
                                                    + n23[i] * m2*m3*r23_inv2 * (- m1*(r12_inv+r13_inv) - (m2+m3)*r23_inv 
                                                               + 0.5*(3*p3_norm2/m32 + 3*p2_norm2/m22 - 7*p2_p3/(m2*m3) - 3*p3_n23*p2_n23/(m2*m3)))
                                                    + 0.5*r13_inv2 * (p1_n13*p3[i] + p3_n13*p1[i]) + 0.5*r23_inv2 * (p3_n23*p2[i] + p2_n23*p3[i]);
            dwdt[array_half + 2 * params->dim + i] /= m3;
        }
    }

    //Add 2PN terms
    if(params->pn_terms[2]){
        update_eom_hamiltonian(w, dwdt, H2PN_threebody, 1e-5, params);
    }

    free_vector(p1);
    free_vector(p2);
    free_vector(p3);
    free_vector(x12);
    free_vector(x13);
    free_vector(x23);
    free_vector(n12);
    free_vector(n13);
    free_vector(n23);
}