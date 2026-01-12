#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "utils.h"
#include "pn_eom.h"
#include "pn_eom_hamiltonians.h"


void rhs_pn_threebody(double t, double* w, struct ode_params* params, double* dwdt)
/* Right-hand side of the first order ODE system describing the evolution of three bodies under the influence of gravity,
including Post-Newtonian correction terms up to 2.5PN order.

t           time
w           pointer to the state vector w = [x1, y1, z1, ..., x3, y3, z3, vx1, vy1, vz1, ..., vx3, vy3, vz3]
params      contains the number of dimensions, the number of bodies, the bodies masses 
            as well the PN terms that should be included in the computation of the accelerations of the bodies
dwdt        pointer to the array of values that will be updated with the right-hand side of the ODE*/
{
    int num_bodies = params->num_bodies;
    int num_dim = params->num_dim;
    int array_half = num_bodies * num_dim; 

    // Masses
    double* m;
    allocate_vector(&m, num_bodies);
    for (int a = 0; a < num_bodies; a++)
        m[a] = params->masses[a];
    
    // Momenta
    double** p;
    allocate_2d_array(&p, num_bodies, num_dim);
    for (int a = 0; a < num_bodies; a++)
        for (int i = 0; i < num_dim; i++)
            p[a][i] = w[array_half + a * num_dim + i];
    
    // Relative positions and distances
    double*** x_rel, **r, ***n;
    allocate_3d_array(&x_rel, num_bodies, num_bodies, num_dim);
    allocate_3d_array(&n, num_bodies, num_bodies, num_dim);
    allocate_2d_array(&r, num_bodies, num_bodies);
    for (int a = 0; a < num_bodies; a++) {
        for (int b = a; b < num_bodies; b++) {
            for (int i = 0; i < num_dim; i++){
                x_rel[a][b][i] = w[a * num_dim + i] - w[b * num_dim + i];
                x_rel[b][a][i] = -x_rel[a][b][i];
            } 
            r[a][b] = norm(x_rel[a][b], num_dim);
            r[b][a] = r[a][b];
            for (int i = 0; i < num_dim; i++){
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
    allocate_2d_array(&p_dot, num_bodies, num_dim);
    allocate_2d_array(&x_dot, num_bodies, num_dim);
    for (int a = 0; a < num_bodies; a++) {
        for (int i = 0; i < num_dim; i++) {
            p_dot[a][i] = 0.0;
            x_dot[a][i] = 0.0;
        }
    }
            
    // Initialize all velocities and accelerations to zero
    for (int i = 0; i < 2 * array_half; i++)   
        dwdt[i] = 0.0;

    // Add 0PN (Newtonian) terms
    if (params->pn_terms[0]) {
        // Velocities
        for (int a = 0; a < num_bodies; a++)
            for (int i = 0; i < num_dim; i++) {
                dwdt[a * num_dim + i] += w[array_half + a * num_dim + i] / m[a];
                x_dot[a][i] += w[array_half + a * num_dim + i] / m[a];
            }
    
        //Accelerations
        for (int a = 0; a < num_bodies; a++) {
            for (int b = a+1; b < num_bodies; b++) {
                for (int i = 0; i < num_dim; i++) {
                    dwdt[array_half + a * num_dim + i] += -m[a] * m[b] * 1/pow(r[a][b], 2) * n[a][b][i];
                    dwdt[array_half + b * num_dim + i] += -m[a] * m[b] * 1/pow(r[a][b], 2) * n[b][a][i];
                    p_dot[a][i] += -m[a] * m[b] * 1/pow(r[a][b], 2) * n[a][b][i];
                    p_dot[b][i] += -m[a] * m[b] * 1/pow(r[a][b], 2) * n[b][a][i];
                }
            }
        }
    }

    // Add 1PN terms
    if (params->pn_terms[1]) {
        // Velocities
        for (int a = 0; a < num_bodies; a++) {
            for (int i = 0; i < num_dim; i++) {
                x_dot[a][i] += -0.5 * dot_product(p[a], p[a], num_dim) / pow(m[a], 3) * p[a][i];
                dwdt[a * num_dim + i] += -0.5 * dot_product(p[a], p[a], num_dim) / pow(m[a], 3) * p[a][i];
                for (int b = 0; b < num_bodies; b++) {
                    if (b != a) {
                        x_dot[a][i] += -0.5 * 1/r[a][b] * (6 * m[b]/m[a] * p[a][i] - 7 * p[b][i] - dot_product(n[a][b], p[b], num_dim) * n[a][b][i]);
                        dwdt[a * num_dim + i] += -0.5 * 1/r[a][b] * (6 * m[b]/m[a] * p[a][i] - 7 * p[b][i] - dot_product(n[a][b], p[b], num_dim) * n[a][b][i]);
                    }
                }
            }
        }

        // Accelerations
        for (int a = 0; a < num_bodies; a++) {
            for (int i = 0; i < num_dim; i++) {
                for (int b = 0; b < num_bodies; b++) {
                    if (b != a) {
                        p_dot[a][i] += -0.5 / pow(r[a][b], 2) * (3 * m[b]/m[a] * dot_product(p[a], p[a], num_dim) + 3 * m[a]/m[b] * dot_product(p[b], p[b], num_dim) 
                                               - 7 * dot_product(p[a], p[b], num_dim) - 3 * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim)) * n[a][b][i];
                        dwdt[array_half + a * num_dim + i] += -0.5 / pow(r[a][b], 2) * (3 * m[b]/m[a] * dot_product(p[a], p[a], num_dim) + 3 * m[a]/m[b] * dot_product(p[b], p[b], num_dim) 
                                               - 7 * dot_product(p[a], p[b], num_dim) - 3 * dot_product(n[a][b], p[a], num_dim) * dot_product(n[a][b], p[b], num_dim)) * n[a][b][i];
                        p_dot[a][i] += -0.5 / pow(r[a][b], 2) * (dot_product(n[a][b], p[b], num_dim) * p[a][i] + dot_product(n[a][b], p[a], num_dim) * p[b][i]);
                        dwdt[array_half + a * num_dim + i] += -0.5 / pow(r[a][b], 2) * (dot_product(n[a][b], p[b], num_dim) * p[a][i] + dot_product(n[a][b], p[a], num_dim) * p[b][i]);
                    }
                    if (num_bodies > 0) {
                        for (int c = 0; c < num_bodies; c++) {
                            if (b != a && c != a) {
                                p_dot[a][i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[a][c]) * n[a][b][i];
                                dwdt[array_half + a * num_dim + i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[a][c]) * n[a][b][i];
                            }
                            if (b != a && c != b) {
                                p_dot[a][i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[b][c]) * n[a][b][i];
                                dwdt[array_half + a * num_dim + i] += m[a] * m[b] * m[c] / (pow(r[a][b], 2) * r[b][c]) * n[a][b][i];
                            }
                        }
                    }
                }
            }
        }
    }

    // Add 2PN terms (using finite differencing on the hamiltonian)
    if (params->pn_terms[2]) {
        //update_eom_hamiltonian_fd(w, dwdt, H2PN_threebody, 1e-5, params);
        update_eom_hamiltonian_cs(w, dwdt, H2PN_nbody_complex, 1e-20, params);
    }

    // Add 2.5PN terms
    if (params->pn_terms[3]) {
        double** chi_dot, ****dp_chi, ****dx_chi;
        double*** x_rel_dot;

        allocate_2d_array(&chi_dot, num_dim, num_dim);
        allocate_3d_array(&x_rel_dot, num_bodies, num_bodies, num_dim);
        allocate_4d_array(&dp_chi, num_bodies, num_dim, num_dim, num_dim);
        allocate_4d_array(&dx_chi, num_bodies, num_dim, num_dim, num_dim);

        for (int a = 0; a < num_bodies; a++) {
            for (int b = 0; b < num_bodies; b++) {
                for (int i = 0; i < num_dim; i++) {
                    x_rel_dot[a][b][i] = x_dot[a][i] - x_dot[b][i];
                }
            }
        }

        // Initialize chi arrays to zero
        for (int i = 0; i < num_dim; i++) {
            for (int j = 0; j < num_dim; j++) {
                chi_dot[i][j] = 0.0;
            }
        }
        for (int a = 0; a < num_bodies; a++) {
            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        dp_chi[a][i][j][k] = 0.0;
                        dx_chi[a][i][j][k] = 0.0;
                    }
                }
            }
        }    
        // Compute chi values
        for (int i = 0; i < num_dim; i++) {
            for (int j = 0; j < num_dim; j++) {
                for (int a = 0; a < num_bodies; a++) {
                    chi_dot[i][j] += 2.0 / m[a] * (2 * dot_product(p_dot[a], p[a], num_dim) * kronecker_delta(i, j) - 3 * (p_dot[a][i] * p[a][j] + p[a][i] * p_dot[a][j]));
                }
                for (int a = 0; a < num_bodies; a++) {
                    for (int b = 0; b < num_bodies; b++) {
                        if (b != a) {
                            chi_dot[i][j] += m[a] * m[b] / pow(r[a][b], 2) * (3 * (x_rel_dot[a][b][i] * n[a][b][j] + n[a][b][i] * x_rel_dot[a][b][j]) +
                                                                              dot_product(n[a][b], x_rel_dot[a][b], num_dim) * (kronecker_delta(i, j) - 9 * n[a][b][i] * n[a][b][j]));
                        }
                    }
                }
            }
        }

        for (int c = 0; c < num_bodies; c++) {
            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        dp_chi[c][i][j][k] += 2.0 / m[c] * (2 * p[c][k] * kronecker_delta(i, j) - 3 * (p[c][j] * kronecker_delta(i, k) + p[c][i] * kronecker_delta(j, k)));
                    }
                }
            }
        }
        for (int c = 0; c < num_bodies; c++) {
            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        for (int a = 0; a < num_bodies; a++) {
                            for (int b = 0; b < num_bodies; b++) {
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

        for (int a = 0; a < num_bodies; a++) {
            for (int k = 0; k < num_dim; k++) {
                for (int i = 0; i < num_dim; i++) {
                    for (int j = 0; j < num_dim; j++) {
                        dwdt[a * num_dim + k] += 1/45.0 * chi_dot[i][j] * dp_chi[a][i][j][k];
                        dwdt[array_half + a * num_dim + k] += -1/45.0 * chi_dot[i][j] * dx_chi[a][i][j][k];
                    }
                }
            }
        }

        free_2d_array(chi_dot, num_dim);
        free_3d_array(x_rel_dot, num_bodies, num_bodies);
        free_4d_array(dp_chi, num_bodies, num_dim, num_dim);
        free_4d_array(dx_chi, num_bodies, num_dim, num_dim);
    }

    free_vector(m);
    free_2d_array(p, num_bodies);
    free_2d_array(r, num_bodies);
    free_2d_array(p_dot, num_bodies);
    free_2d_array(x_dot, num_bodies);
    free_3d_array(x_rel, num_bodies, num_bodies);
    free_3d_array(n, num_bodies, num_bodies);
}