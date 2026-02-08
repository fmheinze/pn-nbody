/**
 * @file physical_quantities.c
 * @brief Functions for the computation of physical quantities
 *
 * Functions for the computation of physical quantities, such as total energy, ...
 * (add more here if needed).
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "eom.h"
#include "utils.h"
#include "hamiltonians.h"


/**
 * @brief Computes the total energy of a system.
 * 
 * Computes the total energy of a system based on the specified post-Newtonian approximation.
 * 
 * @param[in]   w              Full state vector, w = [positions, momenta]
 * @param[in]   ode_params     Parameter struct containing general information about the system
 * @returns Total energy of the system.
 */
double total_energy_conservative(double* w, struct ode_params* ode_params) 
{
    double H = 0.0;
    if (ode_params->pn_terms[0] == 1)
        H += H0PN(w, ode_params);
    if (ode_params->pn_terms[1] == 1)
        H += H1PN(w, ode_params);
    if (ode_params->pn_terms[2] == 1)
        H += H2PN(w, ode_params, 1);
    return H;
}