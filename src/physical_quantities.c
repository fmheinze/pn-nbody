/**
 * @file physical_quantities.c
 * @brief Functions for the computation of physical quantities
 *
 * Functions for the computation of physical quantities, such as total energy, ...
 */

#include "eom.h"
#include "hamiltonian.h"
#include "pair_cache.h"


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
    unsigned int levels = PAIR_CACHE_LEVEL_NONE;

    if (ode_params->pn_terms[0])
        levels |= PAIR_CACHE_LEVEL_GEOMETRY | PAIR_CACHE_LEVEL_P2;

    if (ode_params->pn_terms[1] || ode_params->pn_terms[2])
        levels |= PAIR_CACHE_LEVEL_GEOMETRY | PAIR_CACHE_LEVEL_P2 | PAIR_CACHE_LEVEL_PAIR_DOTS;

    PairCache *cache = pair_cache_get_workspace(ode_params);
    pair_cache_refresh(cache, w, ode_params, levels);

    double H = 0.0;
    if (ode_params->pn_terms[0] == 1)
        H += H0PN_cached(cache);

    if (ode_params->pn_terms[1] == 1)
        H += H1PN_cached(cache);
        
    if (ode_params->pn_terms[2] == 1) {
        if (ode_params->include_utt4)
            H += H2PN_cached(w, ode_params, cache, 1);
        else
            H += H2PN_cached(w, ode_params, cache, 0);
    }
    return H;
}
