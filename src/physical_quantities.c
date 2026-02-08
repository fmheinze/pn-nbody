#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pn_eom.h"
#include "utils.h"
#include "pn_eom_hamiltonians.h"


double total_energy_conservative(double* w, struct ode_params* params) 
{
    double H = 0.0;
    if (params->pn_terms[0] == 1)
        H += H0PN(w, params);
    if (params->pn_terms[1] == 1)
        H += H1PN(w, params);
    if (params->pn_terms[2] == 1)
        H += H2PN_nbody(w, params, 0, 1);
    return H;
}