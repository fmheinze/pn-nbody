#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <complex.h>


struct ode_params;
struct PairCache;


typedef complex double (*c_hamiltonian)(complex double*, struct ode_params*, int);

void utt4_ln_integral_cached(double *w, struct ode_params *ode_params, double *value, double *grad);


// Standard variants that refresh the cache before evaluating the Hamiltonian
double H0PN(double* w, struct ode_params* ode_params);
double H1PN(double* w, struct ode_params* ode_params);
double H2PN(double* w, struct ode_params* ode_params, int utt4_flag);


// Cache-aware variants used after one shared refresh by the energy evaluator
double H0PN_cached(const struct PairCache *cache);
double H1PN_cached(const struct PairCache *cache);
double H2PN_cached(double* w, struct ode_params* ode_params, const struct PairCache *cache, int utt4_flag);

#endif
