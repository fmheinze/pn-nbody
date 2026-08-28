#ifndef EOM_INTERNAL_H
#define EOM_INTERNAL_H

#include "cache.h"


// Conservative 0PN and 1PN parts
void add_velocity_0pn(const PairCache *cache, double *velocity, double *lower_xdot);
void add_momentum_0pn(const PairCache *cache, double *momentum, double *lower_pdot);
void add_velocity_1pn(const PairCache *cache, double *velocity, double *lower_xdot);
void add_momentum_1pn(const PairCache *cache, double *momentum, double *lower_pdot);

// Conservative 2PN parts
void add_velocity_2pn_onebody_analytic(const PairCache *cache, double *velocity);
void add_velocity_2pn_pair_analytic(const PairCache *cache, double *velocity);
void add_velocity_2pn_triple_analytic(const PairCache *cache, double *velocity);
void add_momentum_2pn_pair_analytic(const PairCache *cache, double *momentum);
void add_momentum_2pn_triple_analytic(const PairCache *cache, double *momentum);
void add_momentum_2pn_fourbody_analytic(const PairCache *cache, double *momentum);

// Dissipative 2.5PN parts
void pn_compute_25pn_chi_dot(const PairCache *cache, int num_bodies, int num_dim,
    double x_dot[num_bodies][num_dim], double p_dot[num_bodies][num_dim],
    double C[num_dim][num_dim]);
void pn_add_25pn_velocity_contribution(const PairCache *cache, int num_dim,
    double C[num_dim][num_dim], double *velocity);
void pn_add_25pn_momentum_contribution(const PairCache *cache, int num_dim,
    double C[num_dim][num_dim], double *momentum);

#endif
