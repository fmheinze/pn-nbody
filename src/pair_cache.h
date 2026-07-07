#ifndef PAIR_CACHE_H
#define PAIR_CACHE_H

#include "active_list.h"
#include "eom.h"

/**
 * @brief Cached geometric and momentum invariants for one RHS evaluation.
 *
 * This is infrastructure for the analytic 2PN RHS. It does not change the
 * equations of motion by itself. The expensive analytic terms can later reuse:
 *
 *   r_ab, 1/r_ab, 1/r_ab^2, n_ab,
 *   p_a, p_a^2, p_a.p_b, n_ab.p_c.
 *
 * Arrays are stored flat for contiguous memory access.
 */
typedef struct {
    int num_bodies;
    int num_dim;
    int array_half;

    ActiveList active;

    double *m;        // [N]
    double *inv_m;    // [N]

    double *p;        // [N*D]
    double *p2;       // [N]
    double *p_dot;    // [N*N], p_a . p_b

    double *x_rel;    // [N*N*D], x_a - x_b
    double *r;        // [N*N]
    double *inv_r;    // [N*N]
    double *inv_r2;   // [N*N]
    double *n;        // [N*N*D], n_ab

    double *n_dot_p;  // [N*N*N], n_ab . p_c
} PairCache;


static inline int pair_cache_idx2(const PairCache *cache, int a, int b)
{
    return a * cache->num_bodies + b;
}


static inline int pair_cache_idx3(const PairCache *cache, int a, int b, int i)
{
    return (a * cache->num_bodies + b) * cache->num_dim + i;
}


static inline int pair_cache_idx_ndotp(const PairCache *cache, int a, int b, int c)
{
    return (a * cache->num_bodies + b) * cache->num_bodies + c;
}


static inline double pair_cache_p(const PairCache *cache, int a, int i)
{
    return cache->p[a * cache->num_dim + i];
}


static inline double pair_cache_xrel(const PairCache *cache, int a, int b, int i)
{
    return cache->x_rel[pair_cache_idx3(cache, a, b, i)];
}


static inline double pair_cache_n(const PairCache *cache, int a, int b, int i)
{
    return cache->n[pair_cache_idx3(cache, a, b, i)];
}


static inline double pair_cache_r(const PairCache *cache, int a, int b)
{
    return cache->r[pair_cache_idx2(cache, a, b)];
}


static inline double pair_cache_inv_r(const PairCache *cache, int a, int b)
{
    return cache->inv_r[pair_cache_idx2(cache, a, b)];
}


static inline double pair_cache_p_dot(const PairCache *cache, int a, int b)
{
    return cache->p_dot[pair_cache_idx2(cache, a, b)];
}


static inline double pair_cache_n_dot_p(const PairCache *cache, int a, int b, int c)
{
    return cache->n_dot_p[pair_cache_idx_ndotp(cache, a, b, c)];
}


/**
 * @brief Build all cached quantities from the current state vector.
 *
 * Call pair_cache_free after use.
 */
void pair_cache_build(PairCache *cache, const double *w, const struct ode_params *ode_params);

/**
 * @brief Release all memory owned by the cache.
 */
void pair_cache_free(PairCache *cache);

#endif
