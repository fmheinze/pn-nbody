#ifndef PAIR_CACHE_H
#define PAIR_CACHE_H

#include "active_list.h"


/**
 * @brief Independently refreshable parts of the persistent pair cache.
 *
 * Dependencies are added automatically by pair_cache_refresh. In particular, pair momentum
 * invariants imply geometry and all non-empty levels imply the mass/momentum-reference level.
 */
enum PairCacheLevel {
    PAIR_CACHE_LEVEL_NONE          = 0u,
    PAIR_CACHE_LEVEL_MASS_MOMENTUM = 1u << 0,
    PAIR_CACHE_LEVEL_GEOMETRY      = 1u << 1,
    PAIR_CACHE_LEVEL_P2            = 1u << 2,
    PAIR_CACHE_LEVEL_PAIR_DOTS     = 1u << 3,
    PAIR_CACHE_LEVEL_FULL = PAIR_CACHE_LEVEL_MASS_MOMENTUM | PAIR_CACHE_LEVEL_GEOMETRY |
        PAIR_CACHE_LEVEL_P2 | PAIR_CACHE_LEVEL_PAIR_DOTS
};


/**
 * @brief Persistent RHS workspace for pair geometry and momentum invariants.
 *
 * Owned arrays are allocated once by pair_cache_create and refreshed in place. The masses and
 * current momenta are references to ode_params->masses and the momentum half of the current state
 * vector, respectively. A single ode_params instance must not be evaluated concurrently by
 * multiple threads because refreshes intentionally reuse this mutable workspace.
 */
typedef struct PairCache {
    int num_bodies;
    int num_dim;
    int array_half;
    unsigned int built_levels;

    ActiveList active;

    const double *m;       // [N], non-owning reference
    double *inv_m;         // [N]

    const double *p;       // [N*D], non-owning reference
    double *p2;            // [N]
    double *p_dot;         // [N*N], p_a . p_b

    double *r;             // [N*N]
    double *inv_r;         // [N*N]
    double *n;             // [N*N*D], n_ab

    double *n_dot_p_a;     // [N*N], n_ab . p_a
    double *n_dot_p_b;     // [N*N], n_ab . p_b
} PairCache;


static inline int pair_cache_idx2(const PairCache *cache, int a, int b)
{
    return a * cache->num_bodies + b;
}


static inline int pair_cache_idx3(const PairCache *cache, int a, int b, int i)
{
    return (a * cache->num_bodies + b) * cache->num_dim + i;
}


static inline double pair_cache_p(const PairCache *cache, int a, int i)
{
    return cache->p[a * cache->num_dim + i];
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


/**
 * @brief Return n_ab . p_c without an O(N^3) cache.
 *
 * The two common pair-endpoint cases use cached O(N^2) arrays. Arbitrary c is
 * evaluated on demand in O(D), which is cheaper than maintaining N^3 storage.
 */
static inline double pair_cache_n_dot_p(const PairCache *cache, int a, int b, int c)
{
    const int idx = pair_cache_idx2(cache, a, b);

    if (cache->built_levels & PAIR_CACHE_LEVEL_PAIR_DOTS) {
        if (c == a)
            return cache->n_dot_p_a[idx];
        if (c == b)
            return cache->n_dot_p_b[idx];
    }

    double result = 0.0;
    for (int i = 0; i < cache->num_dim; ++i)
        result += pair_cache_n(cache, a, b, i) * pair_cache_p(cache, c, i);
    return result;
}


PairCache *pair_cache_create(const struct ode_params *ode_params);
void pair_cache_destroy(PairCache *cache);
PairCache *pair_cache_get_workspace(struct ode_params *ode_params);
void pair_cache_refresh(PairCache *cache, const double *w, const struct ode_params *ode_params,
    unsigned int levels);
unsigned int pair_cache_levels_for_rhs(const struct ode_params *ode_params);

#endif
