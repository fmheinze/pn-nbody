#ifndef CACHE_H
#define CACHE_H

#include <stddef.h>
#include "active_list.h"

#define PAIR_CACHE_TRIPLE_PARALLEL_MIN_BODIES 10


// ------------------------------------------------------------------------------------------------
// Pair and body geometry cache
// ------------------------------------------------------------------------------------------------

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

    double *triple_accum;
    int triple_accum_threads;
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


// ------------------------------------------------------------------------------------------------
// UTT4 evaluation cache
// ------------------------------------------------------------------------------------------------

#define UTT4_ORDER_MEMORY_MAX_ENTRIES ((size_t)8u << 20)


/**
 * @brief Cached logarithmic UTT4 value and gradient for one exact position state.
 *
 * The cache key stores active-body positions and masses, the active set, and all numerical
 * integral settings that can affect the result. The cached value and gradient are the complete
 * mass-weighted logarithmic UTT4 contribution, including the Hamiltonian symmetry factor 1/4.
 */
typedef struct {
    int valid;

    // Exact key for the position-only, mass-weighted logarithmic UTT4 contribution.
    double *positions;      // [N*D], only active-body entries are relevant
    double *masses;         // [N]
    int *active;            // [N]

    // Numerical-integral settings are part of the key because they can change the result.
    double epsrel;
    double epsabs;
    int min_order;
    int max_order;
    int adaptive;
    int max_depth;
    int parallel;

    double value;
    double *grad;           // [N*D]
} UTT4LnCache;


/**
 * @brief Per-quadruple memory of the quadrature order that last met the UTT4 tolerance.
 *
 * The same four bodies are integrated many times per step at slowly drifting geometries, and the
 * order they need drifts just as slowly. Remembering it lets the evaluator skip the comparison
 * sweep on most evaluations and start from the cheapest order already known to work, instead of
 * rebuilding both facts from the conservative tolerance-derived guess every time.
 *
 * Entries are addressed by body slot, not by position in the active list, so a merger cannot
 * silently repoint an entry at a different quadruple. The active set is stored alongside and the
 * memory is cleared whenever it changes, because a reused slot holds a different body.
 *
 * age counts evaluations since the entry was last verified. Once it reaches the verify interval
 * the next evaluation runs the full two-order check, which both confirms the order and allows it
 * to move.
 */
typedef struct {
    int enabled;
    size_t count;          // C(num_bodies, 4), or 0 when the memory is disabled
    short *order;          // last order known to meet tolerance, 0 = never evaluated;
    short *age;            // evaluations since the last verified check
    int *active;           // [N], active set this memory was built for
} UTT4OrderMemory;


/**
 * @brief Colexicographic rank of the four-body subset {a < b < c < d}.
 *
 * Addresses the order memory in O(1) without materializing a quadruple list. The rank is
 * C(a,1) + C(b,2) + C(c,3) + C(d,4), which covers 0 .. C(N,4)-1 exactly.
 */
static inline size_t utt4_order_memory_index(int a, int b, int c, int d)
{
    return (size_t)a
         + (size_t)b*(b - 1)/2
         + (size_t)c*(c - 1)*(c - 2)/6
         + (size_t)d*(d - 1)*(d - 2)*(d - 3)/24;
}


// Everything memoized about UTT4 evaluation for one system.
typedef struct UTT4Cache {
    int num_bodies;
    int num_dim;
    int array_half;

    ActiveList active;

    UTT4LnCache ln;
    UTT4OrderMemory order;
} UTT4Cache;


UTT4Cache *utt4_cache_create(const struct ode_params *ode_params);
void utt4_cache_destroy(UTT4Cache *cache);
UTT4Cache *utt4_cache_get_workspace(struct ode_params *ode_params);

#endif
