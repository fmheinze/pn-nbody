#ifndef CACHE_H
#define CACHE_H

#include <stddef.h>
#include "active_list.h"

#define PAIR_CACHE_TRIPLE_PARALLEL_MIN_BODIES 10


// ------------------------------------------------------------------------------------------------
// Shared state key
// ------------------------------------------------------------------------------------------------

/**
 * @brief The one record of "what state were the caches built for".
 *
 * Every cache below is valid only for an exact phase-space state, mass set, active set and the
 * settings that can change a derived quantity. Rather than each keeping its own copy of all that
 * and comparing it independently, they share this key and remember only the generation counter it
 * carried when they were filled. A cache is valid exactly while its stored counter still matches.
 *
 * Three counters are kept because the caches depend on different parts of the state:
 *   generation           bumped by any change at all -- used by PairCache and DynamicsCache;
 *   position_generation  bumped when positions, masses, the active set or the UTT4 integral
 *                        settings change, but not when only momenta move -- used by the UTT4
 *                        logarithmic result, which is a function of positions alone;
 *   active_generation    bumped only when the active set changes -- used by the UTT4 order
 *                        memory, which must survive ordinary motion and be dropped on a merger.
 *
 * Counters only ever increase, so a stale cache can never be mistaken for a fresh one.
 */
typedef struct StateKey {
    int num_bodies;
    int num_dim;
    int valid;

    double *state;      // [2*N*D]
    double *masses;     // [N]
    int *active;        // [N]

    // Settings that alter a derived quantity
    int pn_terms[4];
    int use_impulse_method;
    int include_utt4;
    double utt4_epsrel;
    double utt4_epsabs;
    int utt4_min_order;
    int utt4_max_order;
    int utt4_adaptive;
    int utt4_max_depth;
    int utt4_parallel;
    int utt4_verify_interval;

    unsigned long long generation;
    unsigned long long position_generation;
    unsigned long long active_generation;
} StateKey;


StateKey *state_key_sync(struct ode_params *ode_params, const double *w);

void state_key_destroy(StateKey *key);


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
 * Owned arrays are allocated once by pair_cache_create and refreshed in place. Levels already
 * built for an exactly matching state, masses, and active set are reused by subsequent physical
 * quantities. The masses and current momenta are non-owning references. A single ode_params
 * instance must not be evaluated concurrently by multiple threads because refreshes intentionally
 * reuse this mutable workspace.
 */
typedef struct PairCache {
    int num_bodies;
    int num_dim;
    int array_half;
    unsigned int built_levels;

    // Generation of the shared state key these levels were built for.
    unsigned long long generation;

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

    double *triple_accum;   // [N*N*D], shared momentum/velocity triple accumulator
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

void pair_cache_refresh(PairCache *cache, const double *w, struct ode_params *ode_params,
    unsigned int levels);


// ------------------------------------------------------------------------------------------------
// UTT4 evaluation cache
// ------------------------------------------------------------------------------------------------

#define UTT4_ORDER_MEMORY_MAX_ENTRIES ((size_t)8u << 20)


/**
 * @brief Cached logarithmic UTT4 value and gradient for one exact position state.
 *
 * Valid while the shared key's position generation is unchanged: the logarithmic contribution is
 * a function of the active bodies' positions and masses alone, so it survives a step that moves
 * only momenta. The stored value and gradient are the complete mass-weighted contribution,
 * including the Hamiltonian symmetry factor 1/4.
 */
typedef struct {
    unsigned long long position_generation;
    int valid;

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
 * silently repoint an entry at a different quadruple. The active generation of the shared key is
 * recorded alongside, and the memory is cleared whenever it moves, because a reused slot holds a
 * different body.
 *
 * age counts evaluations since the entry was last verified. Once it reaches the verify interval
 * the next evaluation runs the full two-order check, which both confirms the order and allows it
 * to move.
 */
typedef struct {
    int enabled;
    size_t count;          // C(num_bodies, 4), or 0 when the memory is disabled
    short *order;          // last order known to meet tolerance, 0 = never evaluated
    short *age;            // evaluations since the last verified check
    unsigned long long active_generation;   // active set the entries were recorded for
} UTT4OrderMemory;


/**
 * @brief Colexicographic rank of the four-body subset {a < b < c < d}.
 *
 * Addresses the order memory in O(1) without materializing a quadruple list. The rank is
 * C(a,1) + C(b,2) + C(c,3) + C(d,4), which covers 0 .. C(N,4)-1 exactly.
 */
static inline size_t utt4_order_memory_index(int a, int b, int c, int d)
{
    const size_t aa = (size_t)a;
    const size_t bb = (size_t)b;
    const size_t cc = (size_t)c;
    const size_t dd = (size_t)d;

    return aa
         + bb*(bb - 1)/2
         + cc*(cc - 1)*(cc - 2)/6
         + dd*(dd - 1)*(dd - 2)*(dd - 3)/24;
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


// ------------------------------------------------------------------------------------------------
// Dynamics cache
// ------------------------------------------------------------------------------------------------

/**
 * @brief Results of a complete post-Newtonian evaluation at one exact state.
 *
 * The coordinate velocity dx/dt and the momentum derivative dp/dt come out of the same machinery,
 * but a caller may want only one of them: writing velocities to file needs no forces, and a
 * right-hand side evaluated just after a velocity output needs no second velocity pass. Each half
 * is therefore flagged separately, and the 2.5PN quadrupole derivative -- which both halves need
 * and neither owns -- is kept alongside them.
 *
 * Validity is the shared key's generation, so any change of state, masses, active set or the
 * settings that reach the equations of motion drops the whole entry.
 */
typedef struct DynamicsCache {
    int num_bodies;
    int num_dim;
    int array_half;

    unsigned long long generation;

    int velocity_valid;
    int rhs_valid;
    int chi_dot_valid;

    double *velocity;        // [N*D]
    double *rhs;             // [2*N*D]
    double *scratch;         // lazily allocated [2*N*D] 2.5PN/UTT4 scratch
    double chi_dot[3][3];
} DynamicsCache;


DynamicsCache *dynamics_cache_create(const struct ode_params *ode_params);

void dynamics_cache_destroy(DynamicsCache *cache);

DynamicsCache *dynamics_cache_prepare(struct ode_params *ode_params, const double *w);

double *dynamics_cache_get_scratch(DynamicsCache *cache);

#endif
