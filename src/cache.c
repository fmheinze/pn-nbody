/**
 * @file cache.c
 * @brief Persistent workspaces reused across post-Newtonian evaluations, and the key they share.
 *
 * StateKey is the single record of what state the caches below were built for: the phase-space
 * state, masses, active set and the settings that can change a derived quantity. Each cache
 * stores only the generation counter the key carried when it was filled, so validity is an
 * integer comparison and there is one definition of "the same state" rather than one per cache.
 *
 * PairCache holds the per-state geometry the equations of motion and the conservative
 * Hamiltonians both need: inverse masses, momentum norms and scalar products, pair separations,
 * inverse separations, unit separation vectors, and pair-endpoint momentum contractions. Its
 * storage is allocated once and refreshed in place, and levels already built at the current state
 * are kept while missing ones are added incrementally. Masses and momenta are non-owning
 * references into the existing parameter and state arrays. Pair symmetries avoid redundant
 * calculations, while arbitrary triple contractions are evaluated on demand rather than stored in
 * an O(N^3) array.
 *
 * UTT4Cache holds what is memoized about the four-body UTT4 term: the logarithmic value and
 * gradient, keyed on the position generation because that contribution does not depend on
 * momenta, and the per-quadruple quadrature-order memory, keyed on the active generation so it
 * survives ordinary motion and is dropped on a merger.
 *
 * DynamicsCache holds the results of a full evaluation -- coordinate velocities, the complete
 * right-hand side, and the 2.5PN quadrupole derivative -- so that writing velocities and stepping
 * the integrator at the same state do not repeat each other's work.
 *
 * All are mutable workspaces shared across an evaluation, so a single ode_params instance must
 * not be evaluated concurrently by multiple threads.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cache.h"
#include "utils.h"


static void *xcalloc(size_t count, size_t size);


// ------------------------------------------------------------------------------------------------
// Shared state key
// ------------------------------------------------------------------------------------------------

static StateKey *state_key_create(const struct ode_params *ode_params)
{
    StateKey *key = (StateKey *)calloc(1, sizeof(*key));
    if (key == NULL)
        errorexit("Memory allocation failed for StateKey");

    const size_t N = (size_t)ode_params->num_bodies_initial;
    const size_t half = N*(size_t)ode_params->num_dim;

    key->num_bodies = ode_params->num_bodies_initial;
    key->num_dim = ode_params->num_dim;
    key->state = (double *)xcalloc(2*half, sizeof(*key->state));
    key->masses = (double *)xcalloc(N, sizeof(*key->masses));
    key->active = (int *)xcalloc(N, sizeof(*key->active));
    return key;
}


void state_key_destroy(StateKey *key)
{
    if (key == NULL)
        return;

    free(key->state);
    free(key->masses);
    free(key->active);
    free(key);
}


// Settings that reach the UTT4 logarithmic integral, and therefore any position-keyed result.
static int state_key_utt4_settings_match(const StateKey *key, const struct ode_params *p)
{
    return key->include_utt4 == p->include_utt4
        && key->utt4_epsrel == p->utt4_epsrel
        && key->utt4_epsabs == p->utt4_epsabs
        && key->utt4_min_order == p->utt4_min_order
        && key->utt4_max_order == p->utt4_max_order
        && key->utt4_adaptive == p->utt4_adaptive
        && key->utt4_max_depth == p->utt4_max_depth
        && key->utt4_parallel == p->utt4_parallel
        && key->utt4_verify_interval == p->utt4_verify_interval;
}


StateKey *state_key_sync(struct ode_params *ode_params, const double *w)
{
    if (ode_params == NULL)
        errorexit("state_key_sync received NULL ode_params");
    if (w == NULL)
        errorexit("state_key_sync received NULL state vector");
    if (ode_params->num_bodies_initial > 0
        && (ode_params->masses == NULL || ode_params->active == NULL))
        errorexit("ode_params masses or active list is NULL");

    if (ode_params->state_key == NULL)
        ode_params->state_key = state_key_create(ode_params);

    StateKey *key = ode_params->state_key;
    if (key->num_bodies != ode_params->num_bodies_initial || key->num_dim != ode_params->num_dim)
        errorexit("StateKey dimensions do not match ode_params");

    const size_t N = (size_t)key->num_bodies;
    const size_t half = N*(size_t)key->num_dim;
    const size_t pos_bytes = half*sizeof(*key->state);
    const size_t mom_bytes = half*sizeof(*key->state);

    const int positions_same = key->valid
        && memcmp(key->state, w, pos_bytes) == 0
        && memcmp(key->masses, ode_params->masses, N*sizeof(*key->masses)) == 0;
    const int active_same = key->valid
        && memcmp(key->active, ode_params->active, N*sizeof(*key->active)) == 0;
    const int momenta_same = key->valid
        && memcmp(key->state + half, w + half, mom_bytes) == 0;
    const int utt4_same = key->valid && state_key_utt4_settings_match(key, ode_params);
    const int dynamics_same = key->valid
        && memcmp(key->pn_terms, ode_params->pn_terms, sizeof(key->pn_terms)) == 0
        && key->use_impulse_method == ode_params->use_impulse_method;

    if (positions_same && active_same && momenta_same && utt4_same && dynamics_same)
        return key;

    if (!active_same)
        ++key->active_generation;
    if (!positions_same || !active_same || !utt4_same)
        ++key->position_generation;
    ++key->generation;

    memcpy(key->state, w, pos_bytes + mom_bytes);
    memcpy(key->masses, ode_params->masses, N*sizeof(*key->masses));
    memcpy(key->active, ode_params->active, N*sizeof(*key->active));
    memcpy(key->pn_terms, ode_params->pn_terms, sizeof(key->pn_terms));
    key->use_impulse_method = ode_params->use_impulse_method;
    key->include_utt4 = ode_params->include_utt4;
    key->utt4_epsrel = ode_params->utt4_epsrel;
    key->utt4_epsabs = ode_params->utt4_epsabs;
    key->utt4_min_order = ode_params->utt4_min_order;
    key->utt4_max_order = ode_params->utt4_max_order;
    key->utt4_adaptive = ode_params->utt4_adaptive;
    key->utt4_max_depth = ode_params->utt4_max_depth;
    key->utt4_parallel = ode_params->utt4_parallel;
    key->utt4_verify_interval = ode_params->utt4_verify_interval;
    key->valid = 1;
    return key;
}


// Checked calloc helper function
static void *xcalloc(size_t count, size_t size)
{
    if (count == 0)
        return NULL;

    void *ptr = calloc(count, size);
    if (ptr == NULL)
        errorexit("Memory allocation failed for cache workspace");
    return ptr;
}


// Add level dependencies helper function
static unsigned int pair_cache_add_dependencies(unsigned int levels)
{
    if (levels & PAIR_CACHE_LEVEL_PAIR_DOTS)
        levels |= PAIR_CACHE_LEVEL_GEOMETRY;

    if (levels != PAIR_CACHE_LEVEL_NONE)
        levels |= PAIR_CACHE_LEVEL_MASS_MOMENTUM;

    return levels;
}


// Allocate a persistent cache workspace sized for ode_params.
PairCache *pair_cache_create(const struct ode_params *ode_params)
{
    if (ode_params == NULL)
        errorexit("pair_cache_create received NULL ode_params");
    if (ode_params->num_bodies_initial < 0)
        errorexit("num_bodies_initial must be non-negative");
    if (ode_params->num_dim <= 0)
        errorexit("num_dim must be positive");

    PairCache *cache = (PairCache *)calloc(1, sizeof(*cache));
    if (cache == NULL)
        errorexit("Memory allocation failed for PairCache workspace");

    const size_t N = (size_t)ode_params->num_bodies_initial;
    const size_t D = (size_t)ode_params->num_dim;
    const size_t N2 = N * N;


    cache->num_bodies = ode_params->num_bodies_initial;
    cache->num_dim = ode_params->num_dim;
    cache->array_half = cache->num_bodies * cache->num_dim;

    active_list_init(&cache->active, ode_params);

    cache->inv_m     = (double *)xcalloc(N, sizeof(*cache->inv_m));
    cache->p2        = (double *)xcalloc(N, sizeof(*cache->p2));
    cache->p_dot     = (double *)xcalloc(N2, sizeof(*cache->p_dot));
    cache->r         = (double *)xcalloc(N2, sizeof(*cache->r));
    cache->inv_r     = (double *)xcalloc(N2, sizeof(*cache->inv_r));
    cache->n         = (double *)xcalloc(N2 * D, sizeof(*cache->n));
    cache->n_dot_p_a = (double *)xcalloc(N2, sizeof(*cache->n_dot_p_a));
    cache->n_dot_p_b = (double *)xcalloc(N2, sizeof(*cache->n_dot_p_b));

    // Momentum or velocity scratch: one half-state accumulator row per outer body.
    cache->triple_accum = (double *)xcalloc(N*(size_t)cache->array_half,
        sizeof(*cache->triple_accum));
    cache->triple_accum_threads = 1;
#ifdef _OPENMP
    if (cache->num_bodies >= PAIR_CACHE_TRIPLE_PARALLEL_MIN_BODIES) {
        const int threads = omp_get_max_threads();
        if (threads > 1)
            cache->triple_accum_threads = threads;
    }
#endif

    cache->m = ode_params->masses;
    return cache;
}


// Release a workspace created by pair_cache_create.
void pair_cache_destroy(PairCache *cache)
{
    if (cache == NULL)
        return;

    active_list_free(&cache->active);
    free(cache->inv_m);
    free(cache->p2);
    free(cache->p_dot);
    free(cache->r);
    free(cache->inv_r);
    free(cache->n);
    free(cache->n_dot_p_a);
    free(cache->n_dot_p_b);
    free(cache->triple_accum);
    free(cache);
}


// Return the pair cache workspace.
PairCache *pair_cache_get_workspace(struct ode_params *ode_params)
{
    if (ode_params == NULL)
        errorexit("pair_cache_get_workspace received NULL ode_params");

    if (ode_params->pair_cache == NULL)
        ode_params->pair_cache = pair_cache_create(ode_params);

    return ode_params->pair_cache;
}


// Refresh only missing requested levels for the exact state, without allocating or freeing memory.
void pair_cache_refresh(PairCache *cache, const double *w, struct ode_params *ode_params,
    unsigned int levels)
{
    if (cache == NULL)
        errorexit("pair_cache_refresh received NULL cache");
    if (w == NULL)
        errorexit("pair_cache_refresh received NULL state vector");
    if (ode_params == NULL)
        errorexit("pair_cache_refresh received NULL ode_params");
    if (cache->num_bodies != ode_params->num_bodies_initial || cache->num_dim != ode_params->num_dim)
        errorexit("PairCache dimensions do not match ode_params");
    if (levels & ~((unsigned int)PAIR_CACHE_LEVEL_FULL))
        errorexit("pair_cache_refresh received an invalid cache level");
    if (cache->num_bodies > 0 && ode_params->masses == NULL)
        errorexit("ode_params->masses is NULL");
    if (cache->num_bodies > 0 && ode_params->active == NULL)
        errorexit("ode_params->active is NULL");

    levels = pair_cache_add_dependencies(levels);
    cache->m = ode_params->masses;
    cache->p = w + cache->array_half;

    const StateKey *key = state_key_sync(ode_params, w);

    if (cache->generation != key->generation) {
        cache->generation = key->generation;
        cache->built_levels = PAIR_CACHE_LEVEL_NONE;
        active_list_refresh(&cache->active, ode_params);
    }

    const unsigned int missing_levels = levels & ~cache->built_levels;

    if (missing_levels == PAIR_CACHE_LEVEL_NONE)
        return;

    const int D = cache->num_dim;
    const ActiveList *active = &cache->active;

    if (missing_levels & PAIR_CACHE_LEVEL_MASS_MOMENTUM) {
        for (int ia = 0; ia < active->num_active; ++ia) {
            const int a = active->ids[ia];
            if (cache->m[a] == 0.0)
                errorexit("Encountered zero active-body mass in pair_cache_refresh");
            cache->inv_m[a] = 1.0 / cache->m[a];
        }
    }

    if (missing_levels & PAIR_CACHE_LEVEL_P2) {
        for (int ia = 0; ia < active->num_active; ++ia) {
            const int a = active->ids[ia];
            double p2 = 0.0;
            for (int i = 0; i < D; ++i) {
                const double pai = pair_cache_p(cache, a, i);
                p2 += pai * pai;
            }
            cache->p2[a] = p2;
        }
    }

    if (missing_levels & PAIR_CACHE_LEVEL_GEOMETRY) {
        for (int ia = 0; ia < active->num_active; ++ia) {
            const int a = active->ids[ia];

            for (int ib = ia + 1; ib < active->num_active; ++ib) {
                const int b = active->ids[ib];
                double dx[D];
                double r2 = 0.0;

                for (int i = 0; i < D; ++i) {
                    dx[i] = w[a * D + i] - w[b * D + i];
                    r2 += dx[i] * dx[i];
                }

                if (r2 <= 0.0)
                    errorexit("Encountered zero separation in pair_cache_refresh");

                const double r = sqrt(r2);
                const double inv_r = 1.0 / r;
                const int ab = pair_cache_idx2(cache, a, b);
                const int ba = pair_cache_idx2(cache, b, a);

                cache->r[ab] = r;
                cache->r[ba] = r;
                cache->inv_r[ab] = inv_r;
                cache->inv_r[ba] = inv_r;

                for (int i = 0; i < D; ++i) {
                    const double ni = dx[i] * inv_r;
                    cache->n[pair_cache_idx3(cache, a, b, i)] = ni;
                    cache->n[pair_cache_idx3(cache, b, a, i)] = -ni;
                }
            }
        }
    }

    if (missing_levels & PAIR_CACHE_LEVEL_PAIR_DOTS) {
        for (int ia = 0; ia < active->num_active; ++ia) {
            const int a = active->ids[ia];

            for (int ib = ia; ib < active->num_active; ++ib) {
                const int b = active->ids[ib];
                double pdot = 0.0;

                for (int i = 0; i < D; ++i)
                    pdot += pair_cache_p(cache, a, i) * pair_cache_p(cache, b, i);

                cache->p_dot[pair_cache_idx2(cache, a, b)] = pdot;
                cache->p_dot[pair_cache_idx2(cache, b, a)] = pdot;
            }
        }

        for (int ia = 0; ia < active->num_active; ++ia) {
            const int a = active->ids[ia];

            for (int ib = ia + 1; ib < active->num_active; ++ib) {
                const int b = active->ids[ib];
                double n_dot_pa = 0.0;
                double n_dot_pb = 0.0;

                for (int i = 0; i < D; ++i) {
                    const double ni = pair_cache_n(cache, a, b, i);
                    n_dot_pa += ni * pair_cache_p(cache, a, i);
                    n_dot_pb += ni * pair_cache_p(cache, b, i);
                }

                const int ab = pair_cache_idx2(cache, a, b);
                const int ba = pair_cache_idx2(cache, b, a);
                cache->n_dot_p_a[ab] = n_dot_pa;
                cache->n_dot_p_b[ab] = n_dot_pb;
                cache->n_dot_p_a[ba] = -n_dot_pb;
                cache->n_dot_p_b[ba] = -n_dot_pa;
            }
        }
    }

    cache->built_levels |= levels;
}


// ------------------------------------------------------------------------------------------------
// UTT4 evaluation cache
// ------------------------------------------------------------------------------------------------

/**
 * @brief Allocate the UTT4 evaluation cache for ode_params.
 *
 * The order memory is only allocated when UTT4 is actually evaluated, and only up to a bounded
 * size: C(N,4) grows steeply, so past the cap the memory is left disabled and every evaluation
 * takes the fully verified path.
 */
UTT4Cache *utt4_cache_create(const struct ode_params *ode_params)
{
    if (ode_params == NULL)
        errorexit("utt4_cache_create received NULL ode_params");
    if (ode_params->num_bodies_initial < 0)
        errorexit("num_bodies_initial must be non-negative");
    if (ode_params->num_dim <= 0)
        errorexit("num_dim must be positive");

    UTT4Cache *cache = (UTT4Cache *)calloc(1, sizeof(*cache));
    if (cache == NULL)
        errorexit("Memory allocation failed for UTT4Cache workspace");

    const size_t N = (size_t)ode_params->num_bodies_initial;
    const size_t D = (size_t)ode_params->num_dim;

    cache->num_bodies = ode_params->num_bodies_initial;
    cache->num_dim = ode_params->num_dim;
    cache->array_half = cache->num_bodies * cache->num_dim;

    active_list_init(&cache->active, ode_params);

    cache->ln.grad = (double *)xcalloc(N * D, sizeof(*cache->ln.grad));
    cache->ln.valid = 0;

    cache->order.enabled = 0;
    cache->order.count = 0;
    if (ode_params->include_utt4 && N >= 4) {
        const size_t quadruples = N*(N - 1)*(N - 2)*(N - 3)/24;

        if (quadruples <= UTT4_ORDER_MEMORY_MAX_ENTRIES) {
            UTT4OrderMemory *mem = &cache->order;
            mem->order  = (short *)xcalloc(quadruples, sizeof(*mem->order));
            mem->age    = (short *)xcalloc(quadruples, sizeof(*mem->age));
            mem->count = quadruples;
            mem->enabled = 1;
        }
    }

    return cache;
}


// Release a workspace created by utt4_cache_create.
void utt4_cache_destroy(UTT4Cache *cache)
{
    if (cache == NULL)
        return;

    active_list_free(&cache->active);
    free(cache->ln.grad);
    free(cache->order.order);
    free(cache->order.age);
    free(cache);
}


// Return the UTT4 evaluation cache, allocating it on first use.
UTT4Cache *utt4_cache_get_workspace(struct ode_params *ode_params)
{
    if (ode_params == NULL)
        errorexit("utt4_cache_get_workspace received NULL ode_params");

    if (ode_params->utt4_cache == NULL)
        ode_params->utt4_cache = utt4_cache_create(ode_params);

    return ode_params->utt4_cache;
}


// ------------------------------------------------------------------------------------------------
// Dynamics cache
// ------------------------------------------------------------------------------------------------

DynamicsCache *dynamics_cache_create(const struct ode_params *ode_params)
{
    if (ode_params == NULL)
        errorexit("dynamics_cache_create received NULL ode_params");

    DynamicsCache *cache = (DynamicsCache *)calloc(1, sizeof(*cache));
    if (cache == NULL)
        errorexit("Memory allocation failed for DynamicsCache");

    const size_t N = (size_t)ode_params->num_bodies_initial;
    const size_t half = N*(size_t)ode_params->num_dim;

    cache->num_bodies = ode_params->num_bodies_initial;
    cache->num_dim = ode_params->num_dim;
    cache->array_half = (int)half;
    cache->velocity = (double *)xcalloc(half, sizeof(*cache->velocity));
    cache->rhs = (double *)xcalloc(2*half, sizeof(*cache->rhs));
    return cache;
}


void dynamics_cache_destroy(DynamicsCache *cache)
{
    if (cache == NULL)
        return;

    free(cache->velocity);
    free(cache->rhs);
    free(cache->scratch);
    free(cache);
}


DynamicsCache *dynamics_cache_prepare(struct ode_params *ode_params, const double *w)
{
    if (ode_params == NULL)
        errorexit("dynamics_cache_prepare received NULL ode_params");

    if (ode_params->dynamics_cache == NULL)
        ode_params->dynamics_cache = dynamics_cache_create(ode_params);

    DynamicsCache *cache = ode_params->dynamics_cache;
    if (cache->num_bodies != ode_params->num_bodies_initial
        || cache->num_dim != ode_params->num_dim)
        errorexit("DynamicsCache dimensions do not match ode_params");

    const StateKey *key = state_key_sync(ode_params, w);

    if (cache->generation != key->generation) {
        cache->generation = key->generation;
        cache->velocity_valid = 0;
        cache->rhs_valid = 0;
        cache->chi_dot_valid = 0;
    }

    return cache;
}


// Lazily allocate scratch used only by 2.5PN lower-order derivatives or a UTT4 gradient.
double *dynamics_cache_get_scratch(DynamicsCache *cache)
{
    if (cache == NULL)
        errorexit("dynamics_cache_get_scratch received NULL cache");

    if (cache->scratch == NULL)
        cache->scratch = (double *)xcalloc((size_t)2 * (size_t)cache->array_half,
            sizeof(*cache->scratch));

    return cache->scratch;
}
