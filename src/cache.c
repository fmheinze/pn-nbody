/**
 * @file cache.c
 * @brief Persistent workspaces reused across post-Newtonian RHS and Hamiltonian evaluations.
 *
 * Two independent caches live here. PairCache holds the per-state geometry the equations of
 * motion and the conservative Hamiltonians both need: inverse masses, momentum norms and scalar
 * products, pair separations, inverse separations, unit separation vectors, and pair-endpoint
 * momentum contractions. Its storage is allocated once and refreshed in place, with only the
 * levels required by the enabled post-Newtonian terms being recomputed, and masses and momenta
 * retained as non-owning references to the existing parameter and state arrays. Pair symmetries
 * avoid redundant calculations, while arbitrary triple contractions are evaluated on demand
 * rather than stored in an O(N^3) array.
 *
 * UTT4Cache holds what is memoized about the four-body UTT4 term: the logarithmic value and
 * gradient for one exact position state, so energy and force evaluations at the same positions
 * share one numerical quadrature, and the per-quadruple quadrature-order memory. None of that is
 * pair data, which is why it is a separate workspace rather than a member of PairCache.
 *
 * Both are mutable workspaces shared across an evaluation, so a single ode_params instance must
 * not be evaluated concurrently by multiple threads.
 */

#include <math.h>
#include <stdlib.h>
#include "cache.h"
#include "utils.h"


// Checked calloc helper function
static void *xcalloc(size_t count, size_t size)
{
    if (count == 0)
        return NULL;

    void *ptr = calloc(count, size);
    if (ptr == NULL)
        errorexit("Memory allocation failed in pair_cache");
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


// Refresh only the requested levels, without allocating or freeing memory.
void pair_cache_refresh(PairCache *cache, const double *w, const struct ode_params *ode_params,
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

    levels = pair_cache_add_dependencies(levels);
    cache->built_levels = PAIR_CACHE_LEVEL_NONE;
    cache->m = ode_params->masses;
    cache->p = w + cache->array_half;
    active_list_refresh(&cache->active, ode_params);

    if (levels == PAIR_CACHE_LEVEL_NONE)
        return;
    if (cache->m == NULL)
        errorexit("ode_params->masses is NULL");

    const int D = cache->num_dim;
    const ActiveList *active = &cache->active;

    for (int ia = 0; ia < active->num_active; ++ia) {
        const int a = active->ids[ia];
        if (cache->m[a] == 0.0)
            errorexit("Encountered zero active-body mass in pair_cache_refresh");
        cache->inv_m[a] = 1.0 / cache->m[a];
    }

    if (levels & PAIR_CACHE_LEVEL_P2) {
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

    if (levels & PAIR_CACHE_LEVEL_GEOMETRY) {
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

    if (levels & PAIR_CACHE_LEVEL_PAIR_DOTS) {
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

    cache->built_levels = levels;
}


// Determine the minimum cache levels needed by the enabled RHS PN terms.
unsigned int pair_cache_levels_for_rhs(const struct ode_params *ode_params)
{
    if (ode_params == NULL)
        errorexit("pair_cache_levels_for_rhs received NULL ode_params");
    if (ode_params->pn_terms == NULL)
        errorexit("ode_params->pn_terms is NULL");

    unsigned int levels = PAIR_CACHE_LEVEL_NONE;

    if (ode_params->pn_terms[0] || ode_params->pn_terms[1]
        || ode_params->pn_terms[2] || ode_params->pn_terms[3]) {
        levels |= PAIR_CACHE_LEVEL_MASS_MOMENTUM;
        levels |= PAIR_CACHE_LEVEL_GEOMETRY;
    }

    if (ode_params->pn_terms[1] || ode_params->pn_terms[2]) {
        levels |= PAIR_CACHE_LEVEL_P2;
        levels |= PAIR_CACHE_LEVEL_PAIR_DOTS;
    }

    return levels;
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

    cache->ln.positions = (double *)xcalloc(N * D, sizeof(*cache->ln.positions));
    cache->ln.masses    = (double *)xcalloc(N, sizeof(*cache->ln.masses));
    cache->ln.active    = (int *)xcalloc(N, sizeof(*cache->ln.active));
    cache->ln.grad      = (double *)xcalloc(N * D, sizeof(*cache->ln.grad));
    cache->ln.valid = 0;

    cache->order.enabled = 0;
    cache->order.count = 0;
    if (ode_params->include_utt4 && N >= 4) {
        const size_t quadruples = N*(N - 1)*(N - 2)*(N - 3)/24;

        if (quadruples <= UTT4_ORDER_MEMORY_MAX_ENTRIES) {
            UTT4OrderMemory *mem = &cache->order;
            mem->order  = (short *)xcalloc(quadruples, sizeof(*mem->order));
            mem->age    = (short *)xcalloc(quadruples, sizeof(*mem->age));
            mem->active = (int *)xcalloc(N, sizeof(*mem->active));
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
    free(cache->ln.positions);
    free(cache->ln.masses);
    free(cache->ln.active);
    free(cache->ln.grad);
    free(cache->order.order);
    free(cache->order.age);
    free(cache->order.active);
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
