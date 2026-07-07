/**
 * @file pair_cache.c
 * @brief Cached pair geometry and scalar products for PN RHS evaluations.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pair_cache.h"
#include "utils.h"


static void *xcalloc(size_t count, size_t size)
{
    void *ptr = calloc(count, size);
    if (ptr == NULL)
        errorexit("Memory allocation failed in pair_cache");
    return ptr;
}


static void pair_cache_zero(PairCache *cache)
{
    memset(cache, 0, sizeof(*cache));
}


void pair_cache_build(PairCache *cache, const double *w, const struct ode_params *ode_params)
{
    if (cache == NULL)
        errorexit("pair_cache_build received NULL cache");
    if (w == NULL)
        errorexit("pair_cache_build received NULL state vector");
    if (ode_params == NULL)
        errorexit("pair_cache_build received NULL ode_params");

    pair_cache_zero(cache);

    const int N = ode_params->num_bodies_initial;
    const int D = ode_params->num_dim;
    const int P0 = N * D;

    if (N < 0)
        errorexit("num_bodies_initial must be non-negative");
    if (D <= 0)
        errorexit("num_dim must be positive");

    cache->num_bodies = N;
    cache->num_dim = D;
    cache->array_half = P0;

    active_list_init(&cache->active, ode_params);

    cache->m       = (double *)xcalloc((size_t)N, sizeof(double));
    cache->inv_m   = (double *)xcalloc((size_t)N, sizeof(double));
    cache->p       = (double *)xcalloc((size_t)N * (size_t)D, sizeof(double));
    cache->p2      = (double *)xcalloc((size_t)N, sizeof(double));
    cache->p_dot   = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    cache->x_rel   = (double *)xcalloc((size_t)N * (size_t)N * (size_t)D, sizeof(double));
    cache->r       = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    cache->inv_r   = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    cache->inv_r2  = (double *)xcalloc((size_t)N * (size_t)N, sizeof(double));
    cache->n       = (double *)xcalloc((size_t)N * (size_t)N * (size_t)D, sizeof(double));
    cache->n_dot_p = (double *)xcalloc((size_t)N * (size_t)N * (size_t)N, sizeof(double));

    // Masses and momenta.
    for (int ia = 0; ia < cache->active.num_active; ++ia) {
        const int a = cache->active.ids[ia];

        cache->m[a] = ode_params->masses[a];
        if (cache->m[a] == 0.0)
            errorexit("Encountered zero active-body mass in pair_cache_build");
        cache->inv_m[a] = 1.0 / cache->m[a];

        double p2 = 0.0;
        for (int i = 0; i < D; ++i) {
            const double pai = w[P0 + a * D + i];
            cache->p[a * D + i] = pai;
            p2 += pai * pai;
        }
        cache->p2[a] = p2;
    }

    // Momentum scalar products p_a . p_b.
    for (int ia = 0; ia < cache->active.num_active; ++ia) {
        const int a = cache->active.ids[ia];
        for (int ib = ia; ib < cache->active.num_active; ++ib) {
            const int b = cache->active.ids[ib];

            double pdot = 0.0;
            for (int i = 0; i < D; ++i)
                pdot += cache->p[a * D + i] * cache->p[b * D + i];

            cache->p_dot[pair_cache_idx2(cache, a, b)] = pdot;
            cache->p_dot[pair_cache_idx2(cache, b, a)] = pdot;
        }
    }

    // Relative positions, distances, inverse distances and unit vectors.
    for (int ia = 0; ia < cache->active.num_active; ++ia) {
        const int a = cache->active.ids[ia];
        for (int ib = ia; ib < cache->active.num_active; ++ib) {
            const int b = cache->active.ids[ib];

            if (a == b)
                continue;

            double r2 = 0.0;
            for (int i = 0; i < D; ++i) {
                const double dx = w[a * D + i] - w[b * D + i];
                cache->x_rel[pair_cache_idx3(cache, a, b, i)] = dx;
                cache->x_rel[pair_cache_idx3(cache, b, a, i)] = -dx;
                r2 += dx * dx;
            }

            if (r2 <= 0.0)
                errorexit("Encountered zero separation in pair_cache_build");

            const double r = sqrt(r2);
            const double inv_r = 1.0 / r;
            const double inv_r2 = inv_r * inv_r;

            cache->r[pair_cache_idx2(cache, a, b)] = r;
            cache->r[pair_cache_idx2(cache, b, a)] = r;
            cache->inv_r[pair_cache_idx2(cache, a, b)] = inv_r;
            cache->inv_r[pair_cache_idx2(cache, b, a)] = inv_r;
            cache->inv_r2[pair_cache_idx2(cache, a, b)] = inv_r2;
            cache->inv_r2[pair_cache_idx2(cache, b, a)] = inv_r2;

            for (int i = 0; i < D; ++i) {
                const double ni = cache->x_rel[pair_cache_idx3(cache, a, b, i)] * inv_r;
                cache->n[pair_cache_idx3(cache, a, b, i)] = ni;
                cache->n[pair_cache_idx3(cache, b, a, i)] = -ni;
            }
        }
    }

    // n_ab . p_c for active triples. Entries with a == b remain zero.
    for (int ia = 0; ia < cache->active.num_active; ++ia) {
        const int a = cache->active.ids[ia];
        for (int ib = 0; ib < cache->active.num_active; ++ib) {
            const int b = cache->active.ids[ib];
            if (a == b)
                continue;

            for (int ic = 0; ic < cache->active.num_active; ++ic) {
                const int c = cache->active.ids[ic];

                double ndotp = 0.0;
                for (int i = 0; i < D; ++i)
                    ndotp += cache->n[pair_cache_idx3(cache, a, b, i)]
                           * cache->p[c * D + i];

                cache->n_dot_p[pair_cache_idx_ndotp(cache, a, b, c)] = ndotp;
            }
        }
    }
}


void pair_cache_free(PairCache *cache)
{
    if (cache == NULL)
        return;

    active_list_free(&cache->active);

    free(cache->m);
    free(cache->inv_m);
    free(cache->p);
    free(cache->p2);
    free(cache->p_dot);
    free(cache->x_rel);
    free(cache->r);
    free(cache->inv_r);
    free(cache->inv_r2);
    free(cache->n);
    free(cache->n_dot_p);

    pair_cache_zero(cache);
}
