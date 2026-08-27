/**
 * @file hamiltonian.c
 * @brief Functions for the computation of the post-Newtonian N-body Hamiltonian
 *
 * Functions for the computation of the post-Newtonian N-body Hamiltonian.
 * A complicated part of the N-body 2PN Hamiltonian is the four-point correlation function UTT4,
 * which contains a nontrivial four-point logarithmic integral. The logarithmic part is evaluated
 * with a two-dimensional parameter-space quadrature implemented in integrals.c.
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "eom.h"
#include "hamiltonian.h"
#include "utils.h"
#include "pair_cache.h"
#include "integrals.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NUM_LOCAL_BODIES 4


// ------------------------------------------------------------------------------------------------
// Cached UTT4 logarithmic integral
// ------------------------------------------------------------------------------------------------

// Convert the runtime parameter structure to the numerical-integral settings.
static NumericalIntegralSettings utt4_ln_integral_settings_from_ode(const struct ode_params *ode_params)
{
    NumericalIntegralSettings settings;
    settings.rel_tol = (long double)ode_params->utt4_ln_integral_epsrel;
    settings.abs_tol = (long double)ode_params->utt4_ln_integral_epsabs;
    settings.min_order = ode_params->utt4_ln_integral_min_order;
    settings.max_order = ode_params->utt4_ln_integral_max_order;
    settings.adaptive = ode_params->utt4_ln_integral_adaptive;
    settings.max_depth = ode_params->utt4_ln_integral_max_depth;
    settings.use_openmp = ode_params->utt4_ln_integral_parallel;
    return settings;
}


// Report a failed logarithmic-integral accuracy target only once during a run.
static void utt4_ln_integral_warning(const UTT4LnIntegralResult *result,
    const struct ode_params *ode_params)
{
    if (result->diagnostics.target_met)
        return;

    static int warned = 0;
    if (!warned) {
        progress_bar_break_line();
        printf("Warning: the UTT4 logarithmic integral did not reach the requested tolerance "
               "(epsrel=%.3e, last fixed orders %d/%d, worst tolerance ratio %.3Le). "
               "Consider increasing utt4_ln_integral_max_order or utt4_ln_integral_max_depth.\n",
               ode_params->utt4_ln_integral_epsrel, result->diagnostics.low_order,
               result->diagnostics.high_order,
               result->diagnostics.worst_tolerance_ratio);
        warned = 1;
    }
}


// Check whether the cached logarithmic UTT4 result is valid for the current position state.
static int utt4_ln_cache_matches(const PairCache *cache, const double *w,
    const struct ode_params *ode_params)
{
    const UTT4LnCache *utt4_cache = &cache->utt4_ln;
    if (!utt4_cache->valid)
        return 0;

    if (utt4_cache->epsrel != ode_params->utt4_ln_integral_epsrel
        || utt4_cache->epsabs != ode_params->utt4_ln_integral_epsabs
        || utt4_cache->min_order != ode_params->utt4_ln_integral_min_order
        || utt4_cache->max_order != ode_params->utt4_ln_integral_max_order
        || utt4_cache->adaptive != ode_params->utt4_ln_integral_adaptive
        || utt4_cache->max_depth != ode_params->utt4_ln_integral_max_depth
        || utt4_cache->parallel != ode_params->utt4_ln_integral_parallel)
        return 0;

    const int num_bodies = cache->num_bodies;
    const int num_dim = cache->num_dim;

    for (int a = 0; a < num_bodies; ++a) {
        if (utt4_cache->active[a] != ode_params->active[a])
            return 0;

        if (!ode_params->active[a])
            continue;

        if (utt4_cache->masses[a] != ode_params->masses[a])
            return 0;

        for (int axis = 0; axis < num_dim; ++axis) {
            const int idx = a * num_dim + axis;
            if (utt4_cache->positions[idx] != w[idx])
                return 0;
        }
    }

    return 1;
}


// Store the exact position/mass/active-set/settings key for a freshly evaluated UTT4 cache entry.
static void utt4_ln_cache_store_key(PairCache *cache, const double *w,
    const struct ode_params *ode_params)
{
    UTT4LnCache *utt4_cache = &cache->utt4_ln;
    const int num_bodies = cache->num_bodies;
    const int num_dim = cache->num_dim;

    for (int a = 0; a < num_bodies; ++a) {
        utt4_cache->active[a] = ode_params->active[a];
        utt4_cache->masses[a] = ode_params->masses[a];

        for (int axis = 0; axis < num_dim; ++axis) {
            const int idx = a * num_dim + axis;
            utt4_cache->positions[idx] = w[idx];
        }
    }

    utt4_cache->epsrel = ode_params->utt4_ln_integral_epsrel;
    utt4_cache->epsabs = ode_params->utt4_ln_integral_epsabs;
    utt4_cache->min_order = ode_params->utt4_ln_integral_min_order;
    utt4_cache->max_order = ode_params->utt4_ln_integral_max_order;
    utt4_cache->adaptive = ode_params->utt4_ln_integral_adaptive;
    utt4_cache->max_depth = ode_params->utt4_ln_integral_max_depth;
    utt4_cache->parallel = ode_params->utt4_ln_integral_parallel;
    utt4_cache->valid = 1;
}


// Number of unordered four-body subsets of an active list.
static size_t utt4_ln_num_quadruples(int num_active)
{
    if (num_active < NUM_LOCAL_BODIES)
        return 0;

    const size_t n = (size_t)num_active;
    return n*(n - 1)*(n - 2)*(n - 3)/24;
}


// Decode a lexicographic unordered-quadruple index without materializing a quadruple list.
static void utt4_ln_quadruple_from_index(const ActiveList *active, size_t index,
    int body[NUM_LOCAL_BODIES])
{
    const int n = active->num_active;
    size_t remaining = index;

    int ia = 0;
    for (; ia < n - 3; ++ia) {
        const size_t m = (size_t)(n - ia - 1);
        const size_t count = m*(m - 1)*(m - 2)/6;
        if (remaining < count)
            break;
        remaining -= count;
    }

    int ib = ia + 1;
    for (; ib < n - 2; ++ib) {
        const size_t m = (size_t)(n - ib - 1);
        const size_t count = m*(m - 1)/2;
        if (remaining < count)
            break;
        remaining -= count;
    }

    int ic = ib + 1;
    for (; ic < n - 1; ++ic) {
        const size_t count = (size_t)(n - ic - 1);
        if (remaining < count)
            break;
        remaining -= count;
    }

    const int id = ic + 1 + (int)remaining;
    body[0] = active->ids[ia];
    body[1] = active->ids[ib];
    body[2] = active->ids[ic];
    body[3] = active->ids[id];
}


// Evaluate one unordered quadruple. The caller decides whether the integral may use inner OpenMP.
static int utt4_ln_evaluate_quadruple(const double *w, const struct ode_params *ode_params,
    const int body[NUM_LOCAL_BODIES], const NumericalIntegralSettings *settings,
    UTT4LnIntegralResult *result)
{
    double pos[NUM_LOCAL_BODIES][3];
    for (int local = 0; local < NUM_LOCAL_BODIES; ++local)
        for (int axis = 0; axis < 3; ++axis)
            pos[local][axis] = w[ode_params->num_dim*body[local] + axis];

    return utt4_ln_integral_evaluate(pos, settings, result);
}


// Accumulate one mass-weighted quadruple result into a scalar and system gradient.
static void utt4_ln_accumulate_quadruple(const struct ode_params *ode_params,
    const int body[NUM_LOCAL_BODIES], const UTT4LnIntegralResult *result,
    long double *value, long double *grad)
{
    const long double mass_fac =
        (long double)ode_params->masses[body[0]] * ode_params->masses[body[1]]
      * (long double)ode_params->masses[body[2]] * ode_params->masses[body[3]];
    const long double prefactor = 0.25L*mass_fac;

    *value += prefactor*result->value;

    for (int local = 0; local < NUM_LOCAL_BODIES; ++local) {
        const int global = body[local];
        for (int axis = 0; axis < 3; ++axis) {
            const int idx = ode_params->num_dim*global + axis;
            grad[idx] += prefactor*result->grad[local][axis];
        }
    }
}


// Evaluate the complete mass-weighted logarithmic UTT4 value and gradient for one position state.
static void utt4_ln_cache_refresh(PairCache *cache, double *w, struct ode_params *ode_params)
{
    const int num_dim = ode_params->num_dim;
    const int array_half = cache->array_half;
    const NumericalIntegralSettings settings = utt4_ln_integral_settings_from_ode(ode_params);
    const ActiveList *active = &cache->active;
    UTT4LnCache *utt4_cache = &cache->utt4_ln;
    const size_t num_quadruples = utt4_ln_num_quadruples(active->num_active);

    if (num_dim != 3)
        errorexit("The UTT4 logarithmic integral can only be computed in 3D! Use num_dim = 3");

    for (int i = 0; i < array_half; ++i)
        utt4_cache->grad[i] = 0.0;

    long double sum = 0.0L;

#ifdef _OPENMP
    /*
     * For N > 4, parallelize over independent unordered quadruples rather than inside each
     * two-dimensional quadrature. This gives coarser tasks, avoids repeated nested parallel
     * regions, and scales better as C(N,4) grows. N = 4 retains the inner quadrature parallelism.
     */
    const int use_quadruple_parallelism =
        settings.use_openmp && num_quadruples > 1 && omp_get_max_threads() > 1
        && !omp_in_parallel();

    if (use_quadruple_parallelism) {
        const int max_threads = omp_get_max_threads();
        long double *thread_sum = calloc((size_t)max_threads, sizeof(*thread_sum));
        long double *thread_grad = calloc((size_t)max_threads*(size_t)array_half,
            sizeof(*thread_grad));
        int *thread_failed = calloc((size_t)max_threads, sizeof(*thread_failed));
        int *thread_warn = calloc((size_t)max_threads, sizeof(*thread_warn));
        UTT4LnIntegralResult *warning_result =
            calloc((size_t)max_threads, sizeof(*warning_result));

        if (!thread_sum || !thread_grad || !thread_failed || !thread_warn || !warning_result)
            errorexit("Failed to allocate UTT4 quadruple-parallel workspace");

        NumericalIntegralSettings serial_settings = settings;
        serial_settings.use_openmp = 0;

        #pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            long double *local_grad = thread_grad + (size_t)tid*(size_t)array_half;

            #pragma omp for schedule(static)
            for (size_t q = 0; q < num_quadruples; ++q) {
                int body[NUM_LOCAL_BODIES];
                UTT4LnIntegralResult result;
                utt4_ln_quadruple_from_index(active, q, body);

                if (utt4_ln_evaluate_quadruple(w, ode_params, body, &serial_settings, &result) != 0) {
                    thread_failed[tid] = 1;
                    continue;
                }

                if (!result.diagnostics.target_met && !thread_warn[tid]) {
                    warning_result[tid] = result;
                    thread_warn[tid] = 1;
                }

                utt4_ln_accumulate_quadruple(ode_params, body, &result,
                    &thread_sum[tid], local_grad);
            }
        }

        int failed = 0;
        int warned = 0;
        for (int tid = 0; tid < max_threads; ++tid) {
            failed |= thread_failed[tid];
            sum += thread_sum[tid];

            if (!warned && thread_warn[tid]) {
                utt4_ln_integral_warning(&warning_result[tid], ode_params);
                warned = 1;
            }
        }

        for (int i = 0; i < array_half; ++i) {
            long double component = 0.0L;
            for (int tid = 0; tid < max_threads; ++tid)
                component += thread_grad[(size_t)tid*(size_t)array_half + i];
            utt4_cache->grad[i] = (double)component;
        }

        free(thread_sum);
        free(thread_grad);
        free(thread_failed);
        free(thread_warn);
        free(warning_result);

        if (failed)
            errorexit("UTT4 logarithmic integral evaluation failed");
    } else
#endif
    {
        long double serial_grad[cache->array_half];
        for (int i = 0; i < array_half; ++i)
            serial_grad[i] = 0.0L;

        for (size_t q = 0; q < num_quadruples; ++q) {
            int body[NUM_LOCAL_BODIES];
            UTT4LnIntegralResult result;
            utt4_ln_quadruple_from_index(active, q, body);

            if (utt4_ln_evaluate_quadruple(w, ode_params, body, &settings, &result) != 0)
                errorexit("UTT4 logarithmic integral evaluation failed");
            utt4_ln_integral_warning(&result, ode_params);
            utt4_ln_accumulate_quadruple(ode_params, body, &result, &sum, serial_grad);
        }

        for (int i = 0; i < array_half; ++i)
            utt4_cache->grad[i] = (double)serial_grad[i];
    }

    utt4_cache->value = (double)sum;
    utt4_ln_cache_store_key(cache, w, ode_params);
}


/**
 * @brief Return the cached logarithmic UTT4 value and/or gradient for the current positions.
 *
 * The expensive numerical quadrature is performed only when the active-body positions, masses,
 * active set, or logarithmic-integral settings differ from the cached state. Either output may be
 * NULL. The returned value and gradient already include the four-mass products and the Hamiltonian
 * symmetry factor 1/4.
 */
void utt4_ln_integral_cached(double *w, struct ode_params *ode_params, double *value, double *grad)
{
    if (w == NULL || ode_params == NULL)
        errorexit("utt4_ln_integral_cached received a NULL input");

    PairCache *cache = pair_cache_get_workspace(ode_params);
    active_list_refresh(&cache->active, ode_params);

    if (!utt4_ln_cache_matches(cache, w, ode_params))
        utt4_ln_cache_refresh(cache, w, ode_params);

    if (value != NULL)
        *value = cache->utt4_ln.value;

    if (grad != NULL) {
        for (int i = 0; i < cache->array_half; ++i)
            grad[i] = cache->utt4_ln.grad[i];
    }
}


// ------------------------------------------------------------------------------------------------
// N-Body post-Newtonian Hamiltonians
// ------------------------------------------------------------------------------------------------


// Computes the 0PN (Newtonian) Hamiltonian part from an already-refreshed cache.
double H0PN_cached(const PairCache *cache)
{
    double H = 0.0;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        int a = cache->active.ids[ia];

        // Kinetic energy
        H += cache->p2[a] / (2.0 * cache->m[a]);

        // Potential energy
        for (int ib = ia + 1; ib < cache->active.num_active; ib++) {
            int b = cache->active.ids[ib];
            H -= cache->m[a] * cache->m[b] * pair_cache_inv_r(cache, a, b);
        }
    }

    return H;
}


// Computes the 0PN (Newtonian) Hamiltonian part cleanly with refreshing the cache.
double H0PN(double* w, struct ode_params* ode_params)
{
    PairCache *cache = pair_cache_get_workspace(ode_params);
    const unsigned int levels = PAIR_CACHE_LEVEL_GEOMETRY | PAIR_CACHE_LEVEL_P2;
    pair_cache_refresh(cache, w, ode_params, levels);
    return H0PN_cached(cache);
}


// Computes the 1PN Hamiltonian part from an already-refreshed cache.
double H1PN_cached(const PairCache *cache)
{
    double H = 0.0;

    // Compute kinetic and potential energy
    for (int ia = 0; ia < cache->active.num_active; ia++) {
        int a = cache->active.ids[ia];

        double m_a = cache->m[a];
        double pa_dot_pa = cache->p2[a];

        H += -0.125 * m_a * pa_dot_pa * pa_dot_pa / (m_a * m_a * m_a * m_a);

        for (int ib = 0; ib < cache->active.num_active; ib++) {
            int b = cache->active.ids[ib];
            if (b == a) continue;

            double m_b = cache->m[b];
            double r_ab = pair_cache_r(cache, a, b);
            double pa_dot_pb = pair_cache_p_dot(cache, a, b);
            double na_dot_pa = pair_cache_n_dot_p(cache, a, b, a);
            double na_dot_pb = pair_cache_n_dot_p(cache, a, b, b);

            H += -0.25 * m_a * m_b / r_ab * (6.0 * pa_dot_pa / (m_a * m_a)
                - 7.0 * pa_dot_pb / (m_a * m_b) - (na_dot_pa * na_dot_pb) / (m_a * m_b));

            for (int ic = 0; ic < cache->active.num_active; ic++) {
                int c = cache->active.ids[ic];
                if (c == a) continue;

                double m_c = cache->m[c];
                double r_ac = pair_cache_r(cache, a, c);

                H += 0.5 * m_a * m_b * m_c / (r_ab * r_ac);
            }
        }
    }

    return H;
}


// Computes the 1PN Hamiltonian part cleanly with refreshing the cache.
double H1PN(double* w, struct ode_params* ode_params)
{
    PairCache *cache = pair_cache_get_workspace(ode_params);
    const unsigned int levels = PAIR_CACHE_LEVEL_GEOMETRY | PAIR_CACHE_LEVEL_P2
        | PAIR_CACHE_LEVEL_PAIR_DOTS;
    pair_cache_refresh(cache, w, ode_params, levels);
    return H1PN_cached(cache);
}


// Computes the 2PN Hamiltonian part from an already-refreshed cache.
double H2PN_cached(double* w, struct ode_params* ode_params, const PairCache *cache, int utt4_flag)
{
    int num_bodies = ode_params->num_bodies_initial;
    int num_dim = ode_params->num_dim;
    const ActiveList *active = &cache->active;
    double temp, temp0, temp1, temp2, temp3, temp4, temp5, temp6,
        temp7, temp8, temp9, temp10, temp11, temp12, temp13;
    double ma_inv, mb_inv, mc_inv;

    const double *m = cache->m;
    const double (*p)[num_dim] = (const double (*)[num_dim])cache->p;
    const double (*n)[num_bodies][num_dim] = (const double (*)[num_bodies][num_dim])cache->n;
    const double (*r)[num_bodies] = (const double (*)[num_bodies])cache->r;

    double n_ab_ac[num_dim];
    double n_ab_cb[num_dim];

    // Compute H
    double H = 0.0;
    for (int ia = 0; ia < active->num_active; ++ia) {
        const int a = active->ids[ia];
        ma_inv = 1/m[a];
        temp = cache->p2[a];
        temp2 = temp * ma_inv * ma_inv;

        H += 0.0625 * m[a] * temp2 * temp2 * temp2;

        for (int ib = 0; ib < active->num_active; ++ib) {
            const int b = active->ids[ib];
            mb_inv = 1/m[b];
            temp0 = r[a][b] * r[a][b];
            temp1 = m[a] * m[b];
            temp3 = temp * ma_inv * mb_inv;
            temp4 = pair_cache_n_dot_p(cache, a, b, b);
            temp5 = pair_cache_n_dot_p(cache, a, b, a);
            temp6 = cache->p2[b];
            temp7 = pair_cache_p_dot(cache, a, b) * ma_inv * mb_inv;

            if (b != a){
                H += 0.0625 * 1 / r[a][b] * (10 * temp1 * temp2 * temp2
                    - 11 * temp3 * temp6
                    - 2 * pair_cache_p_dot(cache, a, b) * temp7
                    + 10 * temp3 * temp4 * temp4
                    - 12 * temp7 * temp5 * temp4
                    - 3 * temp5 * temp5 * temp4 * temp4 * ma_inv * mb_inv);
                H += 0.25 * m[a] * temp1 / temp0 * (temp2
                    + temp6 * mb_inv * mb_inv
                    - 2 * temp7);
                H += -0.25 * temp1 * temp1 / (temp0 * r[a][b]);
            }
            for (int ic = 0; ic < active->num_active; ++ic) {
                const int c = active->ids[ic];
                mc_inv = 1/m[c];
                temp8 = pair_cache_n_dot_p(cache, a, c, c);
                temp9 = dot_product(n[a][b], n[a][c], num_dim);
                temp10 = pair_cache_n_dot_p(cache, a, b, c);
                temp11 = r[b][c] * r[b][c];

                if (b != a && c != a) {
                    H += 0.125 * temp1 * m[c] / (r[a][b] * r[a][c]) * (18 * temp2
                        + 14 * temp6 * mb_inv * mb_inv
                        - 2 * temp4 * temp4 * mb_inv * mb_inv
                        - 50 * temp7
                        + 17 * pair_cache_p_dot(cache, b, c) * mb_inv * mc_inv
                        - 14 * temp5 * temp4 * ma_inv * mb_inv
                        + 14 * temp4 * temp10 * mb_inv * mc_inv
                        + temp9 * temp4 * temp8 * mb_inv * mc_inv);
                    H += 0.125 * temp1 * m[c] / (r[a][b] * r[a][b]) * (
                        2 * temp5 * temp8 * ma_inv * mc_inv
                        + 2 * temp4 * temp8 * ma_inv * mc_inv
                        + 5 * temp9 * cache->p2[c] * mc_inv * mc_inv
                        - temp9 * temp8 * temp8 * mc_inv * mc_inv
                        - 14 * temp10 * temp8 * mc_inv * mc_inv
                        );
                }
                if (b != a && c != a && c != b) {
                    for (int i = 0; i < num_dim; i++) {
                        n_ab_ac[i] = n[a][b][i] + n[a][c][i];
                        n_ab_cb[i] = n[a][b][i] + n[c][b][i];
                    }
                    H += 0.5 * temp1 * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * (r[a][b] + r[b][c] + r[c][a])) * (
                        8 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[c], num_dim) * ma_inv * mc_inv
                        - 16 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[a], num_dim) * ma_inv * mc_inv
                        + 3 * dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[b], num_dim) * ma_inv * mb_inv
                        + 4 * dot_product(n_ab_ac, p[c], num_dim) * dot_product(n_ab_cb, p[c], num_dim) * mc_inv * mc_inv
                        + dot_product(n_ab_ac, p[a], num_dim) * dot_product(n_ab_cb, p[a], num_dim) * ma_inv * ma_inv);
                    H += 0.5 * temp1 * m[c] / ((r[a][b] + r[b][c] + r[c][a]) * r[a][b]) * (
                        8 * (pair_cache_p_dot(cache, a, c) - temp5 * temp10) * ma_inv * mc_inv
                        - 3 * (temp7 - temp5 * temp4 * ma_inv * mb_inv)
                        - 4 * (cache->p2[c] - temp10 * temp10) * mc_inv * mc_inv
                        - (temp - temp5 * temp5) * ma_inv * ma_inv);
                    H += -0.015625 * m[a] * temp1 * m[c] / (temp0 * r[a][b] * r[a][c] * r[a][c] * r[a][c] * r[b][c]) * (
                        18 * temp0 * r[a][c] * r[a][c] - 60 * temp0 * temp11
                        - 24 * temp0 * r[a][c] * (r[a][b] + r[b][c])
                        + 60 * r[a][b] * r[a][c] * temp11 + 56 * temp0 * r[a][b] * r[b][c]
                        - 72 * r[a][b] * r[b][c] * temp11 + 35 * temp11 * temp11 + 6 * temp0 * temp0);
                }
                // G^3 quadruple sum terms
                for (int id = 0; id < active->num_active; ++id) {
                    const int d = active->ids[id];
                    temp12 = r[c][d]*r[c][d];
                    temp13 = r[a][d]*r[a][d];

                    if (b != a && c != b && d != c) {
                        H += - 0.375 * temp1 * m[c] * m[d] / (r[a][b] * r[b][c] * r[c][d]);
                    }
                    if (b != a && c != a && d != a) {
                        H += - 0.25 * temp1 * m[c] * m[d] / (r[a][b] * r[a][c] * r[a][d]);
                    }
                    if (utt4_flag) {
                        if (b != a && c != a && c != b && d != a && d != b && d != c) {
                            H += - 0.015625 * m[a] * m[b] * m[c] * m[d] / (temp0*r[a][b] * temp12*r[c][d] * temp13*r[a][d] * temp11*r[b][c]) * (
                                16 * temp0*r[a][b] * temp11*r[b][c] * temp12 * temp13 / r[b][d]
                                - 24 * temp11*r[b][c] * temp13 * temp0 * temp12
                                - 30 * temp13 * temp13 * temp11*r[b][c] * (temp13 + temp11 - r[a][c]*r[a][c] - r[b][d]*r[b][d])
                                + temp0 * (r[b][d]*r[b][d] - temp11 - temp12) * (
                                    -8 * temp13 * r[a][d] * temp11 + 16 * r[a][b] * temp13*r[a][d] * temp11 / (r[a][c] + r[b][c] + r[a][b])
                                        + r[a][b] * temp12 * (r[a][c]*r[a][c] - temp13 - temp12)) );
                        }
                    }
                }
            }
        }
    }
    // Logarithmic integral part of UTT4
    if (utt4_flag) {
        double utt4_ln_value;
        utt4_ln_integral_cached(w, ode_params, &utt4_ln_value, NULL);
        H += utt4_ln_value;
    }
    return H;
}


// Computes the 2PN Hamiltonian part cleanly with refreshing the cache.
double H2PN(double* w, struct ode_params* ode_params, int utt4_flag)
{
    PairCache *cache = pair_cache_get_workspace(ode_params);
    const unsigned int levels = PAIR_CACHE_LEVEL_GEOMETRY | PAIR_CACHE_LEVEL_P2 | PAIR_CACHE_LEVEL_PAIR_DOTS;
    pair_cache_refresh(cache, w, ode_params, levels);
    return H2PN_cached(w, ode_params, cache, utt4_flag);
}
