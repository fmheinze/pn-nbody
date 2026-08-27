/**
 * @file integrals.c
 * @brief Numerical integral implementations used by the N-body Hamiltonian.
 *
 * Generic numerical quadrature infrastructure as well as problem-specific kernels, e.g.
 * utt4_ln_integral_evaluate(), which computes the differentiated logarithmic four-body integral
 * entering UTT4 together with all four Cartesian position gradients.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "integrals.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif


// ------------------------------------------------------------------------------------------------
// Generic quadrature machinery
// ------------------------------------------------------------------------------------------------

typedef struct {
    int n;
    long double *x;
    long double *weight;
} GaussLegendreRule;


// Build an n-point Gauss-Legendre rule on [a, b].
static int gauss_legendre_rule_init(int n, long double a, long double b, GaussLegendreRule *rule)
{
    rule->n = n;
    rule->x = calloc((size_t)n, sizeof(long double));
    rule->weight = calloc((size_t)n, sizeof(long double));
    if (!rule->x || !rule->weight)
        return -1;

    const long double eps = 8.0L * LDBL_EPSILON;
    const long double midpoint = 0.5L * (a + b);
    const long double half_width = 0.5L * (b - a);
    const int m = (n + 1) / 2;

    for (int i = 0; i < m; ++i) {
        long double z = cosl(M_PI * ((long double)i + 0.75L) / ((long double)n + 0.5L));
        long double z1;
        long double pp = 0.0L;
        do {
            long double p1 = 1.0L;
            long double p2 = 0.0L;
            for (int j = 1; j <= n; ++j) {
                const long double p3 = p2;
                p2 = p1;
                p1 = (((2.0L * j - 1.0L) * z * p2) - (j - 1.0L) * p3) / (long double)j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0L);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabsl(z - z1) > eps);

        const long double w = 2.0L / ((1.0L - z * z) * pp * pp);
        const int lo = i;
        const int hi = n - 1 - i;
        rule->x[lo] = midpoint - half_width * z;
        rule->x[hi] = midpoint + half_width * z;
        rule->weight[lo] = half_width * w;
        rule->weight[hi] = half_width * w;
    }
    return 0;
}


static void gauss_legendre_rule_free(GaussLegendreRule *rule)
{
    free(rule->x);
    free(rule->weight);
    memset(rule, 0, sizeof(*rule));
}


// ------------------------------------------------------------------------------------------------
// UTT4 logarithmic integral
// ------------------------------------------------------------------------------------------------

#define UTT4_LN_NPAIR 6

// Six representatives for the 24 ordered permutations.
static const int utt4_ln_sum_representatives[6][4] =
{
    {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1}
};

// UTT4-specific precomputed transforms of a generic [0, pi] Gauss-Legendre rule.
typedef struct {
    GaussLegendreRule base;
    long double *half_sin_theta;
    long double *mix;
} UTT4LnQuadratureRule;


static int utt4_ln_quadrature_rule_init(int n, UTT4LnQuadratureRule *rule)
{
    memset(rule, 0, sizeof(*rule));
    if (gauss_legendre_rule_init(n, 0.0L, M_PI, &rule->base) != 0)
        return -1;

    rule->half_sin_theta = calloc((size_t)n, sizeof(long double));
    rule->mix = calloc((size_t)n, sizeof(long double));
    if (!rule->half_sin_theta || !rule->mix) {
        gauss_legendre_rule_free(&rule->base);
        free(rule->half_sin_theta);
        free(rule->mix);
        memset(rule, 0, sizeof(*rule));
        return -1;
    }

    for (int i = 0; i < n; ++i) {
        rule->half_sin_theta[i] = 0.5L*sinl(rule->base.x[i]);
        rule->mix[i] = 0.5L*(1.0L + cosl(rule->base.x[i]));
    }
    return 0;
}


/**
 * @brief Process-lifetime cache of immutable Gauss-Legendre rules.
 *
 * The same few quadrature orders are reused at virtually every UTT4 evaluation. Keeping them
 * avoids repeated root finding, trigonometric transforms, and heap allocation. Entries are never
 * modified after construction and intentionally live until process exit.
 */
typedef struct UTT4LnRuleCacheEntry {
    UTT4LnQuadratureRule rule;
    struct UTT4LnRuleCacheEntry *next;
} UTT4LnRuleCacheEntry;

static UTT4LnRuleCacheEntry *utt4_ln_rule_cache = NULL;
static int utt4_ln_rule_cache_cleanup_registered = 0;


static void utt4_ln_quadrature_rule_destroy(UTT4LnQuadratureRule *rule)
{
    gauss_legendre_rule_free(&rule->base);
    free(rule->half_sin_theta);
    free(rule->mix);
    memset(rule, 0, sizeof(*rule));
}


static void utt4_ln_rule_cache_free(void)
{
    UTT4LnRuleCacheEntry *entry = utt4_ln_rule_cache;
    while (entry != NULL) {
        UTT4LnRuleCacheEntry *next = entry->next;
        utt4_ln_quadrature_rule_destroy(&entry->rule);
        free(entry);
        entry = next;
    }
    utt4_ln_rule_cache = NULL;
}


static const UTT4LnQuadratureRule *utt4_ln_quadrature_rule_get(int n)
{
    const UTT4LnQuadratureRule *result = NULL;

#ifdef _OPENMP
    #pragma omp critical(utt4_ln_rule_cache_guard)
#endif
    {
        UTT4LnRuleCacheEntry *entry = utt4_ln_rule_cache;
        while (entry != NULL && entry->rule.base.n != n)
            entry = entry->next;

        if (entry == NULL) {
            entry = malloc(sizeof(*entry));
            if (entry != NULL) {
                memset(entry, 0, sizeof(*entry));
                if (utt4_ln_quadrature_rule_init(n, &entry->rule) != 0) {
                    free(entry);
                    entry = NULL;
                } else {
                    entry->next = utt4_ln_rule_cache;
                    utt4_ln_rule_cache = entry;
                    if (!utt4_ln_rule_cache_cleanup_registered) {
                        if (atexit(utt4_ln_rule_cache_free) == 0)
                            utt4_ln_rule_cache_cleanup_registered = 1;
                    }
                }
            }
        }

        if (entry != NULL)
            result = &entry->rule;
    }

    return result;
}


// Global pair ordering: 01, 02, 03, 12, 13, 23
static const int utt4_ln_pair_body_i[UTT4_LN_NPAIR] = {0,0,0,1,1,2};
static const int utt4_ln_pair_body_j[UTT4_LN_NPAIR] = {1,2,3,2,3,3};


static int utt4_ln_pair_index(int a, int b)
{
    if (a > b) { const int t=a; a=b; b=t; }
    if (a==0 && b==1) return 0;
    if (a==0 && b==2) return 1;
    if (a==0 && b==3) return 2;
    if (a==1 && b==2) return 3;
    if (a==1 && b==3) return 4;
    if (a==2 && b==3) return 5;
    return -1;
}


typedef struct {
    long double v;
    long double g[UTT4_LN_NPAIR];
} Dual6;


static inline Dual6 d6_zero(void)
{
    Dual6 z;
    z.v = 0.0L;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        z.g[k] = 0.0L;
    return z;
}


static inline Dual6 d6_const(long double x)
{
    Dual6 z = d6_zero();
    z.v = x;
    return z;
}


static inline Dual6 d6_seed(long double x, int k)
{
    Dual6 z = d6_const(x);
    z.g[k] = 1.0L;
    return z;
}


static inline Dual6 d6_add(Dual6 a, Dual6 b)
{
    Dual6 r;
    r.v = a.v + b.v;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        r.g[k] = a.g[k] + b.g[k];
    return r;
}


static inline Dual6 d6_sub(Dual6 a, Dual6 b)
{
    Dual6 r;
    r.v = a.v - b.v;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        r.g[k] = a.g[k] - b.g[k];
    return r;
}


static inline Dual6 d6_scale(Dual6 a, long double c)
{
    Dual6 r;
    r.v = a.v * c;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        r.g[k] = a.g[k] * c;
    return r;
}


static inline Dual6 d6_mul(Dual6 a, Dual6 b)
{
    Dual6 r;
    r.v = a.v * b.v;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        r.g[k] = a.g[k] * b.v + a.v * b.g[k];
    return r;
}


static inline Dual6 d6_inv(Dual6 a)
{
    const long double iv = 1.0L / a.v;
    Dual6 r;
    r.v = iv;
    const long double fac = -iv*iv;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        r.g[k] = fac * a.g[k];
    return r;
}


static inline Dual6 d6_sqrt(Dual6 a)
{
    const long double sv = sqrtl(a.v);
    Dual6 r;
    r.v = sv;
    const long double fac = 0.5L / sv;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        r.g[k] = fac * a.g[k];
    return r;
}



#define UTT4_LN_NUM_GENERATOR_DERIVATIVES 21
#define UTT4_LN_D_INDEX(i, j) (6*(j) - ((j)*((j) - 1))/2 + (i))
#define UTT4_LN_D(D, i, j) ((D)[UTT4_LN_D_INDEX((i), (j))])


// ------------------------------------------------------------------------------------------------
// Reduced generator derivatives
// ------------------------------------------------------------------------------------------------

/**
 * @brief Evaluate H'(y), ..., H^(5)(y).
 *
 * H itself is not required by the contracted UTT4 kernel. For small y, the five derivatives are
 * accumulated simultaneously from the convergent power series. Away from y = 0, the auxiliary
 *
 *     a(y) = atan(sqrt(y))/sqrt(y)
 *
 * is differentiated with
 *
 *     2 y a^(n+1) + (2 n + 1) a^n = (-1)^n n! / (1 + y)^(n+1).
 *
 * Avoiding H(y) itself also removes an unnecessary log1p() from every quadrature node.
 */
static void utt4_ln_H_derivatives(long double y, long double h[6])
{
    h[0] = 0.0L;

    if (y < 0.10L) {
        for (int k = 1; k <= 5; ++k)
            h[k] = 0.0L;

        long double power[6] = {0.0L, 1.0L, 0.0L, 0.0L, 0.0L, 0.0L};
        for (int n = 1; n <= 128; ++n) {
            if (n <= 5 && n >= 2)
                power[n] = 1.0L;

            const long double nn = (long double)n;
            long double cn = 1.0L/(nn*(4.0L*nn*nn - 1.0L));
            if (n & 1)
                cn = -cn;

            long double falling = nn;
            int converged = (n > 17);
            const int kmax = (n < 5) ? n : 5;
            for (int k = 1; k <= kmax; ++k) {
                const long double term = cn*falling*power[k];
                h[k] += term;
                if (converged &&
                    fabsl(term) >= 2.0L*LDBL_EPSILON*fmaxl(1.0L, fabsl(h[k])))
                    converged = 0;
                falling *= (long double)(n - k);
            }

            if (converged)
                break;

            for (int k = 1; k <= kmax; ++k)
                power[k] *= y;
        }
        return;
    }

    const long double inv_one_plus_y = 1.0L/(1.0L + y);
    const long double two_y = 2.0L*y;
    const long double t = sqrtl(y);
    long double a[6];
    a[0] = atanl(t)/t;

    long double rhs_scale = inv_one_plus_y;
    long double factorial = 1.0L;
    long double sign = 1.0L;
    for (int n = 0; n < 5; ++n) {
        a[n + 1] = (sign*factorial*rhs_scale - (2.0L*n + 1.0L)*a[n])/two_y;
        sign = -sign;
        factorial *= (long double)(n + 1);
        rhs_scale *= inv_one_plus_y;
    }

    long double log_derivative = inv_one_plus_y;
    for (int n = 1; n <= 5; ++n) {
        h[n] = (1.0L - y)*a[n] - (long double)n*a[n - 1] + log_derivative;
        log_derivative *= -(long double)n*inv_one_plus_y;
    }
}


/**
 * @brief Compute the generator derivatives needed by the contracted UTT4 kernel.
 *
 * Only 18 of the 21 derivatives with i+j <= 5 are ever used. Their recurrence can be written
 * explicitly in terms of y=q/s^2 and H'(y),...,H^(5)(y), which removes the generic coefficient
 * recurrence, three unused derivatives, and the logarithm that appears only in D00 and D10.
 */
static void utt4_ln_generator_derivatives(long double s, long double q,
    long double D[UTT4_LN_NUM_GENERATOR_DERIVATIVES])
{
    const long double inv_s = 1.0L/s;
    const long double inv_s2 = inv_s*inv_s;
    const long double inv_s3 = inv_s2*inv_s;
    const long double inv_s4 = inv_s2*inv_s2;
    const long double inv_s5 = inv_s4*inv_s;
    const long double inv_s6 = inv_s3*inv_s3;
    const long double inv_s7 = inv_s6*inv_s;
    const long double inv_s8 = inv_s4*inv_s4;
    const long double inv_s9 = inv_s8*inv_s;

    const long double y = q*inv_s2;
    const long double y2 = y*y;
    const long double y3 = y2*y;
    const long double y4 = y2*y2;

    long double h[6];
    utt4_ln_H_derivatives(y, h);
    const long double h1 = h[1];
    const long double h2 = h[2];
    const long double h3 = h[3];
    const long double h4 = h[4];
    const long double h5 = h[5];

    UTT4_LN_D(D, 0, 2) = h2*inv_s3;
    UTT4_LN_D(D, 0, 3) = h3*inv_s5;
    UTT4_LN_D(D, 0, 4) = h4*inv_s7;
    UTT4_LN_D(D, 0, 5) = h5*inv_s9;

    UTT4_LN_D(D, 1, 1) = (-h1 - 2.0L*h2*y)*inv_s2;
    UTT4_LN_D(D, 1, 2) = (-3.0L*h2 - 2.0L*h3*y)*inv_s4;
    UTT4_LN_D(D, 1, 3) = (-5.0L*h3 - 2.0L*h4*y)*inv_s6;
    UTT4_LN_D(D, 1, 4) = (-7.0L*h4 - 2.0L*h5*y)*inv_s8;

    UTT4_LN_D(D, 2, 0) = (2.0L*y*(h1 + 2.0L*h2*y) + 2.0L)*inv_s;
    UTT4_LN_D(D, 2, 1) =
        2.0L*(h1 + 5.0L*h2*y + 2.0L*h3*y2)*inv_s3;
    UTT4_LN_D(D, 2, 2) =
        2.0L*(6.0L*h2 + 9.0L*h3*y + 2.0L*h4*y2)*inv_s5;
    UTT4_LN_D(D, 2, 3) =
        2.0L*(15.0L*h3 + 13.0L*h4*y + 2.0L*h5*y2)*inv_s7;

    UTT4_LN_D(D, 3, 0) =
        (-2.0L*y*(3.0L*h1 + 12.0L*h2*y + 4.0L*h3*y2) - 2.0L)*inv_s2;
    UTT4_LN_D(D, 3, 1) =
        -2.0L*(3.0L*h1 + 27.0L*h2*y + 24.0L*h3*y2 + 4.0L*h4*y3)*inv_s4;
    UTT4_LN_D(D, 3, 2) =
        -2.0L*(30.0L*h2 + 75.0L*h3*y + 36.0L*h4*y2 + 4.0L*h5*y3)*inv_s6;

    UTT4_LN_D(D, 4, 0) =
        (4.0L*y*(6.0L*h1 + 39.0L*h2*y + 28.0L*h3*y2 + 4.0L*h4*y3)
         + 4.0L)*inv_s3;
    UTT4_LN_D(D, 4, 1) =
        4.0L*(6.0L*h1 + 84.0L*h2*y + 123.0L*h3*y2 + 44.0L*h4*y3
              + 4.0L*h5*y4)*inv_s5;

    UTT4_LN_D(D, 5, 0) =
        (-4.0L*y*(30.0L*h1 + 285.0L*h2*y + 330.0L*h3*y2 + 100.0L*h4*y3
                   + 8.0L*h5*y4) - 12.0L)*inv_s4;
}


// ------------------------------------------------------------------------------------------------
// Pair geometry and distance derivatives
// ------------------------------------------------------------------------------------------------

typedef struct {
    int pab, pcd, pac, pad, pbc, pbd;
    Dual6 r;
    Dual6 rho;
    Dual6 inv_r;
    Dual6 inv_rho;
    Dual6 inv_r_rho;
    Dual6 xy;
    Dual6 xw;
    Dual6 wy;
    Dual6 nm;
} UTT4LnRepresentativeGeometry;


typedef struct {
    long double d2[UTT4_LN_NPAIR];
    Dual6 d2dual[UTT4_LN_NPAIR];
    Dual6 distance[UTT4_LN_NPAIR];
    Dual6 inv_distance[UTT4_LN_NPAIR];
    UTT4LnRepresentativeGeometry representative[6];
} UTT4LnPairGeometry;


static void utt4_ln_build_pair_geometry(const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3],
    UTT4LnPairGeometry *geo)
{
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const int a = utt4_ln_pair_body_i[p], b = utt4_ln_pair_body_j[p];
        long double d2 = 0.0L;
        for (int k = 0; k < 3; ++k) {
            const long double dx = (long double)pos[a][k] - (long double)pos[b][k];
            d2 += dx*dx;
        }
        geo->d2[p] = d2;
        geo->d2dual[p] = d6_seed(d2, p);
        geo->distance[p] = d6_sqrt(geo->d2dual[p]);
        geo->inv_distance[p] = d6_inv(geo->distance[p]);
    }

    for (int p = 0; p < 6; ++p) {
        const int *body = utt4_ln_sum_representatives[p];
        const int a = body[0], b = body[1], c = body[2], d = body[3];
        UTT4LnRepresentativeGeometry *rep = &geo->representative[p];

        rep->pab = utt4_ln_pair_index(a, b);
        rep->pcd = utt4_ln_pair_index(c, d);
        rep->pac = utt4_ln_pair_index(a, c);
        rep->pad = utt4_ln_pair_index(a, d);
        rep->pbc = utt4_ln_pair_index(b, c);
        rep->pbd = utt4_ln_pair_index(b, d);

        const Dual6 dab = geo->d2dual[rep->pab];
        const Dual6 dcd = geo->d2dual[rep->pcd];
        const Dual6 dac = geo->d2dual[rep->pac];
        const Dual6 dad = geo->d2dual[rep->pad];
        const Dual6 dbc = geo->d2dual[rep->pbc];
        const Dual6 dbd = geo->d2dual[rep->pbd];

        rep->r = geo->distance[rep->pab];
        rep->rho = geo->distance[rep->pcd];
        rep->inv_r = geo->inv_distance[rep->pab];
        rep->inv_rho = geo->inv_distance[rep->pcd];
        rep->inv_r_rho = d6_mul(rep->inv_r, rep->inv_rho);
        rep->xy = d6_scale(d6_add(d6_sub(dad, dac), d6_sub(dbc, dbd)), 0.5L);
        rep->xw = d6_scale(d6_sub(d6_sub(dad, dab), dbd), 0.5L);
        rep->wy = d6_scale(d6_sub(d6_add(dbd, dcd), dbc), 0.5L);
        rep->nm = d6_mul(rep->xy, rep->inv_r_rho);
    }
}


// ------------------------------------------------------------------------------------------------
// Distance-only UTT4 logarithmic-integral kernel
// ------------------------------------------------------------------------------------------------

/**
 * @brief Scalar value and analytic partial derivatives of one contracted ordered node.
 *
 * The contracted fourth-derivative expression is evaluated once as a scalar. Its partial
 * derivatives with respect to the reduced generator variables and node geometry were derived
 * algebraically and common-subexpression eliminated. The six squared-distance derivatives are
 * reconstructed afterward by one final chain rule, avoiding Dual6 propagation through the full
 * contraction.
 */
typedef struct {
    long double alpha;
    long double beta;
    long double gamma;
    long double delta;
    long double A;
    long double B;

    long double B2;
    long double A2;
    long double delta_gamma;
    long double alpha_beta;
    long double alpha_beta_delta_gamma;
    long double sixty_alpha_beta_delta_gamma;
    long double A2_delta_gamma;
    long double alpha_beta_B2;
    long double two_A;
    long double two_A_B;
    long double A_B;
    long double two_B;
    long double three_gamma;
    long double three_delta;
} UTT4LnNodeParameters;


static inline UTT4LnNodeParameters utt4_ln_node_parameters(
    long double alpha, long double A, long double gamma, long double B)
{
    UTT4LnNodeParameters p;
    p.alpha = alpha;
    p.beta = 1.0L - alpha;
    p.gamma = gamma;
    p.delta = 1.0L - gamma;
    p.A = A;
    p.B = B;
    p.B2 = B*B;
    p.A2 = A*A;
    p.delta_gamma = p.delta*gamma;
    p.alpha_beta = alpha*p.beta;
    p.alpha_beta_delta_gamma = p.alpha_beta*p.delta_gamma;
    p.sixty_alpha_beta_delta_gamma = 60.0L*p.alpha_beta_delta_gamma;
    p.A2_delta_gamma = p.A2*p.delta_gamma;
    p.alpha_beta_B2 = p.alpha_beta*p.B2;
    p.two_A = 2.0L*A;
    p.two_A_B = B*p.two_A;
    p.A_B = A*B;
    p.two_B = 2.0L*B;
    p.three_gamma = 3.0L*gamma;
    p.three_delta = 3.0L*p.delta;
    return p;
}


typedef struct {
    long double value;
    long double d_s;
    long double d_q;
    long double d_nm;
    long double d_np;
    long double d_pm;
    long double d_inv_r;
    long double d_inv_rho;
} UTT4LnNodePartials;


static inline UTT4LnNodePartials utt4_ln_contracted_node_partials(
    const long double D[UTT4_LN_NUM_GENERATOR_DERIVATIVES],
    const UTT4LnNodeParameters *node, long double nm, long double npv, long double pm,
    long double q, long double inv_r, long double inv_rho)
{
    const long double alpha = node->alpha;
    const long double beta = node->beta;
    const long double gamma = node->gamma;
    const long double delta = node->delta;
    const long double A = node->A;
    const long double B = node->B;
    const long double D02 = UTT4_LN_D(D, 0, 2);
    const long double D03 = UTT4_LN_D(D, 0, 3);
    const long double D04 = UTT4_LN_D(D, 0, 4);
    const long double D05 = UTT4_LN_D(D, 0, 5);
    const long double D11 = UTT4_LN_D(D, 1, 1);
    const long double D12 = UTT4_LN_D(D, 1, 2);
    const long double D13 = UTT4_LN_D(D, 1, 3);
    const long double D14 = UTT4_LN_D(D, 1, 4);
    const long double D20 = UTT4_LN_D(D, 2, 0);
    const long double D21 = UTT4_LN_D(D, 2, 1);
    const long double D22 = UTT4_LN_D(D, 2, 2);
    const long double D23 = UTT4_LN_D(D, 2, 3);
    const long double D30 = UTT4_LN_D(D, 3, 0);
    const long double D31 = UTT4_LN_D(D, 3, 1);
    const long double D32 = UTT4_LN_D(D, 3, 2);
    const long double D40 = UTT4_LN_D(D, 4, 0);
    const long double D41 = UTT4_LN_D(D, 4, 1);
    const long double D50 = UTT4_LN_D(D, 5, 0);

    UTT4LnNodePartials out;
    const long double t0 = A*inv_rho;
    const long double t1 = B*inv_r;
    const long double t2 = D20*t1;
    const long double t3 = node->B2;
    const long double t4 = D30*t3;
    const long double t5 = A*inv_r;
    const long double t6 = node->A2;
    const long double t7 = D30*t6;
    const long double t8 = B*inv_rho;
    const long double t9 = node->delta_gamma;
    const long double t10 = t5*t9;
    const long double t11 = 4*D11;
    const long double t12 = node->alpha_beta;
    const long double t13 = t12*t8;
    const long double t14 = node->alpha_beta_delta_gamma;
    const long double t15 = node->sixty_alpha_beta_delta_gamma;
    const long double t16 = node->A2_delta_gamma;
    const long double t17 = 2*D21;
    const long double t18 = node->alpha_beta_B2;
    const long double t19 = D21*alpha;
    const long double t20 = node->two_A;
    const long double t21 = node->two_A_B;
    const long double t22 = nm*t21;
    const long double t23 = node->A_B;
    const long double t24 = nm*t23;
    const long double t25 = gamma*t24;
    const long double t26 = node->two_B;
    const long double t27 = npv*t0;
    const long double t28 = t26*t27;
    const long double t29 = D21*beta;
    const long double t30 = delta*t24;
    const long double t31 = D21*delta;
    const long double t32 = pm*t1;
    const long double t33 = t20*t32;
    const long double t34 = D21*gamma;
    const long double t35 = A*npv;
    const long double t36 = t35*t9;
    const long double t37 = alpha*t36;
    const long double t38 = 20*D12;
    const long double t39 = beta*t36;
    const long double t40 = q*t10;
    const long double t41 = 4*D12;
    const long double t42 = B*pm;
    const long double t43 = t12*t42;
    const long double t44 = delta*t43;
    const long double t45 = gamma*t43;
    const long double t46 = q*t13;
    const long double t47 = q*t9;
    const long double t48 = 80*t12*t47;
    const long double t49 = (nm*nm);
    const long double t50 = t0*t49;
    const long double t51 = alpha*t42;
    const long double t52 = nm*t0;
    const long double t53 = t51*t52;
    const long double t54 = beta*t42;
    const long double t55 = t52*t54;
    const long double t56 = delta*t35;
    const long double t57 = nm*t1;
    const long double t58 = t56*t57;
    const long double t59 = gamma*t35;
    const long double t60 = t57*t59;
    const long double t61 = 4*D22;
    const long double t62 = t51*t56;
    const long double t63 = alpha*q;
    const long double t64 = t25*t63;
    const long double t65 = beta*q;
    const long double t66 = t30*t65;
    const long double t67 = t54*t59;
    const long double t68 = 8*q;
    const long double t69 = D13*t68;
    const long double t70 = nm*pm*t20;
    const long double t71 = D31*t3*t70;
    const long double t72 = (npv*npv);
    const long double t73 = t10*t72;
    const long double t74 = nm*npv*t26;
    const long double t75 = D31*t6*t74;
    const long double t76 = (pm*pm);
    const long double t77 = t13*t76;
    const long double t78 = (q*q);
    const long double t79 = 16*t14*t78;
    const long double t80 = t3*t49*t6;
    const long double t81 = t4*t49;
    const long double t82 = t49*t7;
    const long double t83 = t16*t72;
    const long double t84 = t18*t76;
    const long double t85 = t49 + 1;
    const long double t86 = 4*t10 + 4*t13;
    const long double t87 = t51 - t54 + t56 - t59;
    const long double t88 = t23*(-t0 + t1*t49 - t1 + t50);
    const long double t89 = t68*(t37 - t39 + t44 - t45);
    const long double t90 = -4*t62 + 4*t64 + 4*t66 - 4*t67 + 4*t83 + 4*t84;
    const long double t91 = 20*t37 - 20*t39 + 4*t40 + 20*t44 - 20*t45 + 4*t46 - 4*t73 - 4*t77;
    const long double t92 = B*t27;
    const long double t93 = A*t32;
    const long double t94 = node->three_gamma;
    const long double t95 = node->three_delta;
    const long double t96 = 2*alpha*t24*t94 - 2*alpha*t30 - 2*alpha*t92 + 2*beta*t24*t95 - 2*beta*t25 + 2*beta*t92 - 2*delta*t93 + 2*gamma*t93 + 2*t16 + 2*t18 + 2*t53 - 2*t55 + 2*t58 - 2*t60;
    const long double t97 = 8*D13;
    const long double t98 = D30*nm;
    const long double t99 = inv_rho*pm;
    const long double t100 = inv_r*npv;
    const long double t101 = 2*D22;
    const long double t102 = gamma*t101;
    const long double t103 = delta*t101;
    const long double t104 = 10*D12;
    const long double t105 = t104*t9;
    const long double t106 = t41*t9;
    const long double t107 = 4*D13;
    const long double t108 = t107*t9;
    const long double t109 = t104*t12;
    const long double t110 = D31*t24;
    const long double t111 = t12*t41;
    const long double t112 = q*t107*t12;
    const long double t113 = D20*t8;
    const long double t114 = 2*t42;
    const long double t115 = D20*t5;
    const long double t116 = 2*t35;
    out.value = D02*t15 + D03*t48 + D04*t79 + D40*t80 + alpha*t71 - beta*t71 + delta*t19*t22 + delta*t75 + gamma*t22*t29 - gamma*t75 + t0*t2 - t10*t11 - t11*t13 - t16*t17 - t17*t18 - t17*t53 + t17*t55 - t17*t58 + t17*t60 - 6*t19*t25 + t19*t28 + t2*t50 - t28*t29 - 6*t29*t30 + t31*t33 - t33*t34 - t37*t38 - t37*t69 + t38*t39 - t38*t44 + t38*t45 + t39*t69 + t4*t5 - t40*t41 - t41*t46 + t41*t73 + t41*t77 - t44*t69 + t45*t69 - t5*t81 + t61*t62 - t61*t64 - t61*t66 + t61*t67 - t61*t83 - t61*t84 + t7*t8 - t8*t82;
    out.d_s = D12*t15 + D13*t48 + D14*t79 - D21*t86 - D22*t91 - D23*t89 + D30*t0*t1*t85 - D31*t96 - D32*t90 - D40*t88 + D41*t22*t87 + D50*t80;
    out.d_q = A*B*D21*inv_r*inv_rho*t85 + 2*A*B*D32*nm*t87 + 8*A*D13*beta*delta*gamma*npv + 8*B*D13*alpha*beta*gamma*pm + 140*D03*alpha*beta*delta*gamma + 112*D04*alpha*beta*delta*gamma*q + 16*D05*alpha*beta*delta*gamma*t78 - D12*t86 - D13*t91 - D14*t89 - D22*t96 - D23*t90 - D31*t88 + D41*t3*t49*t6 - alpha*t25*t61 - beta*t30*t61 - t10*t41 - t13*t41 - t37*t97 - t44*t97;
    out.d_nm = t21*(A*B*D40*nm + A*D31*delta*npv + B*D31*alpha*pm + D20*inv_r*inv_rho*nm + D21*alpha*delta + D21*beta*gamma + D21*beta*inv_rho*pm + D21*gamma*inv_r*npv - D31*t54 - D31*t59 - t0*t98 - t1*t98 - t100*t31 - t102*t63 - t103*t65 - t19*t94 - t19*t99 - t29*t95);
    out.d_np = t20*(-D31*t25 + D31*t30 - alpha*t105 + beta*t105 + t100*t106 + t102*t54 + t103*t51 - t108*t63 + t108*t65 + t19*t8 - t29*t8 - t31*t57 + t34*t57 - t36*t61);
    out.d_pm = t26*(alpha*t101*t56 + alpha*t110 + beta*t101*t59 - beta*t110 - delta*t109 - delta*t112 + gamma*t109 + gamma*t112 + t111*t99 - t19*t52 + t29*t52 + t31*t5 - t34*t5 - t43*t61);
    out.d_inv_r = A*(t106*t72 - t11*t9 + t113*t49 + t113 + t114*t31 - t114*t34 - t31*t74 + t34*t74 + t4 - t41*t47 - t81);
    out.d_inv_rho = B*(-q*t111 - t11*t12 + t111*t76 + t115*t49 + t115 + t116*t19 - t116*t29 - t19*t70 + t29*t70 + t7 - t82);
    return out;
}


/**
 * @brief Evaluate one ordered logarithmic-integral representative at a quadrature node.
 *
 * Position-only geometry and its six squared-distance derivatives are precomputed once per
 * quadruple. At each node the contracted scalar expression is evaluated analytically, followed by
 * a single chain rule from (s,q,nm,np,pm,1/r,1/rho) to the six squared pair distances.
 */
static Dual6 utt4_ln_ordered_node(const UTT4LnPairGeometry *geo,
    const UTT4LnRepresentativeGeometry *rep, const UTT4LnNodeParameters *node)
{
    const long double alpha = node->alpha;
    const long double beta = node->beta;
    const long double gamma = node->gamma;
    const long double delta = node->delta;
    const long double A = node->A;
    const long double B = node->B;

    const long double dab = geo->d2[rep->pab];
    const long double dcd = geo->d2[rep->pcd];
    const long double dac = geo->d2[rep->pac];
    const long double dad = geo->d2[rep->pad];
    const long double dbc = geo->d2[rep->pbc];
    const long double dbd = geo->d2[rep->pbd];

    const long double q_ac = alpha*gamma;
    const long double q_ad = alpha*delta;
    const long double q_bc = beta*gamma;
    const long double q_bd = beta*delta;
    const long double q_ab = -alpha*beta;
    const long double q_cd = -gamma*delta;

    long double q = q_ac*dac + q_ad*dad + q_bc*dbc + q_bd*dbd
        + q_ab*dab + q_cd*dcd;
    const long double s = A*rep->r.v + B*rep->rho.v;

    if (!(s > 0.0L) || q < 0.0L) {
        if (q > -1e-26L*fmaxl(1.0L, s*s))
            q = 0.0L;
        else {
            Dual6 bad = d6_const(NAN);
            for (int k = 0; k < UTT4_LN_NPAIR; ++k)
                bad.g[k] = NAN;
            return bad;
        }
    }

    const long double xp = rep->xw.v + alpha*dab - gamma*rep->xy.v;
    const long double py = rep->wy.v + alpha*rep->xy.v - gamma*dcd;
    const long double inv_r = rep->inv_r.v;
    const long double inv_rho = rep->inv_rho.v;
    const long double npv = xp*inv_r;
    const long double pm = py*inv_rho;

    long double D[UTT4_LN_NUM_GENERATOR_DERIVATIVES];
    utt4_ln_generator_derivatives(s, q, D);

    const UTT4LnNodePartials p = utt4_ln_contracted_node_partials(
        D, node, rep->nm.v, npv, pm, q, inv_r, inv_rho);

    long double dq[UTT4_LN_NPAIR] = {0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
    dq[rep->pac] = q_ac;
    dq[rep->pad] = q_ad;
    dq[rep->pbc] = q_bc;
    dq[rep->pbd] = q_bd;
    dq[rep->pab] = q_ab;
    dq[rep->pcd] = q_cd;

    Dual6 out;
    out.v = p.value;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k) {
        long double dxp = rep->xw.g[k] - gamma*rep->xy.g[k];
        long double dpy = rep->wy.g[k] + alpha*rep->xy.g[k];
        if (k == rep->pab)
            dxp += alpha;
        if (k == rep->pcd)
            dpy -= gamma;

        const long double ds = A*rep->r.g[k] + B*rep->rho.g[k];
        const long double dnp = dxp*inv_r + xp*rep->inv_r.g[k];
        const long double dpm = dpy*inv_rho + py*rep->inv_rho.g[k];

        out.g[k] = p.d_s*ds + p.d_q*dq[k] + p.d_nm*rep->nm.g[k]
                 + p.d_np*dnp + p.d_pm*dpm
                 + p.d_inv_r*rep->inv_r.g[k] + p.d_inv_rho*rep->inv_rho.g[k];
    }
    return out;
}

/**
 * @brief Sum the six inequivalent representatives of the 24 ordered permutations.
 *
 * The factor four reconstructs the full ordered sum using I_ab;cd = I_ba;dc = I_cd;ab = I_dc;ba.
 * The 1/pi normalization is included here.
 */
static Dual6 utt4_ln_sum_node(const UTT4LnPairGeometry *geo, long double alpha, long double A,
    long double gamma, long double B)
{
    const UTT4LnNodeParameters node = utt4_ln_node_parameters(alpha, A, gamma, B);
    Dual6 sum = d6_zero();
    for (int p = 0; p < 6; ++p)
        sum = d6_add(sum, utt4_ln_ordered_node(geo, &geo->representative[p], &node));
    return d6_scale(sum, -4.0L/M_PI);
}


// ------------------------------------------------------------------------------------------------
// Gauss-Legendre quadrature
// ------------------------------------------------------------------------------------------------

typedef struct {
    long double value;
    long double d_d2[UTT4_LN_NPAIR];
    long double grad[UTT4_LN_INTEGRAL_NUM_BODIES][3];
    long double error_estimate;
    long double max_component_error;
    long long node_evals;
} UTT4LnWorkResult;


static UTT4LnWorkResult utt4_ln_fixed_quadrature(const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3],
    const UTT4LnQuadratureRule *q, int use_openmp)
{
#ifndef _OPENMP
    (void)use_openmp;
#endif
    UTT4LnPairGeometry geo;
    utt4_ln_build_pair_geometry(pos, &geo);
    Dual6 total = d6_zero();
#ifdef _OPENMP
    if (use_openmp && omp_get_max_threads() > 1 && !omp_in_parallel()) {
        #pragma omp parallel
        {
            Dual6 local = d6_zero();
            #pragma omp for collapse(2) schedule(static)
            for (int i = 0; i < q->base.n; ++i) {
                for (int j = 0; j < q->base.n; ++j) {
                    Dual6 f = utt4_ln_sum_node(&geo,
                                q->mix[i], q->half_sin_theta[i],
                                q->mix[j], q->half_sin_theta[j]);
                    local = d6_add(local, d6_scale(f, q->base.weight[i]*q->base.weight[j]));
                }
            }
            #pragma omp critical(prod_force_sum)
            total = d6_add(total, local);
        }
    } else
#endif
    {
        for (int i = 0; i < q->base.n; ++i)
            for (int j = 0; j < q->base.n; ++j) {
                Dual6 f = utt4_ln_sum_node(&geo,
                            q->mix[i], q->half_sin_theta[i],
                            q->mix[j], q->half_sin_theta[j]);
                total = d6_add(total, d6_scale(f, q->base.weight[i]*q->base.weight[j]));
            }
    }

    UTT4LnWorkResult out;
    memset(&out, 0, sizeof(out));
    out.value = total.v;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p)
        out.d_d2[p] = total.g[p];
    out.node_evals = (long long)q->base.n*(long long)q->base.n;

    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const int a = utt4_ln_pair_body_i[p], b = utt4_ln_pair_body_j[p];
        for (int k = 0; k < 3; ++k) {
            const long double dx = (long double)pos[a][k] - (long double)pos[b][k];
            const long double c = 2.0L*out.d_d2[p]*dx;
            out.grad[a][k] += c;
            out.grad[b][k] -= c;
        }
    }
    return out;
}


// ------------------------------------------------------------------------------------------------
// Adaptive Gauss-Kronrod quadrature
// ------------------------------------------------------------------------------------------------

// Kronrod-15 nodes and weights on [-1,1], with the embedded Gauss-7 weights.
static const long double gk15_x[15] = {
    -0.991455371120812639206854697526329L,
    -0.949107912342758524526189684047851L,
    -0.864864423359769072789712788640926L,
    -0.741531185599394439863864773280788L,
    -0.586087235467691130294144838258730L,
    -0.405845151377397166906606412076961L,
    -0.207784955007898467600689403773245L,
     0.0L,
     0.207784955007898467600689403773245L,
     0.405845151377397166906606412076961L,
     0.586087235467691130294144838258730L,
     0.741531185599394439863864773280788L,
     0.864864423359769072789712788640926L,
     0.949107912342758524526189684047851L,
     0.991455371120812639206854697526329L
};

static const long double gk15_wk[15] = {
    0.022935322010529224963732008058970L,
    0.063092092629978553290700663189204L,
    0.104790010322250183839876322541518L,
    0.140653259715525918745189590510238L,
    0.169004726639267902826583426598550L,
    0.190350578064785409913256402421014L,
    0.204432940075298892414161999234649L,
    0.209482141084727828012999174891714L,
    0.204432940075298892414161999234649L,
    0.190350578064785409913256402421014L,
    0.169004726639267902826583426598550L,
    0.140653259715525918745189590510238L,
    0.104790010322250183839876322541518L,
    0.063092092629978553290700663189204L,
    0.022935322010529224963732008058970L
};

static const long double gk15_wg[15] = {
    0.0L,
    0.129484966168869693270611432679082L,
    0.0L,
    0.279705391489276667901467771423780L,
    0.0L,
    0.381830050505118944950369775488975L,
    0.0L,
    0.417959183673469387755102040816327L,
    0.0L,
    0.381830050505118944950369775488975L,
    0.0L,
    0.279705391489276667901467771423780L,
    0.0L,
    0.129484966168869693270611432679082L,
    0.0L
};


typedef struct {
    Dual6 kronrod;
    Dual6 gauss;
    long long node_evals;
} UTT4LnCellEstimate;


typedef struct {
    UTT4LnPairGeometry geo;
    long double rel_tol;
    int max_depth;
    int parallel_root;
    long double length2_scale;
} UTT4LnAdaptiveContext;


typedef struct {
    Dual6 result;
    long double abs_error;
    long double max_component_error;
    long long node_evals;
} UTT4LnAdaptiveAccum;


static UTT4LnCellEstimate utt4_ln_cell_estimate(const UTT4LnAdaptiveContext *ctx,
    long double ta, long double tb, long double pa, long double pb)
{
    UTT4LnCellEstimate out;
    out.kronrod = d6_zero();
    out.gauss = d6_zero();
    out.node_evals = 225;
    const long double tm = 0.5L*(ta + tb), th = 0.5L*(tb - ta);
    const long double pm = 0.5L*(pa + pb), ph = 0.5L*(pb - pa);

    long double alpha[15], A[15], wk_t[15], wg_t[15];
    long double gamma[15], B[15], wk_p[15], wg_p[15];
    for (int i = 0; i < 15; ++i) {
        const long double theta = tm + th*gk15_x[i];
        const long double phi = pm + ph*gk15_x[i];
        alpha[i] = 0.5L*(1.0L + cosl(theta));
        A[i] = 0.5L*sinl(theta);
        wk_t[i] = th*gk15_wk[i];
        wg_t[i] = th*gk15_wg[i];
        gamma[i] = 0.5L*(1.0L + cosl(phi));
        B[i] = 0.5L*sinl(phi);
        wk_p[i] = ph*gk15_wk[i];
        wg_p[i] = ph*gk15_wg[i];
    }

    Dual6 ksum = d6_zero(), gsum = d6_zero();
    for (int i = 0; i < 15; ++i) {
        for (int j = 0; j < 15; ++j) {
            Dual6 f = utt4_ln_sum_node(&ctx->geo, alpha[i], A[i], gamma[j], B[j]);
            ksum = d6_add(ksum, d6_scale(f, wk_t[i]*wk_p[j]));
            if (wg_t[i] != 0.0L && wg_p[j] != 0.0L)
                gsum = d6_add(gsum, d6_scale(f, wg_t[i]*wg_p[j]));
        }
    }
    out.kronrod = ksum;
    out.gauss = gsum;
    return out;
}


static int utt4_ln_cell_accept(const UTT4LnAdaptiveContext *ctx, const UTT4LnCellEstimate *e,
    long double local_abs, long double *value_err, long double *max_comp_err)
{
    const long double ev = fabsl(e->kronrod.v - e->gauss.v);
    long double maxe = ev;
    int ok = (ev <= local_abs + ctx->rel_tol*fabsl(e->kronrod.v));
    const long double grad_abs = local_abs/fmaxl(ctx->length2_scale, 1e-30L);
    for (int k = 0; k < UTT4_LN_NPAIR; ++k) {
        const long double eg = fabsl(e->kronrod.g[k] - e->gauss.g[k]);
        if (eg > maxe) maxe = eg;
        if (eg > grad_abs + ctx->rel_tol*fabsl(e->kronrod.g[k])) ok = 0;
    }
    *value_err=ev;
    *max_comp_err=maxe;
    return ok;
}

static UTT4LnAdaptiveAccum utt4_ln_adaptive_cell(const UTT4LnAdaptiveContext *ctx, long double ta,
    long double tb, long double pa, long double pb, long double local_abs, int depth)
{
    UTT4LnAdaptiveAccum acc;
    memset(&acc, 0, sizeof(acc));
    UTT4LnCellEstimate e = utt4_ln_cell_estimate(ctx, ta, tb, pa, pb);
    acc.node_evals = e.node_evals;
    long double verr = 0.0L, maxerr = 0.0L;
    if (utt4_ln_cell_accept(ctx, &e, local_abs, &verr, &maxerr) || depth >= ctx->max_depth) {
        acc.result = e.kronrod;
        acc.abs_error = verr;
        acc.max_component_error = maxerr;
        return acc;
    }

    const long double tm = 0.5L*(ta + tb), pm = 0.5L*(pa + pb);
    const long double bounds[4][4] = {
        {ta, tm, pa, pm}, {tm, tb, pa, pm}, {ta, tm, pm, pb}, {tm, tb, pm, pb}
    };
    UTT4LnAdaptiveAccum child[4];
    memset(child, 0, sizeof(child));

#ifdef _OPENMP
    if (ctx->parallel_root && depth == 0 && omp_get_max_threads() > 1 && !omp_in_parallel()) {
        #pragma omp parallel for schedule(static)
        for (int c = 0; c < 4; ++c) {
            UTT4LnAdaptiveContext sub = *ctx;
            sub.parallel_root = 0;
            child[c] = utt4_ln_adaptive_cell(&sub, bounds[c][0], bounds[c][1],bounds[c][2],
                bounds[c][3], local_abs*0.25L, depth + 1);
        }
    } else
#endif
    {
        for (int c = 0; c < 4; ++c)
            child[c] = utt4_ln_adaptive_cell(ctx, bounds[c][0], bounds[c][1], bounds[c][2],
                bounds[c][3], local_abs*0.25L, depth + 1);
    }

    acc.result = d6_zero();
    acc.node_evals = e.node_evals;
    for (int c = 0; c < 4; ++c) {
        acc.result = d6_add(acc.result, child[c].result);
        acc.abs_error += child[c].abs_error;
        if (child[c].max_component_error > acc.max_component_error)
            acc.max_component_error = child[c].max_component_error;
        acc.node_evals += child[c].node_evals;
    }
    return acc;
}


static UTT4LnWorkResult utt4_ln_adaptive_quadrature(
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3], long double rel_tol, long double abs_tol,
    int max_depth, int parallel_root)
{
    UTT4LnAdaptiveContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    utt4_ln_build_pair_geometry(pos, &ctx.geo);
    ctx.rel_tol = rel_tol;
    ctx.max_depth = max_depth;
    ctx.parallel_root = parallel_root;
    long double maxd2 = 0.0L;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) if (ctx.geo.d2[p] > maxd2) maxd2 = ctx.geo.d2[p];
    ctx.length2_scale = maxd2;

    UTT4LnAdaptiveAccum a = utt4_ln_adaptive_cell(&ctx, 0.0L, M_PI, 0.0L, M_PI, abs_tol, 0);
    UTT4LnWorkResult out;
    memset(&out, 0, sizeof(out));
    out.value = a.result.v;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p)
        out.d_d2[p] = a.result.g[p];
    out.error_estimate = a.abs_error;
    out.max_component_error = a.max_component_error;
    out.node_evals = a.node_evals;

    // Reconstruct all Cartesian derivatives from dS/d(r_ij^2)
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const int i = utt4_ln_pair_body_i[p], j = utt4_ln_pair_body_j[p];
        for (int k = 0; k < 3; ++k) {
            const long double dx = (long double)pos[i][k] - (long double)pos[j][k];
            const long double contribution = 2.0L*out.d_d2[p]*dx;
            out.grad[i][k] += contribution;
            out.grad[j][k] -= contribution;
        }
    }

    return out;
}


// ------------------------------------------------------------------------------------------------
// Error estimation and automatic order refinement
// ------------------------------------------------------------------------------------------------

/**
 * @brief Compare two fixed-order results using the value and all six squared-distance derivatives.
 *
 * A refinement is requested when any component exceeds its absolute-plus-relative tolerance.
 */
static int utt4_ln_needs_refinement(const UTT4LnWorkResult *low, const UTT4LnWorkResult *high,
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3], long double rtol, long double atol,
    long double *value_error, long double *max_component_error, long double *value_rel_est,
    long double *max_component_rel_est, long double *worst_tol_ratio)
{
    const long double ev = fabsl(high->value-low->value);
    const long double value_scale = fmaxl(fabsl(high->value), 1e-300L);
    const long double vr = ev/value_scale;
    const long double value_allowed = atol + rtol*fabsl(high->value);
    long double worst = (value_allowed > 0.0L) ? ev/value_allowed : (ev == 0.0L ? 0.0L : HUGE_VALL);
    long double maxe = ev, maxrel = vr;
    int refine = worst > 1.0L;
    long double maxd2 = 0.0L;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const int i = utt4_ln_pair_body_i[p], j = utt4_ln_pair_body_j[p];
        long double d2 = 0.0L;
        for (int k = 0; k < 3; ++k) {const long double x = (long double)pos[i][k] - pos[j][k]; d2 += x*x; }
        if (d2 > maxd2)maxd2 = d2;
    }
    const long double datol = atol/fmaxl(maxd2, 1e-30L);
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const long double e = fabsl(high->d_d2[p] - low->d_d2[p]);
        const long double scale = fmaxl(fabsl(high->d_d2[p]), 1e-300L);
        const long double rel = e/scale;
        const long double allowed = datol + rtol*fabsl(high->d_d2[p]);
        const long double ratio = (allowed > 0.0L) ? e/allowed : (e == 0.0L ? 0.0L : HUGE_VALL);
        if (e > maxe)maxe = e;
        if (rel > maxrel)maxrel = rel;
        if (ratio > worst)worst = ratio;
        if (ratio > 1.0L)refine = 1;
    }
    if (value_error)*value_error = ev;
    if (max_component_error)*max_component_error = maxe;
    if (value_rel_est)*value_rel_est = vr;
    if (max_component_rel_est)*max_component_rel_est = maxrel;
    if (worst_tol_ratio)*worst_tol_ratio = worst;
    return refine;
}


static int integral_round_up_even_order(long double x)
{
    int n = (int)ceill(x);
    if (n < 4)n = 4;
    if (n & 1)++n;
    return n;
}

static int integral_initial_high_order(long double rtol, int min_order, int max_order)
{
    long double digits = -log10l(rtol);
    if (digits < 1.0L)digits = 1.0L;
    int high = integral_round_up_even_order(2.7L*digits);
    if (high < min_order + 2)high = min_order + 2;
    if (high > max_order)high = max_order;
    return high;
}

static int integral_previous_order(int high, int min_order)
{
    int low = integral_round_up_even_order(0.88L*(long double)high);
    if (low >= high)low = high - 2;
    if (low < min_order)low = min_order;
    if (low >= high)low = high - 2;
    return low;
}

static int integral_next_order(int high, int max_order)
{
    int next = (high < 64) ? high + 4 : integral_round_up_even_order(1.20L*(long double)high);
    if (next <= high)next = high + 2;
    if (next > max_order)next = max_order;
    return next;
}


static int utt4_ln_adaptive_meets_target(const UTT4LnWorkResult *r,
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3], long double rtol, long double atol,
    long double *worst_ratio)
{
    const long double value_allowed = atol + rtol*fabsl(r->value);
    long double worst = (value_allowed > 0.0L) ? r->error_estimate/value_allowed : 0.0L;
    long double maxd2 = 0.0L;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const int i = utt4_ln_pair_body_i[p], j = utt4_ln_pair_body_j[p];
        long double d2 = 0.0L;
        for (int k = 0; k < 3; ++k) {
            const long double dx = (long double)pos[i][k] - (long double)pos[j][k];
            d2 += dx*dx;
        }
        if (d2 > maxd2)maxd2 = d2;
    }
    const long double datol = atol/fmaxl(maxd2, 1e-30L);
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const long double allowed = datol + rtol*fabsl(r->d_d2[p]);
        const long double ratio = (allowed > 0.0L) ? r->max_component_error/allowed : 0.0L;
        if (ratio > worst)worst = ratio;
    }
    if (worst_ratio)*worst_ratio = worst;
    return worst <= 1.0L;
}


// ------------------------------------------------------------------------------------------------
// Public UTT4 logarithmic-integral evaluator
// ------------------------------------------------------------------------------------------------

/**
 * @brief Evaluate the logarithmic four-body integral entering UTT4 and all position derivatives.
 *
 * The routine first uses automatically selected fixed Gauss-Legendre rules and compares two
 * successive orders for both the integral value and the six squared-distance derivatives. If
 * the requested tolerance is still not reached at max_order and adaptive evaluation is enabled,
 * an adaptive Gauss-Kronrod fallback is used.
 *
 * @param[in]  pos       Four local particle positions in three dimensions.
 * @param[in]  settings  Accuracy, order-refinement, adaptive, and OpenMP settings.
 * @param[out] result    Integral value, Cartesian derivatives, and quadrature diagnostics.
 * @return 0 on success, nonzero on invalid input or quadrature setup failure.
 */
int utt4_ln_integral_evaluate(const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3],
    const NumericalIntegralSettings *settings, UTT4LnIntegralResult *out)
{
    if (!pos || !settings || !out) return -1;
    memset(out, 0, sizeof(*out));
    if (!(settings->rel_tol > 0.0L) || !(settings->abs_tol > 0.0L) ||
       settings->min_order < 4 || settings->max_order < settings->min_order ||
       settings->max_depth < 0) return -1;

    int high_order = integral_initial_high_order(settings->rel_tol, settings->min_order,
        settings->max_order);
    int low_order = integral_previous_order(high_order, settings->min_order);
    if (low_order < 4 || high_order <= low_order) {
        low_order = settings->min_order;
        high_order = low_order + 2;
        if (high_order > settings->max_order)
            high_order = settings->max_order;
    }
    if (high_order <= low_order)
        return -1;

    const UTT4LnQuadratureRule *qlow = utt4_ln_quadrature_rule_get(low_order);
    const UTT4LnQuadratureRule *qhigh = utt4_ln_quadrature_rule_get(high_order);
    if (qlow == NULL || qhigh == NULL)
        return -1;

    UTT4LnWorkResult low = utt4_ln_fixed_quadrature(pos, qlow, settings->use_openmp);
    UTT4LnWorkResult high = utt4_ln_fixed_quadrature(pos, qhigh, settings->use_openmp);
    out->diagnostics.node_evals = low.node_evals + high.node_evals;

    long double value_error = 0.0L, max_component_error = 0.0L;
    long double value_rel = 0.0L, max_component_rel = 0.0L, worst = 0.0L;
    int refine = utt4_ln_needs_refinement(&low, &high, pos, settings->rel_tol, settings->abs_tol,
        &value_error, &max_component_error, &value_rel, &max_component_rel, &worst);
    while (refine && high_order < settings->max_order) {
        low_order = high_order;
        low = high;
        high_order = integral_next_order(high_order, settings->max_order);
        if (high_order <= low_order)
            break;
        const UTT4LnQuadratureRule *qnext = utt4_ln_quadrature_rule_get(high_order);
        if (qnext == NULL)
            return -1;
        high = utt4_ln_fixed_quadrature(pos, qnext, settings->use_openmp);
        out->diagnostics.node_evals += high.node_evals;
        refine = utt4_ln_needs_refinement(&low, &high, pos, settings->rel_tol, settings->abs_tol,
            &value_error, &max_component_error, &value_rel, &max_component_rel, &worst);
    }

    UTT4LnWorkResult final = high;
    out->diagnostics.low_order = low_order;
    out->diagnostics.high_order = high_order;
    out->diagnostics.target_met = !refine;
    out->diagnostics.value_error = value_error;
    out->diagnostics.max_derivative_error = max_component_error;
    out->diagnostics.value_rel_est = value_rel;
    out->diagnostics.max_derivative_rel_est = max_component_rel;
    out->diagnostics.worst_tolerance_ratio = worst;

    if (refine && settings->adaptive) {
        final = utt4_ln_adaptive_quadrature(pos, settings->rel_tol, settings->abs_tol,
            settings->max_depth, settings->use_openmp);
        out->diagnostics.adaptive_used = 1;
        out->diagnostics.node_evals += final.node_evals;
        out->diagnostics.value_error = final.error_estimate;
        out->diagnostics.max_derivative_error = final.max_component_error;
        out->diagnostics.value_rel_est = final.error_estimate/fmaxl(fabsl(final.value), 1e-300L);
        long double maxrel = 0.0L;
        for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
            const long double rel = final.max_component_error/
                fmaxl(fabsl(final.d_d2[p]), 1e-300L);
            if (rel > maxrel)maxrel = rel;
        }
        out->diagnostics.max_derivative_rel_est = maxrel;
        out->diagnostics.target_met = utt4_ln_adaptive_meets_target(&final, pos, settings->rel_tol,
            settings->abs_tol, &out->diagnostics.worst_tolerance_ratio);
    }

    out->value = final.value;
    for (int b = 0; b < UTT4_LN_INTEGRAL_NUM_BODIES; ++b)
        for (int k = 0; k < 3; ++k)
            out->grad[b][k] = final.grad[b][k];
    return 0;
}


#undef UTT4_LN_NPAIR
