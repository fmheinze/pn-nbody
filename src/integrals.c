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

#define UTT4_PI ((utt4_real)3.141592653589793238462643383279502884)

#if UTT4_USE_LONG_DOUBLE
#define UTT4_SQRT   sqrtl
#define UTT4_FABS   fabsl
#define UTT4_COS    cosl
#define UTT4_SIN    sinl
#define UTT4_ATAN   atanl
#define UTT4_CEIL   ceill
#define UTT4_LOG10  log10l
#define UTT4_FMAX   fmaxl
#define UTT4_EPS    LDBL_EPSILON
#define UTT4_HUGE   HUGE_VALL
#else
#define UTT4_SQRT   sqrt
#define UTT4_FABS   fabs
#define UTT4_COS    cos
#define UTT4_SIN    sin
#define UTT4_ATAN   atan
#define UTT4_CEIL   ceil
#define UTT4_LOG10  log10
#define UTT4_FMAX   fmax
#define UTT4_EPS    DBL_EPSILON
#define UTT4_HUGE   HUGE_VAL
#endif


// ------------------------------------------------------------------------------------------------
// Generic quadrature machinery
// ------------------------------------------------------------------------------------------------

typedef struct {
    int n;
    utt4_real *x;
    utt4_real *weight;
} GaussLegendreRule;


// Build an n-point Gauss-Legendre rule on [a, b].
static int gauss_legendre_rule_init(int n, utt4_real a, utt4_real b, GaussLegendreRule *rule)
{
    rule->n = n;
    rule->x = calloc((size_t)n, sizeof(utt4_real));
    rule->weight = calloc((size_t)n, sizeof(utt4_real));
    if (!rule->x || !rule->weight)
        return -1;

    const utt4_real eps = 8.0 * UTT4_EPS;
    const utt4_real midpoint = 0.5 * (a + b);
    const utt4_real half_width = 0.5 * (b - a);
    const int m = (n + 1) / 2;

    for (int i = 0; i < m; ++i) {
        utt4_real z = UTT4_COS(UTT4_PI * ((utt4_real)i + 0.75) / ((utt4_real)n + 0.5));
        utt4_real z1;
        utt4_real pp = 0.0;
        do {
            utt4_real p1 = 1.0;
            utt4_real p2 = 0.0;
            for (int j = 1; j <= n; ++j) {
                const utt4_real p3 = p2;
                p2 = p1;
                p1 = (((2.0 * j - 1.0) * z * p2) - (j - 1.0) * p3) / (utt4_real)j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (UTT4_FABS(z - z1) > eps);

        const utt4_real w = 2.0 / ((1.0 - z * z) * pp * pp);
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

// Deepest geometric bisection toward an endpoint; 2^-24 covers any physical ratio.
#define UTT4_LN_MAX_GRADING_LEVELS 24

// Separation ratio at which grading starts to pay.
#define UTT4_LN_GRADING_MIN_RATIO 16.0

//Gauss-Legendre points per graded panel, as a fraction of the smooth-axis order.
#define UTT4_LN_PANEL_ORDER_DIVISOR 4
#define UTT4_LN_PANEL_ORDER_MIN 8

// Order past which the evaluator stops raising the degree and switches to graded panels.
#define UTT4_LN_PLAIN_SWITCH_ORDER 48

// Number of pair representatives
#define UTT4_LN_NPAIR 6


// Six representatives for the 24 ordered permutations.
static const int utt4_ln_sum_representatives[6][4] =
{
    {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1}, {0,3,1,2}, {0,3,2,1}
};


// One axis of the (theta, phi) square as a composite Gauss-Legendre rule on [0, pi].
typedef struct {
    int n;
    int base_order;
    int levels;
    utt4_real *half_sin_theta;
    utt4_real *mix;
    utt4_real *weight;
} UTT4LnQuadratureRule;


// Panel breakpoints on [0, pi]: pi/2 halved `levels` times toward 0, then mirrored about pi/2.
static int utt4_ln_axis_panels(int levels, utt4_real *bounds)
{
    bounds[0] = 0.0;

    // No layer to resolve: one panel spanning the whole range, identical to a plain rule.
    if (levels <= 0) {
        bounds[1] = UTT4_PI;
        return 1;
    }

    utt4_real edge = 0.5*UTT4_PI;
    for (int k = 0; k < levels; ++k)
        edge *= 0.5;

    for (int k = 0; k <= levels; ++k) {
        bounds[k + 1] = edge;
        edge *= 2.0;
    }

    const int half = levels + 1;
    for (int k = 1; k <= half; ++k)
        bounds[half + k] = UTT4_PI - bounds[half - k];

    return 2*half;
}


static int utt4_ln_quadrature_rule_init(int base_order, int levels, UTT4LnQuadratureRule *rule)
{
    memset(rule, 0, sizeof(*rule));
    if (base_order < 1 || levels < 0)
        return -1;

    const int panels = (levels > 0) ? 2*(levels + 1) : 1;
    const int n = panels*base_order;

    rule->half_sin_theta = calloc((size_t)n, sizeof(utt4_real));
    rule->mix = calloc((size_t)n, sizeof(utt4_real));
    rule->weight = calloc((size_t)n, sizeof(utt4_real));
    utt4_real *bounds = calloc((size_t)panels + 1, sizeof(utt4_real));
    if (!rule->half_sin_theta || !rule->mix || !rule->weight || !bounds) {
        free(rule->half_sin_theta);
        free(rule->mix);
        free(rule->weight);
        free(bounds);
        memset(rule, 0, sizeof(*rule));
        return -1;
    }

    utt4_ln_axis_panels(levels, bounds);

    int at = 0;
    for (int panel = 0; panel < panels; ++panel) {
        GaussLegendreRule gl;
        if (gauss_legendre_rule_init(base_order, bounds[panel], bounds[panel + 1], &gl) != 0) {
            free(rule->half_sin_theta);
            free(rule->mix);
            free(rule->weight);
            free(bounds);
            memset(rule, 0, sizeof(*rule));
            return -1;
        }

        for (int i = 0; i < base_order; ++i, ++at) {
            rule->half_sin_theta[at] = 0.5*UTT4_SIN(gl.x[i]);
            rule->mix[at] = 0.5*(1.0 + UTT4_COS(gl.x[i]));
            rule->weight[at] = gl.weight[i];
        }
        gauss_legendre_rule_free(&gl);
    }

    free(bounds);
    rule->n = n;
    rule->base_order = base_order;
    rule->levels = levels;
    return 0;
}


/**
 * @brief Process-lifetime cache of immutable axis rules, keyed on (base_order, levels).
 *
 * The same few rules are reused at virtually every UTT4 evaluation. Keeping them avoids repeated
 * root finding, trigonometric transforms, and heap allocation. Entries are never modified after
 * construction and intentionally live until process exit.
 */
typedef struct UTT4LnRuleCacheEntry {
    UTT4LnQuadratureRule rule;
    struct UTT4LnRuleCacheEntry *next;
} UTT4LnRuleCacheEntry;

static UTT4LnRuleCacheEntry *utt4_ln_rule_cache = NULL;
static int utt4_ln_rule_cache_cleanup_registered = 0;


static void utt4_ln_quadrature_rule_destroy(UTT4LnQuadratureRule *rule)
{
    free(rule->half_sin_theta);
    free(rule->mix);
    free(rule->weight);
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


static const UTT4LnQuadratureRule *utt4_ln_quadrature_rule_get(int base_order, int levels)
{
    const UTT4LnQuadratureRule *result = NULL;

#ifdef _OPENMP
    #pragma omp critical(utt4_ln_rule_cache_guard)
#endif
    {
        UTT4LnRuleCacheEntry *entry = utt4_ln_rule_cache;
        while (entry != NULL &&
               (entry->rule.base_order != base_order || entry->rule.levels != levels))
            entry = entry->next;

        if (entry == NULL) {
            entry = malloc(sizeof(*entry));
            if (entry != NULL) {
                memset(entry, 0, sizeof(*entry));
                if (utt4_ln_quadrature_rule_init(base_order, levels, &entry->rule) != 0) {
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
    utt4_real v;
    utt4_real g[UTT4_LN_NPAIR];
} Dual6;


static inline Dual6 d6_zero(void)
{
    Dual6 z;
    z.v = 0.0;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        z.g[k] = 0.0;
    return z;
}


static inline Dual6 d6_const(utt4_real x)
{
    Dual6 z = d6_zero();
    z.v = x;
    return z;
}


static inline Dual6 d6_seed(utt4_real x, int k)
{
    Dual6 z = d6_const(x);
    z.g[k] = 1.0;
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


static inline Dual6 d6_scale(Dual6 a, utt4_real c)
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
    const utt4_real iv = 1.0 / a.v;
    Dual6 r;
    r.v = iv;
    const utt4_real fac = -iv*iv;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k)
        r.g[k] = fac * a.g[k];
    return r;
}


static inline Dual6 d6_sqrt(Dual6 a)
{
    const utt4_real sv = UTT4_SQRT(a.v);
    Dual6 r;
    r.v = sv;
    const utt4_real fac = 0.5 / sv;
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
static void utt4_ln_H_derivatives(utt4_real y, utt4_real h[6])
{
    h[0] = 0.0;

    if (y < 0.10) {
        for (int k = 1; k <= 5; ++k)
            h[k] = 0.0;

        utt4_real power[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        for (int n = 1; n <= 128; ++n) {
            if (n <= 5 && n >= 2)
                power[n] = 1.0;

            const utt4_real nn = (utt4_real)n;
            utt4_real cn = 1.0/(nn*(4.0*nn*nn - 1.0));
            if (n & 1)
                cn = -cn;

            utt4_real falling = nn;
            int converged = (n > 17);
            const int kmax = (n < 5) ? n : 5;
            for (int k = 1; k <= kmax; ++k) {
                const utt4_real term = cn*falling*power[k];
                h[k] += term;
                if (converged &&
                    UTT4_FABS(term) >= 2.0*UTT4_EPS*UTT4_FMAX(1.0, UTT4_FABS(h[k])))
                    converged = 0;
                falling *= (utt4_real)(n - k);
            }

            if (converged)
                break;

            for (int k = 1; k <= kmax; ++k)
                power[k] *= y;
        }
        return;
    }

    const utt4_real inv_one_plus_y = 1.0/(1.0 + y);
    const utt4_real two_y = 2.0*y;
    const utt4_real t = UTT4_SQRT(y);
    utt4_real a[6];
    a[0] = UTT4_ATAN(t)/t;

    utt4_real rhs_scale = inv_one_plus_y;
    utt4_real factorial = 1.0;
    utt4_real sign = 1.0;
    for (int n = 0; n < 5; ++n) {
        a[n + 1] = (sign*factorial*rhs_scale - (2.0*n + 1.0)*a[n])/two_y;
        sign = -sign;
        factorial *= (utt4_real)(n + 1);
        rhs_scale *= inv_one_plus_y;
    }

    utt4_real log_derivative = inv_one_plus_y;
    for (int n = 1; n <= 5; ++n) {
        h[n] = (1.0 - y)*a[n] - (utt4_real)n*a[n - 1] + log_derivative;
        log_derivative *= -(utt4_real)n*inv_one_plus_y;
    }
}


/**
 * @brief Compute the generator derivatives needed by the contracted UTT4 kernel.
 *
 * Only 18 of the 21 derivatives with i+j <= 5 are ever used. Their recurrence can be written
 * explicitly in terms of y=q/s^2 and H'(y),...,H^(5)(y), which removes the generic coefficient
 * recurrence, three unused derivatives, and the logarithm that appears only in D00 and D10.
 */
static void utt4_ln_generator_derivatives(utt4_real s, utt4_real q,
    utt4_real D[UTT4_LN_NUM_GENERATOR_DERIVATIVES])
{
    const utt4_real inv_s = 1.0/s;
    const utt4_real inv_s2 = inv_s*inv_s;
    const utt4_real inv_s3 = inv_s2*inv_s;
    const utt4_real inv_s4 = inv_s2*inv_s2;
    const utt4_real inv_s5 = inv_s4*inv_s;
    const utt4_real inv_s6 = inv_s3*inv_s3;
    const utt4_real inv_s7 = inv_s6*inv_s;
    const utt4_real inv_s8 = inv_s4*inv_s4;
    const utt4_real inv_s9 = inv_s8*inv_s;

    const utt4_real y = q*inv_s2;
    const utt4_real y2 = y*y;
    const utt4_real y3 = y2*y;
    const utt4_real y4 = y2*y2;

    utt4_real h[6];
    utt4_ln_H_derivatives(y, h);
    const utt4_real h1 = h[1];
    const utt4_real h2 = h[2];
    const utt4_real h3 = h[3];
    const utt4_real h4 = h[4];
    const utt4_real h5 = h[5];

    UTT4_LN_D(D, 0, 2) = h2*inv_s3;
    UTT4_LN_D(D, 0, 3) = h3*inv_s5;
    UTT4_LN_D(D, 0, 4) = h4*inv_s7;
    UTT4_LN_D(D, 0, 5) = h5*inv_s9;

    UTT4_LN_D(D, 1, 1) = (-h1 - 2.0*h2*y)*inv_s2;
    UTT4_LN_D(D, 1, 2) = (-3.0*h2 - 2.0*h3*y)*inv_s4;
    UTT4_LN_D(D, 1, 3) = (-5.0*h3 - 2.0*h4*y)*inv_s6;
    UTT4_LN_D(D, 1, 4) = (-7.0*h4 - 2.0*h5*y)*inv_s8;

    UTT4_LN_D(D, 2, 0) = (2.0*y*(h1 + 2.0*h2*y) + 2.0)*inv_s;
    UTT4_LN_D(D, 2, 1) = 2.0*(h1 + 5.0*h2*y + 2.0*h3*y2)*inv_s3;
    UTT4_LN_D(D, 2, 2) = 2.0*(6.0*h2 + 9.0*h3*y + 2.0*h4*y2)*inv_s5;
    UTT4_LN_D(D, 2, 3) = 2.0*(15.0*h3 + 13.0*h4*y + 2.0*h5*y2)*inv_s7;

    UTT4_LN_D(D, 3, 0) = (-2.0*y*(3.0*h1 + 12.0*h2*y + 4.0*h3*y2) - 2.0)*inv_s2;
    UTT4_LN_D(D, 3, 1) = -2.0*(3.0*h1 + 27.0*h2*y + 24.0*h3*y2 + 4.0*h4*y3)*inv_s4;
    UTT4_LN_D(D, 3, 2) = -2.0*(30.0*h2 + 75.0*h3*y + 36.0*h4*y2 + 4.0*h5*y3)*inv_s6;

    UTT4_LN_D(D, 4, 0) = (4.0*y*(6.0*h1 + 39.0*h2*y + 28.0*h3*y2 + 4.0*h4*y3) + 4.0)*inv_s3;
    UTT4_LN_D(D, 4, 1) = 4.0*(6.0*h1 + 84.0*h2*y + 123.0*h3*y2 + 44.0*h4*y3 + 4.0*h5*y4)*inv_s5;

    UTT4_LN_D(D, 5, 0) = (-4.0*y*(30.0*h1 + 285.0*h2*y + 330.0*h3*y2 + 100.0*h4*y3 + 8.0*h5*y4) - 12.0)*inv_s4;
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
    utt4_real d2[UTT4_LN_NPAIR];
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
        utt4_real d2 = 0.0;
        for (int k = 0; k < 3; ++k) {
            const utt4_real dx = (utt4_real)pos[a][k] - (utt4_real)pos[b][k];
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
        rep->xy = d6_scale(d6_add(d6_sub(dad, dac), d6_sub(dbc, dbd)), 0.5);
        rep->xw = d6_scale(d6_sub(d6_sub(dad, dab), dbd), 0.5);
        rep->wy = d6_scale(d6_sub(d6_add(dbd, dcd), dbc), 0.5);
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
    utt4_real alpha;
    utt4_real beta;
    utt4_real gamma;
    utt4_real delta;
    utt4_real A;
    utt4_real B;

    utt4_real B2;
    utt4_real A2;
    utt4_real delta_gamma;
    utt4_real alpha_beta;
    utt4_real alpha_beta_delta_gamma;
    utt4_real sixty_alpha_beta_delta_gamma;
    utt4_real A2_delta_gamma;
    utt4_real alpha_beta_B2;
    utt4_real two_A;
    utt4_real two_A_B;
    utt4_real A_B;
    utt4_real two_B;
    utt4_real three_gamma;
    utt4_real three_delta;
} UTT4LnNodeParameters;


static inline UTT4LnNodeParameters utt4_ln_node_parameters(
    utt4_real alpha, utt4_real A, utt4_real gamma, utt4_real B)
{
    UTT4LnNodeParameters p;
    p.alpha = alpha;
    p.beta = 1.0 - alpha;
    p.gamma = gamma;
    p.delta = 1.0 - gamma;
    p.A = A;
    p.B = B;
    p.B2 = B*B;
    p.A2 = A*A;
    p.delta_gamma = p.delta*gamma;
    p.alpha_beta = alpha*p.beta;
    p.alpha_beta_delta_gamma = p.alpha_beta*p.delta_gamma;
    p.sixty_alpha_beta_delta_gamma = 60.0*p.alpha_beta_delta_gamma;
    p.A2_delta_gamma = p.A2*p.delta_gamma;
    p.alpha_beta_B2 = p.alpha_beta*p.B2;
    p.two_A = 2.0*A;
    p.two_A_B = B*p.two_A;
    p.A_B = A*B;
    p.two_B = 2.0*B;
    p.three_gamma = 3.0*gamma;
    p.three_delta = 3.0*p.delta;
    return p;
}


typedef struct {
    utt4_real value;
    utt4_real d_s;
    utt4_real d_q;
    utt4_real d_nm;
    utt4_real d_np;
    utt4_real d_pm;
    utt4_real d_inv_r;
    utt4_real d_inv_rho;
} UTT4LnNodePartials;


static inline UTT4LnNodePartials utt4_ln_contracted_node_partials(
    const utt4_real D[UTT4_LN_NUM_GENERATOR_DERIVATIVES],
    const UTT4LnNodeParameters *node, utt4_real nm, utt4_real npv, utt4_real pm,
    utt4_real q, utt4_real inv_r, utt4_real inv_rho)
{
    const utt4_real alpha = node->alpha;
    const utt4_real beta = node->beta;
    const utt4_real gamma = node->gamma;
    const utt4_real delta = node->delta;
    const utt4_real A = node->A;
    const utt4_real B = node->B;
    const utt4_real D02 = UTT4_LN_D(D, 0, 2);
    const utt4_real D03 = UTT4_LN_D(D, 0, 3);
    const utt4_real D04 = UTT4_LN_D(D, 0, 4);
    const utt4_real D05 = UTT4_LN_D(D, 0, 5);
    const utt4_real D11 = UTT4_LN_D(D, 1, 1);
    const utt4_real D12 = UTT4_LN_D(D, 1, 2);
    const utt4_real D13 = UTT4_LN_D(D, 1, 3);
    const utt4_real D14 = UTT4_LN_D(D, 1, 4);
    const utt4_real D20 = UTT4_LN_D(D, 2, 0);
    const utt4_real D21 = UTT4_LN_D(D, 2, 1);
    const utt4_real D22 = UTT4_LN_D(D, 2, 2);
    const utt4_real D23 = UTT4_LN_D(D, 2, 3);
    const utt4_real D30 = UTT4_LN_D(D, 3, 0);
    const utt4_real D31 = UTT4_LN_D(D, 3, 1);
    const utt4_real D32 = UTT4_LN_D(D, 3, 2);
    const utt4_real D40 = UTT4_LN_D(D, 4, 0);
    const utt4_real D41 = UTT4_LN_D(D, 4, 1);
    const utt4_real D50 = UTT4_LN_D(D, 5, 0);

    UTT4LnNodePartials out;
    const utt4_real t0 = A*inv_rho;
    const utt4_real t1 = B*inv_r;
    const utt4_real t2 = D20*t1;
    const utt4_real t3 = node->B2;
    const utt4_real t4 = D30*t3;
    const utt4_real t5 = A*inv_r;
    const utt4_real t6 = node->A2;
    const utt4_real t7 = D30*t6;
    const utt4_real t8 = B*inv_rho;
    const utt4_real t9 = node->delta_gamma;
    const utt4_real t10 = t5*t9;
    const utt4_real t11 = 4*D11;
    const utt4_real t12 = node->alpha_beta;
    const utt4_real t13 = t12*t8;
    const utt4_real t14 = node->alpha_beta_delta_gamma;
    const utt4_real t15 = node->sixty_alpha_beta_delta_gamma;
    const utt4_real t16 = node->A2_delta_gamma;
    const utt4_real t17 = 2*D21;
    const utt4_real t18 = node->alpha_beta_B2;
    const utt4_real t19 = D21*alpha;
    const utt4_real t20 = node->two_A;
    const utt4_real t21 = node->two_A_B;
    const utt4_real t22 = nm*t21;
    const utt4_real t23 = node->A_B;
    const utt4_real t24 = nm*t23;
    const utt4_real t25 = gamma*t24;
    const utt4_real t26 = node->two_B;
    const utt4_real t27 = npv*t0;
    const utt4_real t28 = t26*t27;
    const utt4_real t29 = D21*beta;
    const utt4_real t30 = delta*t24;
    const utt4_real t31 = D21*delta;
    const utt4_real t32 = pm*t1;
    const utt4_real t33 = t20*t32;
    const utt4_real t34 = D21*gamma;
    const utt4_real t35 = A*npv;
    const utt4_real t36 = t35*t9;
    const utt4_real t37 = alpha*t36;
    const utt4_real t38 = 20*D12;
    const utt4_real t39 = beta*t36;
    const utt4_real t40 = q*t10;
    const utt4_real t41 = 4*D12;
    const utt4_real t42 = B*pm;
    const utt4_real t43 = t12*t42;
    const utt4_real t44 = delta*t43;
    const utt4_real t45 = gamma*t43;
    const utt4_real t46 = q*t13;
    const utt4_real t47 = q*t9;
    const utt4_real t48 = 80*t12*t47;
    const utt4_real t49 = (nm*nm);
    const utt4_real t50 = t0*t49;
    const utt4_real t51 = alpha*t42;
    const utt4_real t52 = nm*t0;
    const utt4_real t53 = t51*t52;
    const utt4_real t54 = beta*t42;
    const utt4_real t55 = t52*t54;
    const utt4_real t56 = delta*t35;
    const utt4_real t57 = nm*t1;
    const utt4_real t58 = t56*t57;
    const utt4_real t59 = gamma*t35;
    const utt4_real t60 = t57*t59;
    const utt4_real t61 = 4*D22;
    const utt4_real t62 = t51*t56;
    const utt4_real t63 = alpha*q;
    const utt4_real t64 = t25*t63;
    const utt4_real t65 = beta*q;
    const utt4_real t66 = t30*t65;
    const utt4_real t67 = t54*t59;
    const utt4_real t68 = 8*q;
    const utt4_real t69 = D13*t68;
    const utt4_real t70 = nm*pm*t20;
    const utt4_real t71 = D31*t3*t70;
    const utt4_real t72 = (npv*npv);
    const utt4_real t73 = t10*t72;
    const utt4_real t74 = nm*npv*t26;
    const utt4_real t75 = D31*t6*t74;
    const utt4_real t76 = (pm*pm);
    const utt4_real t77 = t13*t76;
    const utt4_real t78 = (q*q);
    const utt4_real t79 = 16*t14*t78;
    const utt4_real t80 = t3*t49*t6;
    const utt4_real t81 = t4*t49;
    const utt4_real t82 = t49*t7;
    const utt4_real t83 = t16*t72;
    const utt4_real t84 = t18*t76;
    const utt4_real t85 = t49 + 1;
    const utt4_real t86 = 4*t10 + 4*t13;
    const utt4_real t87 = t51 - t54 + t56 - t59;
    const utt4_real t88 = t23*(-t0 + t1*t49 - t1 + t50);
    const utt4_real t89 = t68*(t37 - t39 + t44 - t45);
    const utt4_real t90 = -4*t62 + 4*t64 + 4*t66 - 4*t67 + 4*t83 + 4*t84;
    const utt4_real t91 = 20*t37 - 20*t39 + 4*t40 + 20*t44 - 20*t45 + 4*t46 - 4*t73 - 4*t77;
    const utt4_real t92 = B*t27;
    const utt4_real t93 = A*t32;
    const utt4_real t94 = node->three_gamma;
    const utt4_real t95 = node->three_delta;
    const utt4_real t96 = 2*alpha*t24*t94 - 2*alpha*t30 - 2*alpha*t92 + 2*beta*t24*t95 - 2*beta*t25 + 2*beta*t92 - 2*delta*t93 + 2*gamma*t93 + 2*t16 + 2*t18 + 2*t53 - 2*t55 + 2*t58 - 2*t60;
    const utt4_real t97 = 8*D13;
    const utt4_real t98 = D30*nm;
    const utt4_real t99 = inv_rho*pm;
    const utt4_real t100 = inv_r*npv;
    const utt4_real t101 = 2*D22;
    const utt4_real t102 = gamma*t101;
    const utt4_real t103 = delta*t101;
    const utt4_real t104 = 10*D12;
    const utt4_real t105 = t104*t9;
    const utt4_real t106 = t41*t9;
    const utt4_real t107 = 4*D13;
    const utt4_real t108 = t107*t9;
    const utt4_real t109 = t104*t12;
    const utt4_real t110 = D31*t24;
    const utt4_real t111 = t12*t41;
    const utt4_real t112 = q*t107*t12;
    const utt4_real t113 = D20*t8;
    const utt4_real t114 = 2*t42;
    const utt4_real t115 = D20*t5;
    const utt4_real t116 = 2*t35;
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
    const utt4_real alpha = node->alpha;
    const utt4_real beta = node->beta;
    const utt4_real gamma = node->gamma;
    const utt4_real delta = node->delta;
    const utt4_real A = node->A;
    const utt4_real B = node->B;

    const utt4_real dab = geo->d2[rep->pab];
    const utt4_real dcd = geo->d2[rep->pcd];
    const utt4_real dac = geo->d2[rep->pac];
    const utt4_real dad = geo->d2[rep->pad];
    const utt4_real dbc = geo->d2[rep->pbc];
    const utt4_real dbd = geo->d2[rep->pbd];

    const utt4_real q_ac = alpha*gamma;
    const utt4_real q_ad = alpha*delta;
    const utt4_real q_bc = beta*gamma;
    const utt4_real q_bd = beta*delta;
    const utt4_real q_ab = -alpha*beta;
    const utt4_real q_cd = -gamma*delta;

    utt4_real q = q_ac*dac + q_ad*dad + q_bc*dbc + q_bd*dbd
        + q_ab*dab + q_cd*dcd;
    const utt4_real s = A*rep->r.v + B*rep->rho.v;

    if (!(s > 0.0) || q < 0.0) {
        if (q > -1e-26*UTT4_FMAX(1.0, s*s))
            q = 0.0;
        else {
            Dual6 bad = d6_const(NAN);
            for (int k = 0; k < UTT4_LN_NPAIR; ++k)
                bad.g[k] = NAN;
            return bad;
        }
    }

    const utt4_real xp = rep->xw.v + alpha*dab - gamma*rep->xy.v;
    const utt4_real py = rep->wy.v + alpha*rep->xy.v - gamma*dcd;
    const utt4_real inv_r = rep->inv_r.v;
    const utt4_real inv_rho = rep->inv_rho.v;
    const utt4_real npv = xp*inv_r;
    const utt4_real pm = py*inv_rho;

    utt4_real D[UTT4_LN_NUM_GENERATOR_DERIVATIVES];
    utt4_ln_generator_derivatives(s, q, D);

    const UTT4LnNodePartials p = utt4_ln_contracted_node_partials(
        D, node, rep->nm.v, npv, pm, q, inv_r, inv_rho);

    utt4_real dq[UTT4_LN_NPAIR] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    dq[rep->pac] = q_ac;
    dq[rep->pad] = q_ad;
    dq[rep->pbc] = q_bc;
    dq[rep->pbd] = q_bd;
    dq[rep->pab] = q_ab;
    dq[rep->pcd] = q_cd;

    Dual6 out;
    out.v = p.value;
    for (int k = 0; k < UTT4_LN_NPAIR; ++k) {
        utt4_real dxp = rep->xw.g[k] - gamma*rep->xy.g[k];
        utt4_real dpy = rep->wy.g[k] + alpha*rep->xy.g[k];
        if (k == rep->pab)
            dxp += alpha;
        if (k == rep->pcd)
            dpy -= gamma;

        const utt4_real ds = A*rep->r.g[k] + B*rep->rho.g[k];
        const utt4_real dnp = dxp*inv_r + xp*rep->inv_r.g[k];
        const utt4_real dpm = dpy*inv_rho + py*rep->inv_rho.g[k];

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
static Dual6 utt4_ln_sum_node(const UTT4LnPairGeometry *geo, utt4_real alpha, utt4_real A,
    utt4_real gamma, utt4_real B)
{
    const UTT4LnNodeParameters node = utt4_ln_node_parameters(alpha, A, gamma, B);
    Dual6 sum = d6_zero();
    for (int p = 0; p < 6; ++p)
        sum = d6_add(sum, utt4_ln_ordered_node(geo, &geo->representative[p], &node));
    return d6_scale(sum, -4.0/UTT4_PI);
}


// ------------------------------------------------------------------------------------------------
// Gauss-Legendre quadrature
// ------------------------------------------------------------------------------------------------

typedef struct {
    utt4_real value;
    utt4_real d_d2[UTT4_LN_NPAIR];
    utt4_real grad[UTT4_LN_INTEGRAL_NUM_BODIES][3];
    utt4_real error_estimate;
    utt4_real max_component_error;
    long long node_evals;
} UTT4LnWorkResult;


/**
 * @brief Grading levels for one representative. 
 * 
 * Each level halves the panel nearest the endpoint, so `levels` bisections reach a width of about 
 * 2^-levels; enough to resolve a layer of width eps needs levels ~ log2(1/eps), plus one for 
 * margin. Balanced representatives get 0 and fall back to the plain single-panel rule.
 */
static int utt4_ln_grading_levels_from_geometry(const UTT4LnPairGeometry *geo,
    utt4_real *gmin, utt4_real *gmax)
{
    utt4_real big = 0.0, small = 0.0;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const utt4_real d = geo->distance[p].v;
        if (p == 0 || d > big) big = d;
        if (p == 0 || d < small) small = d;
    }
    if (gmin) *gmin = small;
    if (gmax) *gmax = big;

    if (!(big > 0.0) || !(small > 0.0) || !(big > small))
        return 0;

    const utt4_real ratio = big/small;
    if (ratio < UTT4_LN_GRADING_MIN_RATIO)
        return 0;

    int levels = 1;
    utt4_real reach = 2.0;
    while (reach < ratio && levels < UTT4_LN_MAX_GRADING_LEVELS) {
        reach *= 2.0;
        ++levels;
    }

    return levels;
}


/**
 * @brief Fixed-order quadrature of all six ordered representatives.
 *
 * Each representative is integrated on its own grid rather than on one grid shared by all six.
 * A representative whose two pair separations are comparable gets a plain base_order rule on both
 * axes; one carrying a large separation ratio gets graded panels on both. Grading only the axis
 * conjugate to the larger separation is enough for a distant fourth body, but not for a tight
 * pair inside a wide system, where the small separation also sharpens the other axis through the
 * 1/r factors, so both are graded whenever grading applies at all.
 *
 * Sharing one grid across the six would force every representative onto the worst one's
 * resolution. For a hardening binary that is one pair split out of three, and the other two need
 * nothing beyond the plain rule.
 */
static UTT4LnWorkResult utt4_ln_fixed_quadrature(const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3],
    int base_order, int use_grading, int use_openmp)
{
#ifndef _OPENMP
    (void)use_openmp;
#endif
    UTT4LnPairGeometry geo;
    utt4_ln_build_pair_geometry(pos, &geo);

    UTT4LnWorkResult out;
    memset(&out, 0, sizeof(out));

    Dual6 total = d6_zero();
    long long node_evals = 0;

    const UTT4LnQuadratureRule *smooth = utt4_ln_quadrature_rule_get(base_order, 0);
    if (smooth == NULL) {
        out.value = NAN;
        return out;
    }

    utt4_real gmin = 0.0, gmax = 0.0;
    const int levels = use_grading
        ? utt4_ln_grading_levels_from_geometry(&geo, &gmin, &gmax) : 0;

    int panel_order = base_order/UTT4_LN_PANEL_ORDER_DIVISOR;
    if (panel_order < UTT4_LN_PANEL_ORDER_MIN)
        panel_order = UTT4_LN_PANEL_ORDER_MIN;
    if (panel_order > base_order)
        panel_order = base_order;

    const UTT4LnQuadratureRule *graded =
        (levels > 0) ? utt4_ln_quadrature_rule_get(panel_order, levels) : smooth;
    if (graded == NULL) {
        out.value = NAN;
        return out;
    }

    /*
     * Resolve each representative's two axis rules, then group the representatives that share a
     * rule pair. The node parameters depend only on the node, not on the representative, so a
     * group evaluates them once per node and reuses them across its members; without the grouping
     * a plain evaluation would rebuild them six times per node instead of once.
     */
    const UTT4LnQuadratureRule *group_theta[6];
    const UTT4LnQuadratureRule *group_phi[6];
    int group_member[6][6];
    int group_size[6];
    int num_groups = 0;

    for (int p = 0; p < 6; ++p) {
        const UTT4LnRepresentativeGeometry *rep = &geo.representative[p];

        /*
         * A representative whose own two separations already span the quadruple's full scale
         * range carries the disparity entirely in s = A r + B rho, so grading the axis conjugate
         * to the larger separation resolves it. When the extreme separations instead sit among
         * the representative's cross distances, they reach the integrand through q and sharpen
         * both axes, and grading only one leaves the refinement loop climbing.
         */
        const utt4_real rmin = (rep->r.v < rep->rho.v) ? rep->r.v : rep->rho.v;
        const utt4_real rmax = (rep->r.v < rep->rho.v) ? rep->rho.v : rep->r.v;
        const int spans_range = (levels > 0) && (rmin <= 2.0*gmin) && (2.0*rmax >= gmax);
        const int layer_in_theta = (rep->r.v >= rep->rho.v);

        const UTT4LnQuadratureRule *qt = (!spans_range || layer_in_theta) ? graded : smooth;
        const UTT4LnQuadratureRule *qp = (!spans_range || !layer_in_theta) ? graded : smooth;
        node_evals += (long long)qt->n*(long long)qp->n;

        int g = 0;
        while (g < num_groups && (group_theta[g] != qt || group_phi[g] != qp))
            ++g;
        if (g == num_groups) {
            group_theta[g] = qt;
            group_phi[g] = qp;
            group_size[g] = 0;
            ++num_groups;
        }
        group_member[g][group_size[g]++] = p;
    }

#ifdef _OPENMP
    if (use_openmp && omp_get_max_threads() > 1 && !omp_in_parallel()) {
        #pragma omp parallel
        {
            Dual6 local = d6_zero();
            for (int g = 0; g < num_groups; ++g) {
                const UTT4LnQuadratureRule *qt = group_theta[g];
                const UTT4LnQuadratureRule *qp = group_phi[g];

                #pragma omp for collapse(2) schedule(static)
                for (int i = 0; i < qt->n; ++i) {
                    for (int j = 0; j < qp->n; ++j) {
                        const UTT4LnNodeParameters node = utt4_ln_node_parameters(
                            qt->mix[i], qt->half_sin_theta[i],
                            qp->mix[j], qp->half_sin_theta[j]);
                        const utt4_real w = qt->weight[i]*qp->weight[j];

                        // Sum the group's representatives before weighting: one scaling per
                        // node instead of one per representative.
                        Dual6 node_sum = d6_zero();
                        for (int m = 0; m < group_size[g]; ++m)
                            node_sum = d6_add(node_sum, utt4_ln_ordered_node(&geo,
                                &geo.representative[group_member[g][m]], &node));
                        local = d6_add(local, d6_scale(node_sum, w));
                    }
                }
            }
            #pragma omp critical(prod_force_sum)
            total = d6_add(total, local);
        }
    } else
#endif
    {
        for (int g = 0; g < num_groups; ++g) {
            const UTT4LnQuadratureRule *qt = group_theta[g];
            const UTT4LnQuadratureRule *qp = group_phi[g];

            for (int i = 0; i < qt->n; ++i)
                for (int j = 0; j < qp->n; ++j) {
                    const UTT4LnNodeParameters node = utt4_ln_node_parameters(
                        qt->mix[i], qt->half_sin_theta[i],
                        qp->mix[j], qp->half_sin_theta[j]);
                    const utt4_real w = qt->weight[i]*qp->weight[j];

                    Dual6 node_sum = d6_zero();
                    for (int m = 0; m < group_size[g]; ++m)
                        node_sum = d6_add(node_sum, utt4_ln_ordered_node(&geo,
                            &geo.representative[group_member[g][m]], &node));
                    total = d6_add(total, d6_scale(node_sum, w));
                }
        }
    }

    // The factor four reconstructs the full ordered sum from the six inequivalent
    // representatives, using I_ab;cd = I_ba;dc = I_cd;ab = I_dc;ba. The 1/pi is the
    // normalization of the logarithmic integral.
    total = d6_scale(total, -4.0/UTT4_PI);

    out.value = total.v;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p)
        out.d_d2[p] = total.g[p];
    out.node_evals = node_evals;

    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const int a = utt4_ln_pair_body_i[p], b = utt4_ln_pair_body_j[p];
        for (int k = 0; k < 3; ++k) {
            const utt4_real dx = (utt4_real)pos[a][k] - (utt4_real)pos[b][k];
            const utt4_real c = 2.0*out.d_d2[p]*dx;
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
static const utt4_real gk15_x[15] = {
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

static const utt4_real gk15_wk[15] = {
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

static const utt4_real gk15_wg[15] = {
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
    utt4_real rel_tol;
    int max_depth;
    int parallel_root;
    utt4_real length2_scale;
} UTT4LnAdaptiveContext;


typedef struct {
    Dual6 result;
    utt4_real abs_error;
    utt4_real max_component_error;
    long long node_evals;
} UTT4LnAdaptiveAccum;


static UTT4LnCellEstimate utt4_ln_cell_estimate(const UTT4LnAdaptiveContext *ctx,
    utt4_real ta, utt4_real tb, utt4_real pa, utt4_real pb)
{
    UTT4LnCellEstimate out;
    out.kronrod = d6_zero();
    out.gauss = d6_zero();
    out.node_evals = 225;
    const utt4_real tm = 0.5*(ta + tb), th = 0.5*(tb - ta);
    const utt4_real pm = 0.5*(pa + pb), ph = 0.5*(pb - pa);

    utt4_real alpha[15], A[15], wk_t[15], wg_t[15];
    utt4_real gamma[15], B[15], wk_p[15], wg_p[15];
    for (int i = 0; i < 15; ++i) {
        const utt4_real theta = tm + th*gk15_x[i];
        const utt4_real phi = pm + ph*gk15_x[i];
        alpha[i] = 0.5*(1.0 + UTT4_COS(theta));
        A[i] = 0.5*UTT4_SIN(theta);
        wk_t[i] = th*gk15_wk[i];
        wg_t[i] = th*gk15_wg[i];
        gamma[i] = 0.5*(1.0 + UTT4_COS(phi));
        B[i] = 0.5*UTT4_SIN(phi);
        wk_p[i] = ph*gk15_wk[i];
        wg_p[i] = ph*gk15_wg[i];
    }

    Dual6 ksum = d6_zero(), gsum = d6_zero();
    for (int i = 0; i < 15; ++i) {
        for (int j = 0; j < 15; ++j) {
            Dual6 f = utt4_ln_sum_node(&ctx->geo, alpha[i], A[i], gamma[j], B[j]);
            ksum = d6_add(ksum, d6_scale(f, wk_t[i]*wk_p[j]));
            if (wg_t[i] != 0.0 && wg_p[j] != 0.0)
                gsum = d6_add(gsum, d6_scale(f, wg_t[i]*wg_p[j]));
        }
    }
    out.kronrod = ksum;
    out.gauss = gsum;
    return out;
}


static int utt4_ln_cell_accept(const UTT4LnAdaptiveContext *ctx, const UTT4LnCellEstimate *e,
    utt4_real local_abs, utt4_real *value_err, utt4_real *max_comp_err)
{
    const utt4_real ev = UTT4_FABS(e->kronrod.v - e->gauss.v);
    utt4_real maxe = ev;
    int ok = (ev <= local_abs + ctx->rel_tol*UTT4_FABS(e->kronrod.v));
    const utt4_real grad_abs = local_abs/UTT4_FMAX(ctx->length2_scale, 1e-30);
    for (int k = 0; k < UTT4_LN_NPAIR; ++k) {
        const utt4_real eg = UTT4_FABS(e->kronrod.g[k] - e->gauss.g[k]);
        if (eg > maxe) maxe = eg;
        if (eg > grad_abs + ctx->rel_tol*UTT4_FABS(e->kronrod.g[k])) ok = 0;
    }
    *value_err=ev;
    *max_comp_err=maxe;
    return ok;
}

static UTT4LnAdaptiveAccum utt4_ln_adaptive_cell(const UTT4LnAdaptiveContext *ctx, utt4_real ta,
    utt4_real tb, utt4_real pa, utt4_real pb, utt4_real local_abs, int depth)
{
    UTT4LnAdaptiveAccum acc;
    memset(&acc, 0, sizeof(acc));
    UTT4LnCellEstimate e = utt4_ln_cell_estimate(ctx, ta, tb, pa, pb);
    acc.node_evals = e.node_evals;
    utt4_real verr = 0.0, maxerr = 0.0;
    if (utt4_ln_cell_accept(ctx, &e, local_abs, &verr, &maxerr) || depth >= ctx->max_depth) {
        acc.result = e.kronrod;
        acc.abs_error = verr;
        acc.max_component_error = maxerr;
        return acc;
    }

    const utt4_real tm = 0.5*(ta + tb), pm = 0.5*(pa + pb);
    const utt4_real bounds[4][4] = {
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
                bounds[c][3], local_abs*0.25, depth + 1);
        }
    } else
#endif
    {
        for (int c = 0; c < 4; ++c)
            child[c] = utt4_ln_adaptive_cell(ctx, bounds[c][0], bounds[c][1], bounds[c][2],
                bounds[c][3], local_abs*0.25, depth + 1);
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
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3], utt4_real rel_tol, utt4_real abs_tol,
    int max_depth, int parallel_root)
{
    UTT4LnAdaptiveContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    utt4_ln_build_pair_geometry(pos, &ctx.geo);
    ctx.rel_tol = rel_tol;
    ctx.max_depth = max_depth;
    ctx.parallel_root = parallel_root;
    utt4_real maxd2 = 0.0;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) if (ctx.geo.d2[p] > maxd2) maxd2 = ctx.geo.d2[p];
    ctx.length2_scale = maxd2;

    UTT4LnAdaptiveAccum a = utt4_ln_adaptive_cell(&ctx, 0.0, UTT4_PI, 0.0, UTT4_PI, abs_tol, 0);
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
            const utt4_real dx = (utt4_real)pos[i][k] - (utt4_real)pos[j][k];
            const utt4_real contribution = 2.0*out.d_d2[p]*dx;
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
 * @brief Largest squared pair separation of a quadruple. 
 * 
 * Both tolerance checks below need it to turn an absolute tolerance on the value into one on the
 * squared-distance derivatives.
 */
static utt4_real utt4_ln_max_pair_d2(const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3])
{
    utt4_real maxd2 = 0.0;
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const int i = utt4_ln_pair_body_i[p], j = utt4_ln_pair_body_j[p];
        utt4_real d2 = 0.0;
        for (int k = 0; k < 3; ++k) {
            const utt4_real dx = (utt4_real)pos[i][k] - (utt4_real)pos[j][k];
            d2 += dx*dx;
        }
        if (d2 > maxd2) maxd2 = d2;
    }
    return maxd2;
}


/**
 * @brief Compare two fixed-order results using the value and all six squared-distance derivatives.
 *
 * A refinement is requested when any component exceeds its absolute-plus-relative tolerance.
 */
static int utt4_ln_needs_refinement(const UTT4LnWorkResult *low, const UTT4LnWorkResult *high,
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3], utt4_real rtol, utt4_real atol,
    utt4_real *value_error, utt4_real *max_component_error, utt4_real *value_rel_est,
    utt4_real *max_component_rel_est, utt4_real *worst_tol_ratio)
{
    const utt4_real ev = UTT4_FABS(high->value-low->value);
    const utt4_real value_scale = UTT4_FMAX(UTT4_FABS(high->value), 1e-300);
    const utt4_real vr = ev/value_scale;
    const utt4_real value_allowed = atol + rtol*UTT4_FABS(high->value);
    utt4_real worst = (value_allowed > 0.0) ? ev/value_allowed : (ev == 0.0 ? 0.0 : UTT4_HUGE);
    utt4_real maxe = ev, maxrel = vr;
    int refine = worst > 1.0;
    const utt4_real datol = atol/UTT4_FMAX(utt4_ln_max_pair_d2(pos), 1e-30);
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const utt4_real e = UTT4_FABS(high->d_d2[p] - low->d_d2[p]);
        const utt4_real scale = UTT4_FMAX(UTT4_FABS(high->d_d2[p]), 1e-300);
        const utt4_real rel = e/scale;
        const utt4_real allowed = datol + rtol*UTT4_FABS(high->d_d2[p]);
        const utt4_real ratio = (allowed > 0.0) ? e/allowed : (e == 0.0 ? 0.0 : UTT4_HUGE);
        if (e > maxe)maxe = e;
        if (rel > maxrel)maxrel = rel;
        if (ratio > worst)worst = ratio;
        if (ratio > 1.0)refine = 1;
    }
    if (value_error)*value_error = ev;
    if (max_component_error)*max_component_error = maxe;
    if (value_rel_est)*value_rel_est = vr;
    if (max_component_rel_est)*max_component_rel_est = maxrel;
    if (worst_tol_ratio)*worst_tol_ratio = worst;
    return refine;
}


// Thresholde for the lower comparison rule before its order is suggested for the next evaluation.
#define UTT4_LN_ORDER_DESCEND_MARGIN ((utt4_real)0.05)


static int integral_round_up_even_order(utt4_real x)
{
    int n = (int)UTT4_CEIL(x);
    if (n < 4)n = 4;
    if (n & 1)++n;
    return n;
}


static int integral_initial_high_order(utt4_real rtol, int min_order, int max_order)
{
    utt4_real digits = -UTT4_LOG10(rtol);
    if (digits < 1.0)digits = 1.0;
    int high = integral_round_up_even_order(2.7*digits);
    if (high < min_order + 2)high = min_order + 2;
    if (high > max_order)high = max_order;
    return high;
}


static int integral_previous_order(int high, int min_order)
{
    int low = integral_round_up_even_order(0.88*(utt4_real)high);
    if (low >= high)low = high - 2;
    if (low < min_order)low = min_order;
    if (low >= high)low = high - 2;
    return low;
}


static int integral_next_order(int high, int max_order)
{
    int next = (high < 64) ? high + 4 : integral_round_up_even_order(1.20*(utt4_real)high);
    if (next <= high)next = high + 2;
    if (next > max_order)next = max_order;
    return next;
}


static int utt4_ln_adaptive_meets_target(const UTT4LnWorkResult *r,
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3], utt4_real rtol, utt4_real atol,
    utt4_real *worst_ratio)
{
    const utt4_real value_allowed = atol + rtol*UTT4_FABS(r->value);
    utt4_real worst = (value_allowed > 0.0) ? r->error_estimate/value_allowed : 0.0;
    const utt4_real datol = atol/UTT4_FMAX(utt4_ln_max_pair_d2(pos), 1e-30);
    for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
        const utt4_real allowed = datol + rtol*UTT4_FABS(r->d_d2[p]);
        const utt4_real ratio = (allowed > 0.0) ? r->max_component_error/allowed : 0.0;
        if (ratio > worst)worst = ratio;
    }
    if (worst_ratio)*worst_ratio = worst;
    return worst <= 1.0;
}


// ------------------------------------------------------------------------------------------------
// Public UTT4 logarithmic-integral evaluator
// ------------------------------------------------------------------------------------------------

/**
 * @brief One fixed-order ladder in a single quadrature mode. 
 * 
 * Starts from *high_order, compares against the next rule down, and climbs until the pair agrees
 * or max_order is reached. Returns the converged result and reports through *refined whether the
 * tolerance was actually met.
 */
static int utt4_ln_run_ladder(const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3],
    const NumericalIntegralSettings *settings, int use_grading, int max_order,
    int *high_order, int *low_order,
    UTT4LnWorkResult *result, long long *node_evals, utt4_real *value_error,
    utt4_real *max_component_error, utt4_real *value_rel, utt4_real *max_component_rel,
    utt4_real *worst, int *refine)
{
    if (*high_order > max_order)
        *high_order = max_order;

    *low_order = integral_previous_order(*high_order, settings->min_order);
    if (*low_order < 4 || *high_order <= *low_order) {
        *low_order = settings->min_order;
        *high_order = *low_order + 2;
        if (*high_order > max_order)
            *high_order = max_order;
    }
    if (*high_order <= *low_order)
        return -1;

    UTT4LnWorkResult low = utt4_ln_fixed_quadrature(pos, *low_order, use_grading,
        settings->use_openmp);
    UTT4LnWorkResult high = utt4_ln_fixed_quadrature(pos, *high_order, use_grading,
        settings->use_openmp);
    if (isnan((double)low.value) || isnan((double)high.value))
        return -1;
    *node_evals += low.node_evals + high.node_evals;

    *refine = utt4_ln_needs_refinement(&low, &high, pos, settings->rel_tol, settings->abs_tol,
        value_error, max_component_error, value_rel, max_component_rel, worst);

    while (*refine && *high_order < max_order) {
        *low_order = *high_order;
        low = high;
        *high_order = integral_next_order(*high_order, max_order);
        if (*high_order <= *low_order)
            break;
        high = utt4_ln_fixed_quadrature(pos, *high_order, use_grading, settings->use_openmp);
        if (isnan((double)high.value))
            return -1;
        *node_evals += high.node_evals;
        *refine = utt4_ln_needs_refinement(&low, &high, pos, settings->rel_tol, settings->abs_tol,
            value_error, max_component_error, value_rel, max_component_rel, worst);
    }

    *result = high;
    return 0;
}


/**
 * @brief Evaluate the logarithmic four-body integral entering UTT4 and all position derivatives.
 *
 * A plain tensor-product rule is tried first, since it is the cheapest thing that works for the
 * great majority of quadruples. If its ladder runs past UTT4_LN_PLAIN_SWITCH_ORDER without
 * meeting the tolerance, the geometry has a scale disparity a global polynomial cannot resolve,
 * and the evaluation restarts on graded panels. Choosing between the two by trying rather than by
 * predicting matters: the geometric indicators that separate a hardening binary from a
 * binary-binary encounter disagree with each other, and guessing wrong costs more than probing.
 *
 * The mode is reported back through suggested_order as a sign, so a caller with order memory pays
 * the probe once per quadruple rather than on every evaluation.
 *
 * @param[in]  pos       Four local particle positions in three dimensions.
 * @param[in]  settings  Accuracy, order-refinement, adaptive, and OpenMP settings.
 * @param[out] out       Integral value, Cartesian derivatives, and quadrature diagnostics.
 * @return 0 on success, nonzero on invalid input or quadrature setup failure.
 */
int utt4_ln_integral_evaluate(const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][3],
    const NumericalIntegralSettings *settings, UTT4LnIntegralResult *out)
{
    if (!pos || !settings || !out) return -1;
    memset(out, 0, sizeof(*out));
    if (!(settings->rel_tol > 0.0) || !(settings->abs_tol > 0.0) ||
       settings->min_order < 4 || settings->max_order < settings->min_order ||
       settings->max_depth < 0) return -1;

    // A negative remembered order means that quadruple already needed graded panels
    int use_grading = (settings->start_order < 0);
    const int remembered = (settings->start_order < 0)
        ? -settings->start_order : settings->start_order;

    int high_order;
    if (remembered > 0) {
        high_order = integral_round_up_even_order(remembered);
        if (high_order < settings->min_order) high_order = settings->min_order;
        if (high_order > settings->max_order) high_order = settings->max_order;
    } else {
        high_order = integral_initial_high_order(settings->rel_tol, settings->min_order,
            settings->max_order);
    }

    if (remembered > 0 && !settings->verify) {
        const UTT4LnWorkResult fast = utt4_ln_fixed_quadrature(pos, high_order, use_grading,
            settings->use_openmp);
        if (isnan((double)fast.value))
            return -1;

        out->value = fast.value;
        for (int b = 0; b < UTT4_LN_INTEGRAL_NUM_BODIES; ++b)
            for (int k = 0; k < 3; ++k)
                out->grad[b][k] = fast.grad[b][k];

        out->diagnostics.node_evals = fast.node_evals;
        out->diagnostics.low_order = high_order;
        out->diagnostics.high_order = high_order;
        out->diagnostics.suggested_order = use_grading ? -high_order : high_order;
        out->diagnostics.target_met = 1;
        return 0;
    }

    utt4_real value_error = 0.0, max_component_error = 0.0;
    utt4_real value_rel = 0.0, max_component_rel = 0.0, worst = 0.0;
    int low_order = 0, refine = 1;
    long long node_evals = 0;
    UTT4LnWorkResult final;

    //The plain ladder is capped rather than run to max_order
    int plain_cap = settings->max_order;
    if (!use_grading && plain_cap > UTT4_LN_PLAIN_SWITCH_ORDER)
        plain_cap = UTT4_LN_PLAIN_SWITCH_ORDER;

    if (utt4_ln_run_ladder(pos, settings, use_grading, plain_cap, &high_order, &low_order, &final,
            &node_evals, &value_error, &max_component_error, &value_rel, &max_component_rel,
            &worst, &refine) != 0)
        return -1;

    if (refine && !use_grading) {
        use_grading = 1;
        high_order = integral_initial_high_order(settings->rel_tol, settings->min_order,
            settings->max_order);
        if (utt4_ln_run_ladder(pos, settings, use_grading, settings->max_order, &high_order,
                &low_order, &final, &node_evals, &value_error, &max_component_error, &value_rel,
                &max_component_rel, &worst, &refine) != 0)
            return -1;
    }

    out->diagnostics.node_evals = node_evals;
    out->diagnostics.low_order = low_order;
    out->diagnostics.high_order = high_order;
    out->diagnostics.target_met = !refine;
    out->diagnostics.value_error = value_error;
    out->diagnostics.max_derivative_error = max_component_error;
    out->diagnostics.value_rel_est = value_rel;
    out->diagnostics.max_derivative_rel_est = max_component_rel;
    out->diagnostics.worst_tolerance_ratio = worst;

    /*
     * Order to start from next time this quadruple is evaluated, carrying the mode in its sign.
     * When the comparison sweep showed the lower rule was already well inside tolerance, suggest
     * that lower order: the comparison is direct evidence it sufficed, so a caller feeding this
     * back walks down to the cheapest adequate rule instead of paying the conservative initial
     * guess forever. Nothing is suggested below the current order while a refinement is still
     * outstanding.
     */
    int suggested = high_order;
    if (!refine && low_order >= settings->min_order && worst <= UTT4_LN_ORDER_DESCEND_MARGIN)
        suggested = low_order;
    out->diagnostics.suggested_order = use_grading ? -suggested : suggested;

    if (refine && settings->adaptive) {
        final = utt4_ln_adaptive_quadrature(pos, settings->rel_tol, settings->abs_tol,
            settings->max_depth, settings->use_openmp);
        out->diagnostics.adaptive_used = 1;
        out->diagnostics.node_evals += final.node_evals;
        out->diagnostics.value_error = final.error_estimate;
        out->diagnostics.max_derivative_error = final.max_component_error;
        out->diagnostics.value_rel_est =
            final.error_estimate/UTT4_FMAX(UTT4_FABS(final.value), 1e-300);
        utt4_real maxrel = 0.0;
        for (int p = 0; p < UTT4_LN_NPAIR; ++p) {
            const utt4_real rel = final.max_component_error/
                UTT4_FMAX(UTT4_FABS(final.d_d2[p]), 1e-300);
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
