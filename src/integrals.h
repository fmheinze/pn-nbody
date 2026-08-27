#ifndef INTEGRALS_H
#define INTEGRALS_H


// ------------------------------------------------------------------------------------------------
// Working precision of the numerical-integral kernels
// ------------------------------------------------------------------------------------------------

/*
 * The quadrature kernels are written against utt4_real, which defaults to double. Double reaches
 * a relative accuracy of roughly 1e-13 on the logarithmic four-body integral, comfortably inside
 * the tolerances the evaluator is asked for, and is about three times faster than the extended
 * type on x86-64.
 *
 * Build with -DUTT4_USE_LONG_DOUBLE=1 to restore the extended-precision kernel, for example for
 * convergence studies or to check a suspicious result. Note that the cost of that build is highly
 * platform dependent: long double is 80-bit x87 on x86-64, plain 64-bit double on Apple Silicon,
 * and software-emulated binary128 on AArch64 Linux.
 */
#ifndef UTT4_USE_LONG_DOUBLE
#define UTT4_USE_LONG_DOUBLE 0
#endif

#if UTT4_USE_LONG_DOUBLE
typedef long double utt4_real;
#else
typedef double utt4_real;
#endif


// ------------------------------------------------------------------------------------------------
// Generic numerical-integral configuration
// ------------------------------------------------------------------------------------------------

// Accuracy and refinement settings shared by numerical integral evaluators
typedef struct {
    utt4_real rel_tol;
    utt4_real abs_tol;
    int min_order;
    int max_order;
    int adaptive;
    int max_depth;
    int use_openmp;

    /*
     * Optional per-caller order memory. start_order overrides the tolerance-derived initial
     * Gauss-Legendre order; 0 selects the built-in heuristic. When start_order is positive and
     * verify is zero the evaluator runs a single sweep at that order and skips the comparison
     * sweep, which is what makes a remembered order cheaper than a fresh one. Callers that use
     * this must re-verify periodically (verify nonzero) so the order can be corrected, and must
     * treat target_met from an unverified evaluation as inherited rather than established.
     */
    int start_order;
    int verify;
} NumericalIntegralSettings;


// Diagnostics returned by an automatic or adaptive integral evaluation
typedef struct {
    utt4_real value_error;
    utt4_real max_derivative_error;
    utt4_real value_rel_est;
    utt4_real max_derivative_rel_est;
    utt4_real worst_tolerance_ratio;
    long long node_evals;
    int low_order;
    int high_order;
    int adaptive_used;
    int target_met;

    /*
     * Order a caller should start from next time this quadruple is evaluated at a nearby
     * geometry. It drops below high_order when the comparison sweep showed the lower rule was
     * already well inside tolerance, so a caller that feeds it back converges onto the cheapest
     * sufficient order instead of paying the conservative initial guess forever.
     */
    int suggested_order;
} NumericalIntegralDiagnostics;


// ------------------------------------------------------------------------------------------------
// UTT4 logarithmic integral
// ------------------------------------------------------------------------------------------------

#define UTT4_LN_INTEGRAL_NUM_BODIES 4
#define UTT4_LN_INTEGRAL_NUM_DIM 3


// Result and gradient of the logarithmic four-body integral evaluation entering UTT4
typedef struct {
    utt4_real value;
    utt4_real grad[UTT4_LN_INTEGRAL_NUM_BODIES][UTT4_LN_INTEGRAL_NUM_DIM];
    NumericalIntegralDiagnostics diagnostics;
} UTT4LnIntegralResult;


int utt4_ln_integral_evaluate(
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][UTT4_LN_INTEGRAL_NUM_DIM],
    const NumericalIntegralSettings *settings, UTT4LnIntegralResult *result);


#endif
