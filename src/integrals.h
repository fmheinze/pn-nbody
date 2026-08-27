#ifndef INTEGRALS_H
#define INTEGRALS_H


// ------------------------------------------------------------------------------------------------
// Generic numerical-integral configuration
// ------------------------------------------------------------------------------------------------

// Accuracy and refinement settings shared by numerical integral evaluators
typedef struct {
    long double rel_tol;
    long double abs_tol;
    int min_order;
    int max_order;
    int adaptive;
    int max_depth;
    int use_openmp;
} NumericalIntegralSettings;


// Diagnostics returned by an automatic or adaptive integral evaluation
typedef struct {
    long double value_error;
    long double max_derivative_error;
    long double value_rel_est;
    long double max_derivative_rel_est;
    long double worst_tolerance_ratio;
    long long node_evals;
    int low_order;
    int high_order;
    int adaptive_used;
    int target_met;
} NumericalIntegralDiagnostics;


// ------------------------------------------------------------------------------------------------
// UTT4 logarithmic integral
// ------------------------------------------------------------------------------------------------

#define UTT4_LN_INTEGRAL_NUM_BODIES 4
#define UTT4_LN_INTEGRAL_NUM_DIM 3


// Result and gradient of the logarithmic four-body integral evaluation entering UTT4
typedef struct {
    long double value;
    long double grad[UTT4_LN_INTEGRAL_NUM_BODIES][UTT4_LN_INTEGRAL_NUM_DIM];
    NumericalIntegralDiagnostics diagnostics;
} UTT4LnIntegralResult;


int utt4_ln_integral_evaluate(
    const double pos[UTT4_LN_INTEGRAL_NUM_BODIES][UTT4_LN_INTEGRAL_NUM_DIM],
    const NumericalIntegralSettings *settings, UTT4LnIntegralResult *result);


#endif
