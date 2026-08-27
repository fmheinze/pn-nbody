/**
 * @file utils.c
 * @brief Various small useful functions.
 *
 * Various small useful functions, e.g. for allocation/freeing, mathematics, printing or
 * system calls.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <limits.h>
#include <libgen.h>
#include <errno.h>
#include "utils.h"
#include "eom.h"
#include "cache.h"

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif


// ------------------------------------------------------------------------------------------------
// Allocate / free
// ------------------------------------------------------------------------------------------------

void allocate_vector(double** ptr, int num_elements)
{
    *ptr = (double *)malloc(num_elements * sizeof(double));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}


void allocate_2d_array(double*** ptr, int num_vectors, int num_elements)
{
    *ptr = (double **)malloc(num_vectors * sizeof(double *));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for 2D array\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < num_vectors; i++) {
        allocate_vector(&((*ptr)[i]), num_elements);
    }
}


void allocate_3d_array(double**** ptr, int num_arrays, int num_vectors, int num_elements)
{
    // Allocate memory for the array of 2D arrays
    *ptr = (double ***)malloc(num_arrays * sizeof(double **));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for 3D array\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each 2D array
    for (int i = 0; i < num_arrays; i++) {
        allocate_2d_array(&((*ptr)[i]), num_vectors, num_elements);
    }
}


void allocate_4d_array(double***** ptr, int num_3d_arrays, int num_arrays, int num_vectors,
    int num_elements)
{
    // Allocate memory for the array of 2D arrays
    *ptr = (double ****)malloc(num_3d_arrays * sizeof(double ***));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for 4D array\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each 3D array
    for (int i = 0; i < num_3d_arrays; i++) {
        allocate_3d_array(&((*ptr)[i]), num_arrays, num_vectors, num_elements);
    }
}


void free_vector(double* ptr)
{
    if (ptr != NULL)
        free(ptr);
}


void free_2d_array(double** ptr, int num_vectors)
{
    if (ptr != NULL) {
        for (int i = 0; i < num_vectors; i++)
            free_vector(ptr[i]);
        free(ptr);
    }
}


void free_3d_array(double*** ptr, int num_arrays, int num_vectors)
{
    if (ptr != NULL) {
        for (int i = 0; i < num_arrays; i++)
            free_2d_array(ptr[i], num_vectors);
        free(ptr);
    }
}


void free_4d_array(double**** ptr, int num_3d_arrays, int num_arrays, int num_vectors)
{
    for (int i = 0; i < num_3d_arrays; i++)
        free_3d_array(ptr[i], num_arrays, num_vectors);
    free(ptr);
}


void free_ode_params(struct ode_params* params)
{
    if (params == NULL)
        return;

    pair_cache_destroy(params->pair_cache);
    params->pair_cache = NULL;

    utt4_cache_destroy(params->utt4_cache);
    params->utt4_cache = NULL;

    free_vector(params->masses);
    params->masses = NULL;

    free(params->pn_terms);
    params->pn_terms = NULL;

    free(params->active);
    params->active = NULL;

    free(params->body_id);
    params->body_id = NULL;

    free(params->generation);
    params->generation = NULL;

    free(params->remnant_prescription);
    params->remnant_prescription = NULL;

    params->num_active = 0;
    params->next_body_id = 0;
}


// ------------------------------------------------------------------------------------------------
// Math utils
// ------------------------------------------------------------------------------------------------

double dot_product(const double *a, const double *b, int dim)
{
    double result = 0;
    for (int i = 0; i < dim; i++)
        result += a[i] * b[i];
    return result;
}


complex double dot_product_c(const complex double *a, const complex double *b, int dim)
{
    complex double result = 0;
    for (int i = 0; i < dim; i++)
        result += a[i] * b[i];
    return result;
}


double norm(const double *v, int dim)
{
    return sqrt(dot_product(v, v, dim));
}


complex double norm_c(const complex double *v, int dim)
{
    return csqrt(dot_product_c(v, v, dim));
}


void normalize(double v[3], double result[3])
{
    double mag = norm(v, 3);
    if (!isfinite(mag) || mag < 1e-12)
        errorexit("Cannot normalize: infinite, zero, or near-zero norm");
    else {
        result[0] = v[0] / mag;
        result[1] = v[1] / mag;
        result[2] = v[2] / mag;
    }
}


int safe_normalize(double v[3], double result[3])
{
    double mag = norm(v, 3);

    if (!isfinite(mag) || mag < 1e-12)
        return 0;

    for (int k = 0; k < 3; k++)
        result[k] = v[k] / mag;

    return 1;
}


void cross_product(double a[3], double b[3], double result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}


void vec_project_perp(double v[3], double n[3], double v_perp[3], int num_dim)
{
    double v_dot_n = dot_product(v, n, num_dim);
    for (int k = 0; k < num_dim; k++)
        v_perp[k] = v[k] - v_dot_n * n[k];
}


void create_rotation_matrix(double axis[3], double angle, double R[3][3])
{
    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0 - c;

    double axis_norm[3] = {0.0, 0.0, 0.0};
    normalize(axis, axis_norm);

    R[0][0] = t * axis_norm[0] * axis_norm[0] + c;
    R[0][1] = t * axis_norm[0] * axis_norm[1] - s * axis_norm[2];
    R[0][2] = t * axis_norm[0] * axis_norm[2] + s * axis_norm[1];

    R[1][0] = t * axis_norm[0] * axis_norm[1] + s * axis_norm[2];
    R[1][1] = t * axis_norm[1] * axis_norm[1] + c;
    R[1][2] = t * axis_norm[1] * axis_norm[2] - s * axis_norm[0];

    R[2][0] = t * axis_norm[0] * axis_norm[2] - s * axis_norm[1];
    R[2][1] = t * axis_norm[1] * axis_norm[2] + s * axis_norm[0];
    R[2][2] = t * axis_norm[2] * axis_norm[2] + c;
}


void rotate_vector(double v[3], double R[3][3], double result[3])
{
    double temp[3];
    temp[0] = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
    temp[1] = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
    temp[2] = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];

    result[0] = temp[0];
    result[1] = temp[1];
    result[2] = temp[2];
}


void align_vectors_rotation_matrix(double* v, double* v_target, double R[3][3])
{
    const double parallel_tol = 1e-12;

    double from[3];
    double to[3];

    if (!safe_normalize(v, from) || !safe_normalize(v_target, to))
        errorexit("Cannot align zero, non-finite, or near-zero vectors");

    const double cos_theta =
        clamp_unit(dot_product(from, to, 3));

    // Parallel vectors: no rotation.
    if (cos_theta >= 1.0 - parallel_tol) {
        R[0][0] = 1.0; R[0][1] = 0.0; R[0][2] = 0.0;
        R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
        R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
        return;
    }

    double axis[3];

    // Antiparallel vectors: choose a stable perpendicular axis.
    if (cos_theta <= -1.0 + parallel_tol) {
        double reference[3] = {0.0, 0.0, 0.0};
        int least_aligned = 0;

        if (fabs(from[1]) < fabs(from[least_aligned]))
            least_aligned = 1;
        if (fabs(from[2]) < fabs(from[least_aligned]))
            least_aligned = 2;

        reference[least_aligned] = 1.0;
        cross_product(from, reference, axis);

        create_rotation_matrix(axis, acos(-1.0), R);
        return;
    }

    // General case.
    cross_product(from, to, axis);

    const double sin_theta = norm(axis, 3);
    const double angle = atan2(sin_theta, cos_theta);

    create_rotation_matrix(axis, angle, R);
}


int delta(int i, int j)
{
    return (i == j) ? 1 : 0;
}


double sign_double(double x) {
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}


double clamp0(double x)
{
    return (x < 0.0) ? 0.0 : x;
}


double clamp_unit(double x)
{
    if (x > 1.0) return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}


double clamp_double(double x, double xmin, double xmax)
{
    if (x < xmin) return xmin;
    if (x > xmax) return xmax;
    return x;
}


int almost_equal(double a, double b, double rel_eps)
{
    double diff = fabs(a - b);
    double scale = fmax(fabs(a), fabs(b));
    if (scale < 1.0) scale = 1.0;
    return diff <= rel_eps * scale;
}


double wrap_to_2pi(double x)
{
    const double twopi = 2.0 * M_PI;

    x = fmod(x, twopi);

    if (x < 0.0)
        x += twopi;

    return x;
}


double angle_difference_abs(double a, double b)
{
    double d = wrap_to_2pi(a - b);

    if (d > M_PI)
        d = 2.0 * M_PI - d;

    return fabs(d);
}


void map_from_orbital_basis(const double v_old[3],
    const double e_hat[3], const double q_hat[3], const double h_hat[3], double v_new[3])
{
    /*
     * Old canonical basis:
     *     x_old = periapsis/reference direction
     *     y_old = positive true-anomaly direction
     *     z_old = angular momentum direction
     *
     * New basis:
     *     x_new = e_hat
     *     y_new = q_hat = h_hat x e_hat
     *     z_new = h_hat
     */
    v_new[0] = v_old[0] * e_hat[0] + v_old[1] * q_hat[0] + v_old[2] * h_hat[0];
    v_new[1] = v_old[0] * e_hat[1] + v_old[1] * q_hat[1] + v_old[2] * h_hat[1];
    v_new[2] = v_old[0] * e_hat[2] + v_old[1] * q_hat[2] + v_old[2] * h_hat[2];
}


void choose_perpendicular_reference_vector(double h_hat[3], double e_hat[3])
{
    double tmp[3];

    if (fabs(h_hat[0]) < 0.9) {
        tmp[0] = 1.0;
        tmp[1] = 0.0;
        tmp[2] = 0.0;
    }
    else {
        tmp[0] = 0.0;
        tmp[1] = 1.0;
        tmp[2] = 0.0;
    }

    const double dot = dot_product(tmp, h_hat, 3);

    e_hat[0] = tmp[0] - dot * h_hat[0];
    e_hat[1] = tmp[1] - dot * h_hat[1];
    e_hat[2] = tmp[2] - dot * h_hat[2];

    normalize(e_hat, e_hat);
}


void orientation_vectors_from_angles(double inc, double Omega, double omega,
    double h_hat[3], double e_hat[3])
{
    const double si = sin(inc);
    const double ci = cos(inc);

    const double sO = sin(Omega);
    const double cO = cos(Omega);

    const double so = sin(omega);
    const double co = cos(omega);

    h_hat[0] = si * sO;
    h_hat[1] = -si * cO;
    h_hat[2] = ci;

    e_hat[0] = cO * co - sO * so * ci;
    e_hat[1] = sO * co + cO * so * ci;
    e_hat[2] = so * si;

    normalize(h_hat, h_hat);
    normalize(e_hat, e_hat);
}


void angles_from_orientation_vectors(double h_input[3], double e_input[3],
    double *inc_out, double *Omega_out, double *omega_out)
{
    double h[3] = {h_input[0], h_input[1], h_input[2]};
    double e[3] = {e_input[0], e_input[1], e_input[2]};

    normalize(h, h);

    const double edoth = dot_product(e, h, 3);

    e[0] -= edoth * h[0];
    e[1] -= edoth * h[1];
    e[2] -= edoth * h[2];

    normalize(e, e);

    const double inc = acos(clamp_unit(h[2]));
    const double sin_i = sin(inc);

    double Omega;
    double omega;

    if (fabs(sin_i) > 1e-12) {
        Omega = atan2(h[0], -h[1]);

        double n_hat[3] = {-h[1] / sin_i, h[0] / sin_i, 0.0};

        double n_cross_e[3];
        cross_product(n_hat, e, n_cross_e);

        const double cos_omega = dot_product(n_hat, e, 3);
        const double sin_omega = dot_product(n_cross_e, h, 3);

        omega = atan2(sin_omega, cos_omega);
    }
    else {
        Omega = 0.0;

        if (h[2] >= 0.0) {
            // Prograde coplanar: h = +z, e = (cos omega, sin omega, 0)
            omega = atan2(e[1], e[0]);
        }
        else {
            // Retrograde coplanar: h = -z, e = (cos omega, -sin omega, 0)
            omega = atan2(-e[1], e[0]);
        }
    }

    *inc_out = inc;
    *Omega_out = wrap_to_2pi(Omega);
    *omega_out = wrap_to_2pi(omega);
}


void normalize_and_project_orientation_vectors(double h_hat[3], double e_hat[3],
    double eccentricity)
{
    if (!isfinite(h_hat[0]) || !isfinite(h_hat[1]) || !isfinite(h_hat[2]))
        errorexit("Invalid h_hat: components must be finite");

    if (!isfinite(e_hat[0]) || !isfinite(e_hat[1]) || !isfinite(e_hat[2]))
        errorexit("Invalid e_hat: components must be finite");

    normalize(h_hat, h_hat);

    // Project e_hat into the plane perpendicular to h_hat
    const double edoth = dot_product(e_hat, h_hat, 3);

    e_hat[0] -= edoth * h_hat[0];
    e_hat[1] -= edoth * h_hat[1];
    e_hat[2] -= edoth * h_hat[2];

    if (norm(e_hat, 3) < 1e-12) {
        if (eccentricity < 1e-12)
            choose_perpendicular_reference_vector(h_hat, e_hat);
        else
            errorexit("Invalid eccentric-binary orientation: e_hat is parallel to h_hat");
    }
    else
        normalize(e_hat, e_hat);
}


// ------------------------------------------------------------------------------------------------
// Random
// ------------------------------------------------------------------------------------------------

double rng_uniform(unsigned long long *state)
{
    *state = 6364136223846793005ULL * (*state) + 1442695040888963407ULL;
    return (double)(*state >> 11) * (1.0 / 9007199254740992.0);
}


void random_unit_vector(unsigned long long *rng, double n[3])
{
    const double z = 2.0 * rng_uniform(rng) - 1.0;
    const double phi = 2.0 * M_PI * rng_uniform(rng);
    const double s = sqrt(fmax(0.0, 1.0 - z*z));

    n[0] = s * cos(phi);
    n[1] = s * sin(phi);
    n[2] = z;
}


// ------------------------------------------------------------------------------------------------
// Print output
// ------------------------------------------------------------------------------------------------

void print_divider()
{
    printf("----------------------------------------------------------------------------------\n");
}


void print_state_vector(const double *w0, int num_bodies, int num_dim)
{
    for (int i = 0; i < num_bodies; ++i) {
        // Positions
        printf("pos%-2d = ", i + 1);
        for (int j = 0; j < num_dim; ++j) {
            printf("% .16e", w0[i * num_dim + j]);   // fixed width, sign included
            if (j + 1 < num_dim) printf("  ");
        }
        printf("\n");

        // Momenta
        printf("p%-4d = ", i + 1);
        for (int j = 0; j < num_dim; ++j) {
            printf("% .16e", w0[num_bodies * num_dim + i * num_dim + j]);
            if (j + 1 < num_dim) printf("  ");
        }
        printf("\n");
    }
}


void print_progress_bar(int percent)
{
    const int bar_width = 75;
    int filled = (percent * bar_width) / 100;

    printf("\r[");
    for (int i = 0; i < bar_width; i++) {
        if (i < filled)
            printf("=");
        else
            printf(" ");
    }
    printf("] %d%%", percent);
    fflush(stdout);
}


void progress_bar_break_line()
{
    printf("\r\033[2K\n");
    fflush(stdout);
}


noreturn void errorexit_function(const char *file, int line, const char *s)
{
    fprintf(stderr, "Error: %s (%s:%d)\n", s, file, line);
    fflush(stderr);
    exit(EXIT_FAILURE);
}


// ------------------------------------------------------------------------------------------------
// System calls
// ------------------------------------------------------------------------------------------------

int get_executable_dir(char out_dir[PATH_MAX])
{
    char exe_path[PATH_MAX];

#ifdef __APPLE__
    uint32_t size = (uint32_t)sizeof(exe_path);
    if (_NSGetExecutablePath(exe_path, &size) != 0) {
        // Buffer too small (shouldn't happen with PATH_MAX, but handle anyway)
        errno = ENAMETOOLONG;
        return 1;
    }
    // _NSGetExecutablePath may return a path with ../; normalize it
    char resolved[PATH_MAX];
    if (!realpath(exe_path, resolved)) return 1;
    snprintf(exe_path, sizeof(exe_path), "%s", resolved);
#else
    ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
    if (len < 0) return 1;
    exe_path[len] = '\0';
#endif

    // dirname() may modify its input, so use a temp buffer
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof(tmp), "%s", exe_path);

    char *dir = dirname(tmp);
    if (!dir) return 1;

    // Normalize directory to an absolute path
    char resolved_dir[PATH_MAX];
    if (!realpath(dir, resolved_dir)) return 1;

    snprintf(out_dir, PATH_MAX, "%s", resolved_dir);
    return 0;
}


char* make_filepath(const char* outdir, const char* filename)
{
    size_t len = strlen(outdir);
    size_t extra_slash = (outdir[len - 1] == '/' || outdir[len - 1] == '\\') ? 0 : 1;
    char* filepath = malloc(len + extra_slash + strlen(filename) + 1);
    if (!filepath)
        errorexit("Filepath could not be allocated");
    sprintf(filepath, "%s%s%s", outdir, (extra_slash ? "/" : ""), filename);
    return filepath;
}


void path_join(char out[PATH_MAX], const char *a, const char *b)
{
    int n = snprintf(out, PATH_MAX, "%s/%s", a, b);
    if (n < 0 || n >= PATH_MAX) errorexit("Path too long");
}


void mkdir_or_die(const char *path, mode_t mode)
{
    if (mkdir(path, mode) != 0 && errno != EEXIST) {
        fprintf(stderr, "Error: mkdir('%s') failed: %s\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }
}
