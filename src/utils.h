#ifndef UTILS_H
#define UTILS_H

struct ode_params;
#include <complex.h>
#include <limits.h>
#include <stdnoreturn.h>
#include <sys/types.h>

void allocate_vector(double** ptr, int num_elements);
void allocate_2d_array(double*** ptr, int num_vectors, int num_elements);
void allocate_3d_array(double**** ptr, int num_arrays, int num_vectors, int num_elements);
void allocate_4d_array(double***** ptr, int num_3d_arrays, int num_arrays, int num_vectors,
    int num_elements);
void free_vector(double* ptr);
void free_2d_array(double** ptr, int num_vectors);
void free_3d_array(double*** ptr, int num_arrays, int num_vectors);
void free_4d_array(double**** ptr, int num_3d_arrays, int num_arrays, int num_vectors);
void free_ode_params(struct ode_params* params);

double dot_product(double *a, double *b, int dim);
complex double dot_product_c(complex double *a, complex double *b, int dim);
double norm(double *v, int dim);
complex double norm_c(complex double *v, int dim);
void normalize(double v[3], double result[3]);
int safe_normalize(double v[3], double result[3]);
void cross_product(double a[3], double b[3], double result[3]);
void vec_project_perp(double v[3], double n[3], double v_perp[3], int num_dim);
void create_rotation_matrix(double axis[3], double angle, double R[3][3]);
void rotate_vector(double v[3], double R[3][3], double result[3]);
void align_vectors_rotation_matrix(double* v, double* v_target, double R[3][3]);
int delta(int i, int j);
double sign_double(double x);
double clamp0(double x);
double clamp_unit(double x);
double clamp_double(double x, double xmin, double xmax);
int almost_equal(double a, double b, double rel_eps);
double wrap_to_2pi(double x);
double angle_difference_abs(double a, double b);
void map_from_orbital_basis(const double v_old[3], 
    const double e_hat[3], const double q_hat[3], const double h_hat[3], double v_new[3]);
void orientation_vectors_from_angles(double inc, double Omega, double omega, 
    double h_hat[3], double e_hat[3]);   
void angles_from_orientation_vectors(double h_input[3], double e_input[3],
    double *inc_out, double *Omega_out, double *omega_out);
void normalize_and_project_orientation_vectors(double h_hat[3], double e_hat[3], 
    double eccentricity);

void print_divider(void);
void print_state_vector(const double *w0, int num_bodies, int num_dim);
void print_progress_bar(int percent);
void progress_bar_break_line(void);
noreturn void errorexit_function(const char *file, int line, const char *s);
#define errorexit(s) errorexit_function(__FILE__, __LINE__, (s))

int get_executable_dir(char out_dir[PATH_MAX]);
char* make_filepath(const char* outdir, const char* filename);
void path_join(char out[PATH_MAX], const char *a, const char *b);
void mkdir_or_die(const char *path, mode_t mode);

#endif
