#ifndef UTILS_H
#define UTILS_H

#include "eom.h"
#include <complex.h>
#include <limits.h>

void allocate_vector(double** ptr, int num_elements);
void allocate_2d_array(double*** ptr, int num_vectors, int num_elements);
void allocate_3d_array(double**** ptr, int num_arrays, int num_vectors, int num_elements);
void allocate_4d_array(double***** ptr, int num_3d_arrays, int num_arrays, int num_vectors, int num_elements);
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
void cross_product(double a[3], double b[3], double result[3]);
void create_rotation_matrix(double axis[3], double angle, double R[3][3]);
void rotate_vector(double v[3], double R[3][3], double result[3]);
void align_vectors_rotation_matrix(double* v, double* v_target, double R[3][3]);
int delta(int i, int j);
double clamp0(double x);
int almost_equal(double a, double b, double rel_eps);

void print_divider();
void print_state_vector(const double *w0, int num_bodies, int num_dim);
void print_progress_bar(int percent);
void progress_bar_break_line();
int get_executable_dir(char out_dir[PATH_MAX]);
char* make_filepath(const char* outdir, const char* filename);
void mkdir_or_die(const char *path, mode_t mode);
void errorexit(const char *file, const int line, const char *s);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))

#endif