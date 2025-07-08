#ifndef UTILS_H
#define UTILS_H

void allocate_vector(double** ptr, int num_elements);
void allocate_array(double*** ptr, int num_vectors, int num_elements);
void allocate_3d_array(double**** ptr, int num_arrays, int num_vectors, int num_elements);
void allocate_4d_array(double***** ptr, int num_3d_arrays, int num_arrays, int num_vectors, int num_elements);
void free_vector(double* ptr);
void free_array(double** ptr, int num_vectors);
void free_3d_array(double*** ptr, int num_arrays, int num_vectors);
void free_4d_array(double**** ptr, int num_3d_arrays, int num_arrays, int num_vectors);

double dot_product(double *a, double *b, int dim);
double norm(double *v, int dim);
void normalize(double v[3], double result[3]);
void cross_product(double a[3], double b[3], double result[3]);
void create_rotation_matrix(double axis[3], double angle, double R[3][3]);
void rotate_vector(double v[3], double R[3][3], double result[3]);
void align_vectors_rotation_matrix(double* v, double* v_target, double R[3][3]);
int kronecker_delta(int i, int j);

#endif