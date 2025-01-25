#ifndef UTILS_H
#define UTILS_H

void allocate_vector(double** ptr, size_t num_elements);
void free_vector(double* ptr);

double dot_product(double *a, double *b, int dim);
void normalize(double v[3]);
void cross_product(double a[3], double b[3], double result[3]);
void create_rotation_matrix(double axis[3], double angle, double R[3][3]);
void rotate_vector(double v[3], double R[3][3], double result[3]);

#endif