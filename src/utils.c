#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "math.h"


// ------------------- Allocation/Free ------------------- //

// Allocates memory for a single vector
void allocate_vector(double** ptr, int num_elements) {
    *ptr = malloc(num_elements * sizeof(double));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

// Allocates memory for an array of vectors (a 2D array)
void allocate_array(double*** ptr, int num_vectors, int num_elements) {
    *ptr = (double **)malloc(num_vectors * sizeof(double *));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for 2D array\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < num_vectors; i++) {
        allocate_vector(&((*ptr)[i]), num_elements);
    }
}

// Allocates memory for a 3D array
void allocate_3d_array(double**** ptr, int num_arrays, int num_vectors, int num_elements) {
    // Allocate memory for the array of 2D arrays
    *ptr = (double ***)malloc(num_arrays * sizeof(double **));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for 3D array\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each 2D array
    for (int i = 0; i < num_arrays; i++) {
        allocate_array(&((*ptr)[i]), num_vectors, num_elements);
    }
}

// Frees memory for a single vector
void free_vector(double* ptr) {
    if (ptr != NULL) {
        free(ptr);
    }
}

// Frees memory for an array of vectors (a 2D array)
void free_array(double** ptr, int num_vectors) {
    if (ptr != NULL) {
        for (int i = 0; i < num_vectors; i++) {
            free_vector(ptr[i]);
        }
        free(ptr);
    }
}

// Frees memory for a 3D array
void free_3d_array(double*** ptr, int num_arrays, int num_vectors) {
    if (ptr != NULL) {
        for (int i = 0; i < num_arrays; i++) {
            free_array(ptr[i], num_vectors);
        }
        free(ptr);
    }
}



// --------------------- Math Utils --------------------- //

double dot_product(double *a, double *b, int dim) {
    double result = 0;
    for (int i = 0; i < dim; i++)
        result += a[i] * b[i];
    return result;
}


double norm(double *v, int dim) {
    return sqrt(dot_product(v, v, dim));
}


void normalize(double v[3]) {
    double mag = norm(v, 3);
    if (mag > 0) {
        v[0] /= mag;
        v[1] /= mag;
        v[2] /= mag;
    }
}


void cross_product(double a[3], double b[3], double result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}


void create_rotation_matrix(double axis[3], double angle, double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0 - c;

    normalize(axis);

    R[0][0] = t * axis[0] * axis[0] + c;
    R[0][1] = t * axis[0] * axis[1] - s * axis[2];
    R[0][2] = t * axis[0] * axis[2] + s * axis[1];

    R[1][0] = t * axis[0] * axis[1] + s * axis[2];
    R[1][1] = t * axis[1] * axis[1] + c;
    R[1][2] = t * axis[1] * axis[2] - s * axis[0];

    R[2][0] = t * axis[0] * axis[2] - s * axis[1];
    R[2][1] = t * axis[1] * axis[2] + s * axis[0];
    R[2][2] = t * axis[2] * axis[2] + c;
}


void rotate_vector(double v[3], double R[3][3], double result[3]) {
    result[0] = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
    result[1] = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
    result[2] = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];
}