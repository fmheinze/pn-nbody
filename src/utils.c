#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "utils.h"
#include "pn_eom.h"
#include "math.h"


// ------------------- Allocation/Free ------------------- //

void allocate_vector(double** ptr, int num_elements) {
    *ptr = (double *)malloc(num_elements * sizeof(double));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

void allocate_2d_array(double*** ptr, int num_vectors, int num_elements) {
    *ptr = (double **)malloc(num_vectors * sizeof(double *));
    if (*ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for 2D array\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < num_vectors; i++) {
        allocate_vector(&((*ptr)[i]), num_elements);
    }
}

void allocate_3d_array(double**** ptr, int num_arrays, int num_vectors, int num_elements) {
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

void allocate_4d_array(double***** ptr, int num_3d_arrays, int num_arrays, int num_vectors, int num_elements) {
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

void free_vector(double* ptr) {
    if (ptr != NULL)
        free(ptr);
}


void free_2d_array(double** ptr, int num_vectors) {
    if (ptr != NULL) {
        for (int i = 0; i < num_vectors; i++)
            free_vector(ptr[i]);
        free(ptr);
    }
}

void free_3d_array(double*** ptr, int num_arrays, int num_vectors) {
    if (ptr != NULL) {
        for (int i = 0; i < num_arrays; i++)
            free_2d_array(ptr[i], num_vectors);
        free(ptr);
    }
}

void free_4d_array(double**** ptr, int num_3d_arrays, int num_arrays, int num_vectors) {
    for (int i = 0; i < num_3d_arrays; i++)
        free_3d_array(ptr[i], num_arrays, num_vectors);
    free(ptr);
}

void free_ode_params(struct ode_params* params) {
    if (params->masses) {
        free_vector(params->masses);
        params->masses = NULL;
    }
    if (params->pn_terms) {
        free(params->pn_terms);
        params->pn_terms = NULL;
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


void normalize(double v[3], double result[3]) {
    double mag = norm(v, 3);
    if (mag > 0) {
        result[0] = v[0] / mag;
        result[1] = v[1] / mag;
        result[2] = v[2] / mag;
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

    double axis_norm[3];
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


void rotate_vector(double v[3], double R[3][3], double result[3]) {
    double temp[3];
    temp[0] = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
    temp[1] = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
    temp[2] = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];

    result[0] = temp[0];
    result[1] = temp[1];
    result[2] = temp[2];
}


void align_vectors_rotation_matrix(double* v, double* v_target, double R[3][3]) {
    double v_norm[3];
    double v_target_norm[3];
    normalize(v, v_norm);
    normalize(v_target, v_target_norm);

    double cos_theta = dot_product(v_norm, v_target_norm, 3);
    
    // If vectors are already aligned, return identity matrix
    if (fabs(cos_theta - 1.0) < 1e-6) {
        R[0][0] = R[1][1] = R[2][2] = 1;
        R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0;
        return;
    }

    // If vectors are opposite, find a perpendicular axis
    double axis[3];
    cross_product(v_norm, v_target_norm, axis);
    normalize(axis, axis);

    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    // Rodrigues' formula
    R[0][0] = cos_theta + axis[0] * axis[0] * (1 - cos_theta);
    R[0][1] = axis[0] * axis[1] * (1 - cos_theta) - axis[2] * sin_theta;
    R[0][2] = axis[0] * axis[2] * (1 - cos_theta) + axis[1] * sin_theta;

    R[1][0] = axis[1] * axis[0] * (1 - cos_theta) + axis[2] * sin_theta;
    R[1][1] = cos_theta + axis[1] * axis[1] * (1 - cos_theta);
    R[1][2] = axis[1] * axis[2] * (1 - cos_theta) - axis[0] * sin_theta;

    R[2][0] = axis[2] * axis[0] * (1 - cos_theta) - axis[1] * sin_theta;
    R[2][1] = axis[2] * axis[1] * (1 - cos_theta) + axis[0] * sin_theta;
    R[2][2] = cos_theta + axis[2] * axis[2] * (1 - cos_theta);
}


int kronecker_delta(int i, int j) {
    return (i == j) ? 1 : 0;
}


// --------------------- Miscellaneous --------------------- //

void print_divider() {
    printf("----------------------------------------------------------------------\n");
}


void print_progress_bar(int percent) {
    const int bar_width = 50;
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


void get_executable_path(char* buffer, size_t size) {
    ssize_t len = readlink("/proc/self/exe", buffer, size - 1);
    if (len != -1) {
        buffer[len] = '\0';
    } else {
        errorexit("Could not read executable path");
    }
}


/* Errorexit fuction*/
/* Note that utils.h defines a macro so that the user does not have
   to specify __FILE__ and __LINE__ for the location where the error occured */
#undef errorexit

// Switch between exit(1) (normal termination) and abort() (crash with core dump)
#define ABORT 0
#define EXIT 1
#define THE_END if(EXIT) exit(1);  else abort();

void errorexit(const char *file, const int line, const char *s)
{
    fprintf(stderr, "Error: %s  ", s);
    fprintf(stderr, "(%s, line %d)\n", file, line);
    fflush(stdout);
    fflush(stderr);
    THE_END
}

/* Do not write functions beyond this line because the undef/define 
   method for the errorexit function means that it should go last! */