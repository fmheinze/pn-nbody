#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "integration.h"
#include "pn_eom.h"
#include "orbital_elements.h"


void general_output_init(FILE** file, struct ode_params* params) {
    *file = fopen("output.dat", "w");
    if (*file == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    // Write masses into the first line of the file
    fprintf(*file, "Masses:\t"); 
    for (int i = 0; i < params->num_bodies; i++)
        fprintf(*file, "m%d = %lf\t", i, params->masses[i]);

    // Write column names into the second line
    fprintf(*file, "\nt\t");
    for (int i = 0; i < params->num_bodies; i++) {
        fprintf(*file, "x%d\ty%d\t", i, i);
        if (params->num_dim == 3) fprintf(*file, "z%d\t", i);
    }
    for (int i = 0; i < params->num_bodies; i++) {
        fprintf(*file, "px%d\tpy%d\t", i, i);
        if (params->num_dim == 3) fprintf(*file, "pz%d\t", i);
    }
}


void general_output_write(FILE* file, struct ode_params* params, double* w, double t) {
    fprintf(file, "\n%.20e\t", t);
    for(int i = 0; i < 2 * params->num_dim * params->num_bodies; i++)
        fprintf(file, "%.20e\t", w[i]);
}


void eccentricity_output_init(FILE** file) {
    *file = fopen("eccentricity.dat", "w");
    if (*file == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    // Write column names into the first line
    fprintf(*file, "t\teccentricity");
}


void eccentricity_output_write(FILE* file, struct ode_params* params, double* w, double t) {
    int bodies[] = {0, 1}; // Fix: Bodies hardcoded!
    fprintf(file, "\n%.20e\t%.20e", t, eccentricity(bodies, w, params));
}


void semi_major_axes_output_init(FILE** file) {
    *file = fopen("semi_major_axes.dat", "w");
    if (*file == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    // Write column names into the first line
    fprintf(*file, "t\ta\ta1_com\ta2_com");
}


void semi_major_axes_output_write(FILE* file, struct ode_params* params, double* w, double t) {
    int bodies[] = {0, 1}; // Fix: Bodies hardcoded!
    double a, a1_com, a2_com;
    semi_major_axes(bodies, w, params, &a, &a1_com, &a2_com);
    fprintf(file, "\n%.20e\t%.20e\t%.20e\t%.20e", t, a, a1_com, a2_com);
}