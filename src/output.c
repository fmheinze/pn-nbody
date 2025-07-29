#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "integration.h"
#include "pn_eom.h"
#include "orbital_elements.h"
#include "parameters.h"

/* Add feature to select which quantities to output */


char* make_filepath(const char* outdir, const char* filename) {
    size_t len = strlen(outdir);
    size_t extra_slash = (outdir[len - 1] == '/' || outdir[len - 1] == '\\') ? 0 : 1;
    char* filepath = malloc(len + extra_slash + strlen(filename) + 1);
    if (!filepath)
        errorexit("Filepath could not be allocated");
    sprintf(filepath, "%s%s%s", outdir, (extra_slash ? "/" : ""), filename);
    return filepath;
}


void output_init(FILE** file_masses, FILE** file_pos, FILE** file_mom, FILE** file_energy, struct ode_params* params) {
    char* outdir = get_parameter_string("outdir");
    char* path_masses  = make_filepath(outdir, "output_masses.dat");
    char* path_pos     = make_filepath(outdir, "output_pos.dat");
    char* path_mom     = make_filepath(outdir, "output_mom.dat");
    char* path_energy  = make_filepath(outdir, "output_energy.dat");

    *file_masses = fopen(path_masses, "w");
    *file_pos    = fopen(path_pos, "w");
    *file_mom    = fopen(path_mom, "w");
    *file_energy = fopen(path_energy, "w");

    if (!*file_masses || !*file_pos || !*file_mom || !*file_energy) {
        free(path_masses); free(path_pos); free(path_mom); free(path_energy);
        errorexit("One or more of the output files could not be created");
    }

    // Write masses into the corresponding file
    for (int i = 0; i < params->num_bodies; i++)
        fprintf(*file_masses, "m%d = %lf\n", i, params->masses[i]);

    // Write position column names into the corresponding file
    fprintf(*file_pos, "t\t");
    for (int i = 0; i < params->num_bodies; i++) {
        fprintf(*file_pos, "x%d\ty%d\t", i, i);
        if (params->num_dim == 3) fprintf(*file_pos, "z%d\t", i);
    }

    // Write momentum column names into the corresponding file
    fprintf(*file_mom, "t\t");
    for (int i = 0; i < params->num_bodies; i++) {
        fprintf(*file_mom, "px%d\tpy%d\t", i, i);
        if (params->num_dim == 3) fprintf(*file_mom, "pz%d\t", i);
    }

    // Write energy column names into the corresponding file
    fprintf(*file_energy, "t\tE");

    // Clean up
    free(path_masses);
    free(path_pos);
    free(path_mom);
    free(path_energy);
}


void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_energy, struct ode_params* params, double* w, double t) {
    // Write time
    fprintf(file_pos, "\n%.20e\t", t);
    fprintf(file_mom, "\n%.20e\t", t);
    fprintf(file_energy, "\n%.20e\t", t);

    // Write positions
    for(int i = 0; i < params->num_dim * params->num_bodies; i++)
        fprintf(file_pos, "%.20e\t", w[i]);

    // Write momenta
    for(int i = 0; i < params->num_dim * params->num_bodies; i++)
        fprintf(file_mom, "%.20e\t", w[params->num_dim * params->num_bodies + i]);

    // Write total energy
    fprintf(file_energy, "%.20e", total_energy_conservative(w, params));
}