#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "integration.h"
#include "pn_eom.h"
#include "orbital_elements.h"
#include "parameters.h"
#include "pn_eom_hamiltonians.h"

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
    fprintf(*file_energy, "t\tH_base\tH_2PN_base\tUTT4\tH_2PN_tot\tH_tot\tdUdx\tdH_2PN_base_dx\tdH_2PN_dx\tcos_theta\tdHdx");

    // Clean up
    free(path_masses);
    free(path_pos);
    free(path_mom);
    free(path_energy);
}


void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_energy, struct ode_params* params, double* w, double t) {
    int array_half = params->num_dim * params->num_bodies;

    // Write time
    fprintf(file_pos, "\n%.20e\t", t);
    fprintf(file_mom, "\n%.20e\t", t);
    fprintf(file_energy, "\n%.20e\t", t);

    // Write positions
    for(int i = 0; i < array_half; i++)
        fprintf(file_pos, "%.20e\t", w[i]);

    // Write momenta
    for(int i = 0; i < array_half; i++)
        fprintf(file_mom, "%.20e\t", w[array_half + i]);

    // Write energy
    // double H_base, H_2PN_base, UTT4, norm_dUdx, norm_dH_2PN_base_dx, norm_dH_2PN_dx, cos_theta, norm_dHdx;
    // double dUdx[array_half];
    // double dwdt[2 * array_half];
    // for (int i = 0; i < 2 * array_half; i++)
    //     dwdt[i] = 0.0;
    // for (int i = 0; i < array_half; i++)
    //     dUdx[i] = 0.0;

    // total_energy_conservative(w, params, &H_base, &H_2PN_base, &UTT4);
    // compute_dUTT4_dx(w, params, dUdx);
    // update_eom_hamiltonian_cs(w, dwdt, H2PN_nbody_base_complex, 1e-30, params);

    // norm_dUdx = norm(dUdx, array_half);
    // norm_dH_2PN_base_dx = 0.0;
    // norm_dH_2PN_dx = 0.0;
    // for (int i = 0; i < array_half; i++) {
    //     norm_dH_2PN_base_dx += dwdt[array_half + i] * dwdt[array_half + i];
    //     norm_dH_2PN_dx += (dwdt[array_half + i] - dUdx[i]) * (dwdt[array_half + i] - dUdx[i]);
    // }
    // norm_dH_2PN_base_dx = sqrt(norm_dH_2PN_base_dx);
    // norm_dH_2PN_dx = sqrt(norm_dH_2PN_dx);
    
    // cos_theta = 0.0;
    // for (int i = 0; i < array_half; i++)
    //     cos_theta += (-dUdx[i] * dwdt[array_half + i]) / (norm_dUdx * norm_dH_2PN_base_dx);

    // norm_dHdx = 0.0;
    // update_eom_hamiltonian_cs(w, dwdt, H0PN_complex, 1e-30, params);
    // update_eom_hamiltonian_cs(w, dwdt, H1PN_complex, 1e-30, params);
    // for (int i = 0; i < array_half; i++) {
    //     norm_dHdx += (dwdt[array_half + i] - dUdx[i]) * (dwdt[array_half + i] - dUdx[i]);
    // }
    // norm_dHdx = sqrt(norm_dHdx);

    // fprintf(file_energy, "%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t", 
    //     H_base, H_2PN_base, UTT4, H_2PN_base + UTT4, H_base + UTT4, norm_dUdx, norm_dH_2PN_base_dx, norm_dH_2PN_dx, cos_theta, norm_dHdx);
}