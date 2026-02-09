/**
 * @file output.c
 * @brief Functions for the output of physical quantities
 *
 * Functions for the output of physical quantities, such as the initialization of the output files
 * and writing output values as specified times.
 * 
 * TODO: Allow the user to specify which quantities he wants to output. Also add more possible
 * output quantities (angular momentum, eccentricity, ...).
 */

#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "parameters.h"
#include "physical_quantities.h"
#include "eom.h"


/**
 * @brief Initializes the output files.
 * 
 * Initializes the output files, by creating the files, opening them, and writing the column names.
 * 
 * @param[in]   file_mass      Pointer to the file containing the masses
 * @param[in]   file_pos       Pointer to the file containing the particle positions
 * @param[in]   file_mom       Pointer to the file containing the particle momenta
 * @param[in]   file_energy    Pointer to the file containing the particle energies
 * @param[in]   ode_params     Parameter struct containing general information about the system
 */
void output_init(FILE** file_mass, FILE** file_pos, FILE** file_mom, FILE** file_energy, 
    struct ode_params* ode_params) 
{
    // Create and open files
    char* outdir = get_parameter_string("outdir");
    char* path_masses  = make_filepath(outdir, "output_mass.dat");
    char* path_pos     = make_filepath(outdir, "output_pos.dat");
    char* path_mom     = make_filepath(outdir, "output_mom.dat");
    char* path_energy  = make_filepath(outdir, "output_energy.dat");

    *file_mass   = fopen(path_masses, "w");
    *file_pos    = fopen(path_pos, "w");
    *file_mom    = fopen(path_mom, "w");
    *file_energy = fopen(path_energy, "w");

    if (!*file_mass || !*file_pos || !*file_mom || !*file_energy) {
        free(path_masses); free(path_pos); free(path_mom); free(path_energy);
        errorexit("One or more of the output files could not be created");
    }

    // Write masses into the corresponding file
    for (int i = 0; i < ode_params->num_bodies; i++)
        fprintf(*file_mass, "m%d = %lf\n", i, ode_params->masses[i]);

    // Write position column names into the corresponding file
    fprintf(*file_pos, "t\t");
    for (int i = 0; i < ode_params->num_bodies; i++) {
        fprintf(*file_pos, "x%d\ty%d\t", i, i);
        if (ode_params->num_dim == 3) fprintf(*file_pos, "z%d\t", i);
    }

    // Write momentum column names into the corresponding file
    fprintf(*file_mom, "t\t");
    for (int i = 0; i < ode_params->num_bodies; i++) {
        fprintf(*file_mom, "px%d\tpy%d\t", i, i);
        if (ode_params->num_dim == 3) fprintf(*file_mom, "pz%d\t", i);
    }

    // Write energy column names into the corresponding file
    fprintf(*file_energy, "t\tH");

    // Clean up
    free(path_masses);
    free(path_pos);
    free(path_mom);
    free(path_energy);
}

/**
 * @brief Writes output quantities to a file for a given timestep.
 * 
 * @param[in]   file_pos       Pointer to the file containing the particle positions
 * @param[in]   file_mom       Pointer to the file containing the particle momenta
 * @param[in]   file_energy    Pointer to the file containing the particle energies
 * @param[in]   ode_params     Parameter struct containing general information about the system
 * @param[in]   w              Current state of the full system, w = [positions, momenta]
 * @param[in]   t              Current time
 */
void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_energy, 
    struct ode_params* ode_params, double* w, double t) 
{
    int array_half = ode_params->num_dim * ode_params->num_bodies;

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
    fprintf(file_energy, "%.20e\t", total_energy_conservative(w, ode_params));
}
