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
 * @param[in]   file_merger    Pointer to the file containing the merger information
 * @param[in]   ode_params     Parameter struct containing general information about the system
 */
void output_init(FILE** file_mass, FILE** file_pos, FILE** file_mom, FILE** file_spin, 
    FILE** file_energy, FILE** file_merger, struct ode_params* ode_params) 
{
    // Create and open files
    char* outdir = get_parameter_string("outdir");
    char* path_masses  = make_filepath(outdir, "output_mass.dat");
    char* path_pos     = make_filepath(outdir, "output_pos.dat");
    char* path_mom     = make_filepath(outdir, "output_mom.dat");
    char* path_spin    = make_filepath(outdir, "output_spin.dat");
    char* path_energy  = make_filepath(outdir, "output_energy.dat");
    char* path_merger  = make_filepath(outdir, "output_merger.dat");

    *file_mass   = fopen(path_masses, "w");
    *file_pos    = fopen(path_pos, "w");
    *file_mom    = fopen(path_mom, "w");
    *file_spin   = fopen(path_spin, "w");
    *file_energy = fopen(path_energy, "w");
    *file_merger = fopen(path_merger, "w");

    if (!*file_mass || !*file_pos || !*file_mom || !*file_spin || !*file_energy || !*file_merger) {
        free(path_masses); free(path_pos); free(path_mom); free(path_spin); 
        free(path_energy); free(path_merger);
        errorexit("One or more of the output files could not be created");
    }

    // Write masses into the corresponding file
    for (int i = 0; i < ode_params->num_bodies_initial; i++)
        fprintf(*file_mass, "m%d = %lf\n", i, ode_params->masses[i]);

    // Write position column names into the corresponding file
    fprintf(*file_pos, "t\t");
    for (int i = 0; i < ode_params->num_bodies_initial; i++) {
        fprintf(*file_pos, "x%d\ty%d\t", i, i);
        if (ode_params->num_dim == 3) fprintf(*file_pos, "z%d\t", i);
    }

    // Write momentum column names into the corresponding file
    fprintf(*file_mom, "t\t");
    for (int i = 0; i < ode_params->num_bodies_initial; i++) {
        fprintf(*file_mom, "px%d\tpy%d\t", i, i);
        if (ode_params->num_dim == 3) fprintf(*file_mom, "pz%d\t", i);
    }

    // Write spin column names into the corresponding file
    fprintf(*file_spin, "t\t");
    for (int i = 0; i < ode_params->num_bodies_initial; i++) {
        fprintf(*file_spin, "sx%d\tsy%d\tsz%d\t", i, i, i);
    }

    // Write energy column names into the corresponding file
    fprintf(*file_energy, "t\tH");

    // Write merger column names into the corresponding file
    fprintf(*file_merger,
        "# t "
        "slot_i slot_j slot_remnant "
        "id_i id_j id_remnant "
        "gen_i gen_j gen_remnant "
        "m_i m_j m_remnant "
    );

    for (int n = 0; n < ode_params->num_dim; n++)
        fprintf(*file_merger, "x_i_%d ", n);

    for (int n = 0; n < ode_params->num_dim; n++)
        fprintf(*file_merger, "x_j_%d ", n);

    for (int n = 0; n < ode_params->num_dim; n++)
        fprintf(*file_merger, "x_rem_%d ", n);

    for (int n = 0; n < ode_params->num_dim; n++)
        fprintf(*file_merger, "p_i_%d ", n);

    for (int n = 0; n < ode_params->num_dim; n++)
        fprintf(*file_merger, "p_j_%d ", n);

    for (int n = 0; n < ode_params->num_dim; n++)
        fprintf(*file_merger, "p_rem_%d ", n);

    fprintf(*file_merger, "r_ij\n");

    // Clean up
    free(path_masses);
    free(path_pos);
    free(path_mom);
    free(path_energy);
    free(path_merger);
}


static int component_is_active(int idx, struct ode_params *ode_params)
{
    int num_dim = ode_params->num_dim;
    int num_bodies = ode_params->num_bodies_initial;
    int array_half = num_dim * num_bodies;

    int body;

    if (idx < array_half)
        body = idx / num_dim;
    else
        body = (idx - array_half) / num_dim;

    return ode_params->active[body];
}


/**
 * @brief Writes output quantities to a file for a given timestep.
 * 
 * @param[in]   file_pos       Pointer to the file containing the particle positions
 * @param[in]   file_mom       Pointer to the file containing the particle momenta
 * @param[in]   file_spin      Pointer to the file containing the particle spins
 * @param[in]   file_energy    Pointer to the file containing the particle energies
 * @param[in]   ode_params     Parameter struct containing general information about the system
 * @param[in]   w              Current state of the full system, w = [positions, momenta]
 * @param[in]   t              Current time
 */
void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_spin, FILE* file_energy, 
    struct ode_params* ode_params, double* w, double t) 
{
    int array_half = ode_params->num_dim * ode_params->num_bodies_initial;
    int spin_offset = 2 * array_half;

    // Write time
    fprintf(file_pos, "\n%.20e\t", t);
    fprintf(file_mom, "\n%.20e\t", t);
    fprintf(file_spin, "\n%.20e\t", t);
    fprintf(file_energy, "\n%.20e\t", t);

    // Write positions
    for (int i = 0; i < array_half; i++) {
        if (component_is_active(i, ode_params))
            fprintf(file_pos, "%.20e\t", w[i]);
        else
            fprintf(file_pos, "nan\t");
    }

    // Write momenta
    for (int i = 0; i < array_half; i++) {
        if (component_is_active(i, ode_params))
            fprintf(file_mom, "%.20e\t", w[array_half + i]);
        else
            fprintf(file_mom, "nan\t");
    }

    // Write spins
    for (int i = 0; i < array_half; i++) {
        if (component_is_active(i, ode_params))
            fprintf(file_spin, "%.20e\t", w[spin_offset + i]);
        else
            fprintf(file_spin, "nan\t");
    }

    // Write energy
    fprintf(file_energy, "%.20e\t", total_energy_conservative(w, ode_params));
}


void write_merger_header(FILE *file)
{
    fprintf(file,
        "t "
        "slot_i slot_j slot_remnant "
        "id_i id_j id_remnant "
        "gen_i gen_j gen_remnant "
        "m_i m_j m_remnant "
        "x_i y_i z_i "
        "x_j y_j z_j "
        "x_rem y_rem z_rem "
        "px_i py_i pz_i "
        "px_j py_j pz_j "
        "px_rem py_rem pz_rem "
        "r_ij\n"
    );
}


void output_write_merger_event(
    FILE *file_merger,
    double t,
    const struct ode_params *params,
    int i,
    int j,
    int slot_remnant,
    long long id_i,
    long long id_j,
    long long id_remnant,
    int gen_i,
    int gen_j,
    int gen_remnant,
    double mi,
    double mj,
    double m_remnant,
    const double *x_i_old,
    const double *x_j_old,
    const double *x_rem,
    const double *p_i_old,
    const double *p_j_old,
    const double *p_rem,
    double r_ij
) {
    const int num_dim = params->num_dim;

    fprintf(file_merger,
        "%.16e "
        "%d %d %d "
        "%lld %lld %lld "
        "%d %d %d "
        "%.16e %.16e %.16e ",
        t,
        i, j, slot_remnant,
        id_i, id_j, id_remnant,
        gen_i, gen_j, gen_remnant,
        mi, mj, m_remnant
    );

    for (int n = 0; n < num_dim; n++)
        fprintf(file_merger, "%.16e ", x_i_old[n]);

    for (int n = 0; n < num_dim; n++)
        fprintf(file_merger, "%.16e ", x_j_old[n]);

    for (int n = 0; n < num_dim; n++)
        fprintf(file_merger, "%.16e ", x_rem[n]);

    for (int n = 0; n < num_dim; n++)
        fprintf(file_merger, "%.16e ", p_i_old[n]);

    for (int n = 0; n < num_dim; n++)
        fprintf(file_merger, "%.16e ", p_j_old[n]);

    for (int n = 0; n < num_dim; n++)
        fprintf(file_merger, "%.16e ", p_rem[n]);

    fprintf(file_merger, "%.16e\n", r_ij);
}