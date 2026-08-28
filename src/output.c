/**
 * @file output.c
 * @brief Functions for the output of physical quantities
 *
 * Functions for the output of physical quantities, such as the initialization of the output files
 * and writing output values at specified times.
 *
 * TODO: Add more possible output quantities (angular momentum, eccentricity, ...).
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "parameters.h"
#include "physical_quantities.h"
#include "eom.h"


enum output_quantity {
    OUTPUT_MASS     = 1 << 0,
    OUTPUT_POSITION = 1 << 1,
    OUTPUT_MOMENTUM = 1 << 2,
    OUTPUT_VELOCITY = 1 << 3,
    OUTPUT_SPIN     = 1 << 4,
    OUTPUT_ENERGY   = 1 << 5,
    OUTPUT_MERGER   = 1 << 6
};


static unsigned int get_output_quantities(void)
{
    const char *value = get_parameter_string("output");
    char *copy = strdup(value);
    if (!copy)
        errorexit("Could not allocate memory while parsing parameter \"output\"");

    unsigned int quantities = 0;
    char *token = strtok(copy, " \t");
    while (token) {
        if (strcmp(token, "mass") == 0)
            quantities |= OUTPUT_MASS;
        else if (strcmp(token, "position") == 0)
            quantities |= OUTPUT_POSITION;
        else if (strcmp(token, "momentum") == 0)
            quantities |= OUTPUT_MOMENTUM;
        else if (strcmp(token, "velocity") == 0)
            quantities |= OUTPUT_VELOCITY;
        else if (strcmp(token, "spin") == 0)
            quantities |= OUTPUT_SPIN;
        else if (strcmp(token, "energy") == 0)
            quantities |= OUTPUT_ENERGY;
        else if (strcmp(token, "merger") == 0)
            quantities |= OUTPUT_MERGER;
        else {
            char message[256];
            snprintf(message, sizeof(message),
                "Unknown output quantity \"%s\"; choose from mass position momentum velocity "
                "spin energy merger (separated by spaces)", token);
            free(copy);
            errorexit(message);
        }
        token = strtok(NULL, " \t");
    }

    free(copy);
    if (quantities == 0)
        errorexit("Parameter \"output\" must contain at least one output quantity");
    return quantities;
}


static FILE *open_output_file(const char *outdir, const char *filename)
{
    char *path = make_filepath(outdir, filename);
    FILE *file = fopen(path, "w");
    if (!file) {
        char message[256];
        snprintf(message, sizeof(message), "Output file \"%s\" could not be created", path);
        free(path);
        errorexit(message);
    }
    free(path);
    return file;
}


/**
 * @brief Initializes the output files.
 *
 * Initializes the output files, by creating the files, opening them, and writing the column names.
 *
 * @param[in]   file_mass      Pointer to the file containing the masses
 * @param[in]   file_pos       Pointer to the file containing the particle positions
 * @param[in]   file_mom       Pointer to the file containing the particle momenta
 * @param[in]   file_vel       Pointer to the file containing the particle velocities
 * @param[in]   file_spin      Pointer to the file containing the particle spins
 * @param[in]   file_energy    Pointer to the file containing the particle energies
 * @param[in]   file_merger    Pointer to the file containing the merger information
 * @param[in]   ode_params     Parameter struct containing general information about the system
 */
void output_init(FILE** file_mass, FILE** file_pos, FILE** file_mom, FILE** file_vel,
    FILE** file_spin, FILE** file_energy, FILE** file_merger, struct ode_params* ode_params)
{
    const unsigned int quantities = get_output_quantities();
    const char *outdir = get_parameter_string("outdir");

    *file_mass = NULL;
    *file_pos = NULL;
    *file_mom = NULL;
    *file_vel = NULL;
    *file_spin = NULL;
    *file_energy = NULL;
    *file_merger = NULL;

    if (quantities & OUTPUT_MASS)
        *file_mass = open_output_file(outdir, "output_mass.dat");
    if (quantities & OUTPUT_POSITION)
        *file_pos = open_output_file(outdir, "output_pos.dat");
    if (quantities & OUTPUT_MOMENTUM)
        *file_mom = open_output_file(outdir, "output_mom.dat");
    if (quantities & OUTPUT_VELOCITY)
        *file_vel = open_output_file(outdir, "output_vel.dat");
    if (quantities & OUTPUT_SPIN)
        *file_spin = open_output_file(outdir, "output_spin.dat");
    if (quantities & OUTPUT_ENERGY)
        *file_energy = open_output_file(outdir, "output_energy.dat");
    if (quantities & OUTPUT_MERGER)
        *file_merger = open_output_file(outdir, "output_merger.dat");

    // Write masses into the corresponding file
    if (*file_mass) {
        for (int i = 0; i < ode_params->num_bodies_initial; i++)
            fprintf(*file_mass, "m%d = %lf\n", i, ode_params->masses[i]);
    }

    // Write position column names into the corresponding file
    if (*file_pos) {
        fprintf(*file_pos, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++) {
            fprintf(*file_pos, "x%d\ty%d\t", i, i);
            if (ode_params->num_dim == 3) fprintf(*file_pos, "z%d\t", i);
        }
    }

    // Write momentum column names into the corresponding file
    if (*file_mom) {
        fprintf(*file_mom, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++) {
            fprintf(*file_mom, "px%d\tpy%d\t", i, i);
            if (ode_params->num_dim == 3) fprintf(*file_mom, "pz%d\t", i);
        }
    }

    // Write velocity column names into the corresponding file
    if (*file_vel) {
        fprintf(*file_vel, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++) {
            fprintf(*file_vel, "vx%d\tvy%d\t", i, i);
            if (ode_params->num_dim == 3) fprintf(*file_vel, "vz%d\t", i);
        }
    }

    // Write spin column names into the corresponding file
    if (*file_spin) {
        fprintf(*file_spin, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++)
            fprintf(*file_spin, "sx%d\tsy%d\tsz%d\t", i, i, i);
    }

    // Write energy column names into the corresponding file
    if (*file_energy)
        fprintf(*file_energy, "t\tH");

    // Write merger column names into the corresponding file
    if (*file_merger) {
        fprintf(*file_merger,
            "t "
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

        for (int n = 0; n < 3; n++)
            fprintf(*file_merger, "s_i_%d ", n);

        for (int n = 0; n < 3; n++)
            fprintf(*file_merger, "s_j_%d ", n);

        for (int n = 0; n < 3; n++)
            fprintf(*file_merger, "s_rem_%d ", n);

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(*file_merger, "v_kick_%d ", n);

        fprintf(*file_merger, "r_ij\n");
    }
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
 * @param[in]   file_vel       Pointer to the file containing the particle velocities
 * @param[in]   file_spin      Pointer to the file containing the particle spins
 * @param[in]   file_energy    Pointer to the file containing the particle energies
 * @param[in]   ode_params     Parameter struct containing general information about the system
 * @param[in]   w              Current state of the full system, w = [positions, momenta]
 * @param[in]   t              Current time
 */
void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_vel, FILE* file_spin,
    FILE* file_energy, struct ode_params* ode_params, double* w, double t)
{
    int array_half = ode_params->num_dim * ode_params->num_bodies_initial;
    int spin_offset = 2 * array_half;
    int num_spin_components = 3 * ode_params->num_bodies_initial;

    // Write positions
    if (file_pos) {
        fprintf(file_pos, "\n%.20e\t", t);
        for (int i = 0; i < array_half; i++) {
            if (component_is_active(i, ode_params))
                fprintf(file_pos, "%.20e\t", w[i]);
            else
                fprintf(file_pos, "nan\t");
        }
    }

    // Write momenta
    if (file_mom) {
        fprintf(file_mom, "\n%.20e\t", t);
        for (int i = 0; i < array_half; i++) {
            if (component_is_active(i, ode_params))
                fprintf(file_mom, "%.20e\t", w[array_half + i]);
            else
                fprintf(file_mom, "nan\t");
        }
    }

    // Write coordinate velocities
    if (file_vel) {
        double velocities[array_half];
        compute_coordinate_velocities(w, ode_params, velocities);

        fprintf(file_vel, "\n%.20e\t", t);
        for (int i = 0; i < array_half; i++) {
            if (component_is_active(i, ode_params))
                fprintf(file_vel, "%.20e\t", velocities[i]);
            else
                fprintf(file_vel, "nan\t");
        }
    }

    // Write spins
    if (file_spin) {
        fprintf(file_spin, "\n%.20e\t", t);
        for (int i = 0; i < num_spin_components; i++) {
            if (ode_params->active[i / 3])
                fprintf(file_spin, "%.20e\t", w[spin_offset + i]);
            else
                fprintf(file_spin, "nan\t");
        }
    }

    // Write energy
    if (file_energy)
        fprintf(file_energy, "\n%.20e\t%.20e\t", t,
            total_energy_conservative(w, ode_params));
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
    const double *s_i_old,
    const double *s_j_old,
    const double *s_rem,
    const double *v_kick_kms,
    double r_ij
) {
    if (!file_merger)
        return;

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

    for (int n = 0; n < 3; n++)
        fprintf(file_merger, "%.16e ", s_i_old[n]);

    for (int n = 0; n < 3; n++)
        fprintf(file_merger, "%.16e ", s_j_old[n]);

    for (int n = 0; n < 3; n++)
        fprintf(file_merger, "%.16e ", s_rem[n]);

    for (int n = 0; n < num_dim; n++)
        fprintf(file_merger, "%.16e ", v_kick_kms[n]);

    fprintf(file_merger, "%.16e\n", r_ij);
}
