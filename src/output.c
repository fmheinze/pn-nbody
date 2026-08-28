/**
 * @file output.c
 * @brief Functions for the output of physical quantities
 *
 * Functions for the output of physical quantities, such as the initialization of the output files
 * and writing output values at specified times.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "output.h"
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


typedef enum OrbitQuantity {
    ORBIT_SEMIMAJOR_AXIS_NEWTONIAN,
    ORBIT_ECCENTRICITY_NEWTONIAN,
    ORBIT_PERICENTER_RADIAL,
    ORBIT_APOCENTER_RADIAL,
    ORBIT_SEMIMAJOR_AXIS_RADIAL,
    ORBIT_ECCENTRICITY_RADIAL
} OrbitQuantity;


typedef struct OrbitOutput {
    char *name;
    int body_a;
    int body_b;
    int ended;
    FILE *file;
    OrbitQuantity *quantities;
    size_t num_quantities;
} OrbitOutput;


static const char *orbit_quantity_name(OrbitQuantity quantity)
{
    switch (quantity) {
        case ORBIT_SEMIMAJOR_AXIS_NEWTONIAN:
            return "semimajor_axis_newtonian";
        case ORBIT_ECCENTRICITY_NEWTONIAN:
            return "eccentricity_newtonian";
        case ORBIT_PERICENTER_RADIAL:
            return "pericenter_radial";
        case ORBIT_APOCENTER_RADIAL:
            return "apocenter_radial";
        case ORBIT_SEMIMAJOR_AXIS_RADIAL:
            return "semimajor_axis_radial";
        case ORBIT_ECCENTRICITY_RADIAL:
            return "eccentricity_radial";
    }
    return "unknown";
}


static int parse_orbit_quantity(const char *token, OrbitQuantity *quantity)
{
    if (strcmp(token, "semimajor_axis_newtonian") == 0
        || strcmp(token, "semi_major_axis_newtonian") == 0
        || strcmp(token, "newtonian_semimajor_axis") == 0
        || strcmp(token, "a_newtonian") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_NEWTONIAN;
        return 1;
    }
    if (strcmp(token, "eccentricity_newtonian") == 0
        || strcmp(token, "newtonian_eccentricity") == 0
        || strcmp(token, "e_newtonian") == 0) {
        *quantity = ORBIT_ECCENTRICITY_NEWTONIAN;
        return 1;
    }
    if (strcmp(token, "pericenter_radial") == 0
        || strcmp(token, "radial_pericenter") == 0
        || strcmp(token, "pericenter_radial_pn") == 0
        || strcmp(token, "r_p") == 0) {
        *quantity = ORBIT_PERICENTER_RADIAL;
        return 1;
    }
    if (strcmp(token, "apocenter_radial") == 0
        || strcmp(token, "radial_apocenter") == 0
        || strcmp(token, "apocenter_radial_pn") == 0
        || strcmp(token, "r_a") == 0) {
        *quantity = ORBIT_APOCENTER_RADIAL;
        return 1;
    }
    if (strcmp(token, "semimajor_axis_radial") == 0
        || strcmp(token, "radial_semimajor_axis") == 0
        || strcmp(token, "semimajor_axis_radial_pn") == 0
        || strcmp(token, "a_r") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_RADIAL;
        return 1;
    }
    if (strcmp(token, "eccentricity_radial") == 0
        || strcmp(token, "radial_eccentricity") == 0
        || strcmp(token, "eccentricity_radial_pn") == 0
        || strcmp(token, "e_r") == 0) {
        *quantity = ORBIT_ECCENTRICITY_RADIAL;
        return 1;
    }
    return 0;
}


static int parse_orbit_name(const char *name, int num_bodies, int *body_a, int *body_b)
{
    int user_a = 0;
    int user_b = 0;
    char trailing = '\0';
    if (sscanf(name, "orbit_%d_%d%c", &user_a, &user_b, &trailing) != 2)
        return 0;
    if (user_a < 1 || user_b < 1 || user_a > num_bodies || user_b > num_bodies
        || user_a == user_b)
        return 0;

    *body_a = user_a - 1;
    *body_b = user_b - 1;
    return 1;
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


static void add_orbit_output(OutputContext *output, const char *name, const char *outdir,
    const struct ode_params *ode_params)
{
    for (size_t i = 0; i < output->num_orbits; ++i) {
        if (strcmp(output->orbits[i].name, name) == 0)
            return;
    }

    int body_a = 0;
    int body_b = 0;
    if (!parse_orbit_name(name, ode_params->num_bodies_initial, &body_a, &body_b)) {
        char message[256];
        snprintf(message, sizeof(message),
            "Invalid orbit output \"%s\"; use orbit_X_Y with distinct 1-based body numbers", name);
        errorexit(message);
    }

    OrbitOutput *orbits = realloc(output->orbits,
        (output->num_orbits + 1)*sizeof(*output->orbits));
    if (orbits == NULL)
        errorexit("Could not allocate orbit output configuration");
    output->orbits = orbits;

    OrbitOutput *orbit = &output->orbits[output->num_orbits++];
    memset(orbit, 0, sizeof(*orbit));
    orbit->name = strdup(name);
    orbit->body_a = body_a;
    orbit->body_b = body_b;
    if (orbit->name == NULL)
        errorexit("Could not allocate orbit output name");

    const char *value = get_parameter_string(name);
    char *copy = strdup(value);
    if (copy == NULL)
        errorexit("Could not allocate orbit quantity list");

    char *saveptr = NULL;
    for (char *token = strtok_r(copy, " \t", &saveptr); token != NULL;
         token = strtok_r(NULL, " \t", &saveptr)) {
        OrbitQuantity quantity;
        if (!parse_orbit_quantity(token, &quantity)) {
            char message[384];
            snprintf(message, sizeof(message),
                "Unknown quantity \"%s\" for %s; choose from semimajor_axis_newtonian "
                "eccentricity_newtonian pericenter_radial apocenter_radial "
                "semimajor_axis_radial eccentricity_radial", token, name);
            free(copy);
            errorexit(message);
        }

        int duplicate = 0;
        for (size_t i = 0; i < orbit->num_quantities; ++i)
            duplicate |= orbit->quantities[i] == quantity;
        if (duplicate)
            continue;

        OrbitQuantity *quantities = realloc(orbit->quantities,
            (orbit->num_quantities + 1)*sizeof(*orbit->quantities));
        if (quantities == NULL) {
            free(copy);
            errorexit("Could not allocate orbit quantity list");
        }
        orbit->quantities = quantities;
        orbit->quantities[orbit->num_quantities++] = quantity;
    }
    free(copy);

    if (orbit->num_quantities == 0) {
        char message[256];
        snprintf(message, sizeof(message), "Parameter \"%s\" must list at least one quantity", name);
        errorexit(message);
    }

    const size_t filename_size = strlen("output_") + strlen(name) + strlen(".dat") + 1;
    char *filename = malloc(filename_size);
    if (filename == NULL)
        errorexit("Could not allocate orbit output filename");
    snprintf(filename, filename_size, "output_%s.dat", name);
    orbit->file = open_output_file(outdir, filename);
    free(filename);

    fprintf(orbit->file, "t\t");
    for (size_t i = 0; i < orbit->num_quantities; ++i)
        fprintf(orbit->file, "%s\t", orbit_quantity_name(orbit->quantities[i]));
}


static unsigned int parse_output_quantities(OutputContext *output, const char *outdir,
    const struct ode_params *ode_params)
{
    const char *value = get_parameter_string("output");
    char *copy = strdup(value);
    if (copy == NULL)
        errorexit("Could not allocate memory while parsing parameter \"output\"");

    unsigned int quantities = 0;
    int num_tokens = 0;
    char *saveptr = NULL;
    for (char *token = strtok_r(copy, " \t", &saveptr); token != NULL;
         token = strtok_r(NULL, " \t", &saveptr)) {
        ++num_tokens;
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
        else if (strncmp(token, "orbit_", strlen("orbit_")) == 0)
            add_orbit_output(output, token, outdir, ode_params);
        else {
            char message[320];
            snprintf(message, sizeof(message),
                "Unknown output quantity \"%s\"; choose from mass position momentum velocity "
                "spin energy merger or orbit_X_Y (separated by spaces)", token);
            free(copy);
            errorexit(message);
        }
    }

    free(copy);
    if (num_tokens == 0)
        errorexit("Parameter \"output\" must contain at least one output quantity");
    return quantities;
}


/**
 * @brief Initializes the output files.
 *
 * Initializes the output files, by creating the files, opening them, and writing the column names.
 *
 * @param[out]  output         Files and orbit-tracking state initialized by this function
 * @param[in]   ode_params     Parameter struct containing general information about the system
 */
void output_init(OutputContext *output, struct ode_params *ode_params)
{
    if (output == NULL || ode_params == NULL)
        errorexit("output_init received a NULL input");

    memset(output, 0, sizeof(*output));
    const char *outdir = get_parameter_string("outdir");
    const unsigned int quantities = parse_output_quantities(output, outdir, ode_params);

    if (quantities & OUTPUT_MASS)
        output->file_mass = open_output_file(outdir, "output_mass.dat");
    if (quantities & OUTPUT_POSITION)
        output->file_pos = open_output_file(outdir, "output_pos.dat");
    if (quantities & OUTPUT_MOMENTUM)
        output->file_mom = open_output_file(outdir, "output_mom.dat");
    if (quantities & OUTPUT_VELOCITY)
        output->file_vel = open_output_file(outdir, "output_vel.dat");
    if (quantities & OUTPUT_SPIN)
        output->file_spin = open_output_file(outdir, "output_spin.dat");
    if (quantities & OUTPUT_ENERGY)
        output->file_energy = open_output_file(outdir, "output_energy.dat");
    if (quantities & OUTPUT_MERGER)
        output->file_merger = open_output_file(outdir, "output_merger.dat");

    // Write masses into the corresponding file
    if (output->file_mass) {
        for (int i = 0; i < ode_params->num_bodies_initial; i++)
            fprintf(output->file_mass, "m%d = %lf\n", i, ode_params->masses[i]);
    }

    // Write position column names into the corresponding file
    if (output->file_pos) {
        fprintf(output->file_pos, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++) {
            fprintf(output->file_pos, "x%d\ty%d\t", i, i);
            if (ode_params->num_dim == 3) fprintf(output->file_pos, "z%d\t", i);
        }
    }

    // Write momentum column names into the corresponding file
    if (output->file_mom) {
        fprintf(output->file_mom, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++) {
            fprintf(output->file_mom, "px%d\tpy%d\t", i, i);
            if (ode_params->num_dim == 3) fprintf(output->file_mom, "pz%d\t", i);
        }
    }

    // Write velocity column names into the corresponding file
    if (output->file_vel) {
        fprintf(output->file_vel, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++) {
            fprintf(output->file_vel, "vx%d\tvy%d\t", i, i);
            if (ode_params->num_dim == 3) fprintf(output->file_vel, "vz%d\t", i);
        }
    }

    // Write spin column names into the corresponding file
    if (output->file_spin) {
        fprintf(output->file_spin, "t\t");
        for (int i = 0; i < ode_params->num_bodies_initial; i++)
            fprintf(output->file_spin, "sx%d\tsy%d\tsz%d\t", i, i, i);
    }

    // Write energy column names into the corresponding file
    if (output->file_energy)
        fprintf(output->file_energy, "t\tH");

    // Write merger column names into the corresponding file
    if (output->file_merger) {
        fprintf(output->file_merger,
            "t "
            "slot_i slot_j slot_remnant "
            "id_i id_j id_remnant "
            "gen_i gen_j gen_remnant "
            "m_i m_j m_remnant "
        );

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(output->file_merger, "x_i_%d ", n);

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(output->file_merger, "x_j_%d ", n);

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(output->file_merger, "x_rem_%d ", n);

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(output->file_merger, "p_i_%d ", n);

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(output->file_merger, "p_j_%d ", n);

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(output->file_merger, "p_rem_%d ", n);

        for (int n = 0; n < 3; n++)
            fprintf(output->file_merger, "s_i_%d ", n);

        for (int n = 0; n < 3; n++)
            fprintf(output->file_merger, "s_j_%d ", n);

        for (int n = 0; n < 3; n++)
            fprintf(output->file_merger, "s_rem_%d ", n);

        for (int n = 0; n < ode_params->num_dim; n++)
            fprintf(output->file_merger, "v_kick_%d ", n);

        fprintf(output->file_merger, "r_ij\n");
    }
}


static int component_is_active(int idx, const struct ode_params *ode_params)
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


static int orbit_quantity_needs_velocity(OrbitQuantity quantity)
{
    return quantity == ORBIT_SEMIMAJOR_AXIS_NEWTONIAN
        || quantity == ORBIT_ECCENTRICITY_NEWTONIAN;
}


static int any_active_orbit_needs_velocity(const OutputContext *output)
{
    for (size_t i = 0; i < output->num_orbits; ++i) {
        const OrbitOutput *orbit = &output->orbits[i];
        if (orbit->ended)
            continue;
        for (size_t j = 0; j < orbit->num_quantities; ++j) {
            if (orbit_quantity_needs_velocity(orbit->quantities[j]))
                return 1;
        }
    }
    return 0;
}


static void write_orbit_timestep(OrbitOutput *orbit, const struct ode_params *ode_params,
    const double *w, const double *velocities, double t)
{
    int need_newtonian_state = 0;
    int need_radial_state = 0;
    for (size_t i = 0; i < orbit->num_quantities; ++i) {
        if (orbit_quantity_needs_velocity(orbit->quantities[i]))
            need_newtonian_state = 1;
        else
            need_radial_state = 1;
    }

    NewtonianRelativeOrbitState newtonian_state = {0};
    PNRadialOrbitState radial_state = {0};

    if (!orbit->ended && need_newtonian_state)
        compute_newtonian_relative_orbit_state(w, velocities, ode_params, orbit->body_a,
            orbit->body_b, &newtonian_state);
    if (!orbit->ended && need_radial_state)
        compute_pn_radial_orbit_state(w, ode_params, orbit->body_a, orbit->body_b,
            &radial_state);

    fprintf(orbit->file, "\n%.20e\t", t);
    for (size_t i = 0; i < orbit->num_quantities; ++i) {
        double value = NAN;
        switch (orbit->quantities[i]) {
            case ORBIT_SEMIMAJOR_AXIS_NEWTONIAN:
                value = semimajor_axis_newtonian(&newtonian_state);
                break;
            case ORBIT_ECCENTRICITY_NEWTONIAN:
                value = eccentricity_newtonian(&newtonian_state);
                break;
            case ORBIT_PERICENTER_RADIAL:
                value = pericenter_radial_pn(&radial_state);
                break;
            case ORBIT_APOCENTER_RADIAL:
                value = apocenter_radial_pn(&radial_state);
                break;
            case ORBIT_SEMIMAJOR_AXIS_RADIAL:
                value = semimajor_axis_radial_pn(&radial_state);
                break;
            case ORBIT_ECCENTRICITY_RADIAL:
                value = eccentricity_radial_pn(&radial_state);
                break;
        }
        fprintf(orbit->file, "%.20e\t", value);
    }
}


/**
 * @brief Writes output quantities to a file for a given timestep.
 *
 * @param[in,out] output       Open output files and orbit-tracking state
 * @param[in]   ode_params     Parameter struct containing general information about the system
 * @param[in]   w              Current state of the full system, w = [positions, momenta]
 * @param[in]   t              Current time
 */
void output_write_timestep(OutputContext *output, struct ode_params *ode_params, double *w,
    double t)
{
    if (output == NULL || ode_params == NULL || w == NULL)
        errorexit("output_write_timestep received a NULL input");

    int array_half = ode_params->num_dim * ode_params->num_bodies_initial;
    int spin_offset = 2 * array_half;
    int num_spin_components = 3 * ode_params->num_bodies_initial;
    const int need_velocities =
        output->file_vel != NULL || any_active_orbit_needs_velocity(output);
    double velocities[array_half];
    const double *velocity_values = NULL;

    if (need_velocities) {
        compute_coordinate_velocities(w, ode_params, velocities);
        velocity_values = velocities;
    }

    // Write positions
    if (output->file_pos) {
        fprintf(output->file_pos, "\n%.20e\t", t);
        for (int i = 0; i < array_half; i++) {
            if (component_is_active(i, ode_params))
                fprintf(output->file_pos, "%.20e\t", w[i]);
            else
                fprintf(output->file_pos, "nan\t");
        }
    }

    // Write momenta
    if (output->file_mom) {
        fprintf(output->file_mom, "\n%.20e\t", t);
        for (int i = 0; i < array_half; i++) {
            if (component_is_active(i, ode_params))
                fprintf(output->file_mom, "%.20e\t", w[array_half + i]);
            else
                fprintf(output->file_mom, "nan\t");
        }
    }

    // Write coordinate velocities
    if (output->file_vel) {
        fprintf(output->file_vel, "\n%.20e\t", t);
        for (int i = 0; i < array_half; i++) {
            if (component_is_active(i, ode_params))
                fprintf(output->file_vel, "%.20e\t", velocity_values[i]);
            else
                fprintf(output->file_vel, "nan\t");
        }
    }

    // Write spins
    if (output->file_spin) {
        fprintf(output->file_spin, "\n%.20e\t", t);
        for (int i = 0; i < num_spin_components; i++) {
            if (ode_params->active[i / 3])
                fprintf(output->file_spin, "%.20e\t", w[spin_offset + i]);
            else
                fprintf(output->file_spin, "nan\t");
        }
    }

    // Write energy
    if (output->file_energy)
        fprintf(output->file_energy, "\n%.20e\t%.20e\t", t,
            total_energy_conservative(w, ode_params));

    // Write selected pairwise orbital elements. Newtonian definitions share one velocity evaluation;
    // radial PN definitions use canonical momenta directly and do not trigger that calculation.
    for (size_t i = 0; i < output->num_orbits; ++i)
        write_orbit_timestep(&output->orbits[i], ode_params, w, velocity_values, t);
}


/**
 * @brief Follow merger remnants for every configured orbit output.
 *
 * Each endpoint denotes the lineage of one initially selected body. If an endpoint participates in
 * a merger it is redirected to the remnant slot. If both endpoints enter the same remnant, the
 * binary no longer exists and all subsequent values in that orbit file are NAN.
 */
void output_follow_merger(OutputContext *output, int parent_a, int parent_b, int remnant)
{
    if (output == NULL)
        return;

    for (size_t i = 0; i < output->num_orbits; ++i) {
        OrbitOutput *orbit = &output->orbits[i];
        if (orbit->ended)
            continue;

        if (orbit->body_a == parent_a || orbit->body_a == parent_b)
            orbit->body_a = remnant;
        if (orbit->body_b == parent_a || orbit->body_b == parent_b)
            orbit->body_b = remnant;
        if (orbit->body_a == orbit->body_b)
            orbit->ended = 1;
    }
}


void output_close(OutputContext *output)
{
    if (output == NULL)
        return;

    if (output->file_mass) fclose(output->file_mass);
    if (output->file_pos) fclose(output->file_pos);
    if (output->file_mom) fclose(output->file_mom);
    if (output->file_vel) fclose(output->file_vel);
    if (output->file_spin) fclose(output->file_spin);
    if (output->file_energy) fclose(output->file_energy);
    if (output->file_merger) fclose(output->file_merger);

    for (size_t i = 0; i < output->num_orbits; ++i) {
        if (output->orbits[i].file) fclose(output->orbits[i].file);
        free(output->orbits[i].name);
        free(output->orbits[i].quantities);
    }
    free(output->orbits);
    memset(output, 0, sizeof(*output));
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
