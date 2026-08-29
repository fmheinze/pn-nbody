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
    OUTPUT_MASS                    = 1 << 0,
    OUTPUT_POSITION                = 1 << 1,
    OUTPUT_MOMENTUM                = 1 << 2,
    OUTPUT_VELOCITY                = 1 << 3,
    OUTPUT_SPIN                    = 1 << 4,
    OUTPUT_ENERGY                  = 1 << 5,
    OUTPUT_MERGER                  = 1 << 6,

    OUTPUT_ANGULAR_MOMENTUM                = 1 << 7,
    OUTPUT_ANGULAR_MOMENTUM_VECTOR         = 1 << 8,
    OUTPUT_ORBITAL_ANGULAR_MOMENTUM        = 1 << 9,
    OUTPUT_ORBITAL_ANGULAR_MOMENTUM_VECTOR = 1 << 10,
    OUTPUT_SPIN_ANGULAR_MOMENTUM           = 1 << 11,
    OUTPUT_SPIN_ANGULAR_MOMENTUM_VECTOR    = 1 << 12,

    OUTPUT_ANGULAR_MOMENTUM_COM                = 1 << 13,
    OUTPUT_ANGULAR_MOMENTUM_VECTOR_COM         = 1 << 14,
    OUTPUT_ORBITAL_ANGULAR_MOMENTUM_COM        = 1 << 15,
    OUTPUT_ORBITAL_ANGULAR_MOMENTUM_VECTOR_COM = 1 << 16
};


typedef enum OrbitQuantity {
    ORBIT_SEPARATION,
    ORBIT_RADIAL_VELOCITY,
    ORBIT_FREQUENCY_INSTANTANEOUS,

    ORBIT_PERICENTER_RADIAL,
    ORBIT_APOCENTER_RADIAL,
    ORBIT_SEMILATUS_RECTUM,

    ORBIT_SEMIMAJOR_AXIS_NEWTONIAN,
    ORBIT_SEMIMAJOR_AXIS_RADIAL,
    ORBIT_SEMIMAJOR_AXIS_RADIAL_ANALYTIC,
    ORBIT_SEMIMAJOR_AXIS_ENERGY,
    ORBIT_SEMIMAJOR_AXIS_MEAN_MOTION,
    ORBIT_SEMIMAJOR_AXIS_AZIMUTHAL,

    ORBIT_ECCENTRICITY_NEWTONIAN,
    ORBIT_ECCENTRICITY_RADIAL,
    ORBIT_ECCENTRICITY_RADIAL_ANALYTIC,
    ORBIT_ECCENTRICITY_TIME,
    ORBIT_ECCENTRICITY_ANGULAR,
    ORBIT_ECCENTRICITY_FREQUENCY,

    ORBIT_PEAK_GW_FREQUENCY_WEN_NEWTONIAN,
    ORBIT_PEAK_GW_FREQUENCY_WEN_RADIAL,

    ORBIT_ENERGY_BINARY,
    ORBIT_ENERGY_RATE_BINARY,
    ORBIT_ANGULAR_MOMENTUM_BINARY,
    ORBIT_ANGULAR_MOMENTUM_RATE_BINARY,

    ORBIT_LAPLACE_RUNGE_LENZ_VECTOR,
    ORBIT_ANGULAR_MOMENTUM_VECTOR,
    ORBIT_ECCENTRICITY_VECTOR,

    ORBIT_FREQUENCY_RADIAL,
    ORBIT_FREQUENCY_AZIMUTHAL,
    ORBIT_PERIOD_RADIAL,
    ORBIT_PERIOD_AZIMUTHAL,
    ORBIT_PERICENTER_ADVANCE,
    ORBIT_PRECESSION_RATE
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
        case ORBIT_SEPARATION:
            return "separation";
        case ORBIT_RADIAL_VELOCITY:
            return "radial_velocity";
        case ORBIT_FREQUENCY_INSTANTANEOUS:
            return "orbital_frequency_instantaneous";

        case ORBIT_PERICENTER_RADIAL:
            return "pericenter_radial";
        case ORBIT_APOCENTER_RADIAL:
            return "apocenter_radial";
        case ORBIT_SEMILATUS_RECTUM:
            return "semilatus_rectum";

        case ORBIT_SEMIMAJOR_AXIS_NEWTONIAN:
            return "semimajor_axis_newtonian";
        case ORBIT_SEMIMAJOR_AXIS_RADIAL:
            return "semimajor_axis_radial";
        case ORBIT_SEMIMAJOR_AXIS_RADIAL_ANALYTIC:
            return "semimajor_axis_radial_analytic";
        case ORBIT_SEMIMAJOR_AXIS_ENERGY:
            return "semimajor_axis_energy";
        case ORBIT_SEMIMAJOR_AXIS_MEAN_MOTION:
            return "semimajor_axis_mean_motion";
        case ORBIT_SEMIMAJOR_AXIS_AZIMUTHAL:
            return "semimajor_axis_azimuthal";

        case ORBIT_ECCENTRICITY_NEWTONIAN:
            return "eccentricity_newtonian";
        case ORBIT_ECCENTRICITY_RADIAL:
            return "eccentricity_radial";
        case ORBIT_ECCENTRICITY_RADIAL_ANALYTIC:
            return "eccentricity_radial_analytic";
        case ORBIT_ECCENTRICITY_TIME:
            return "eccentricity_time";
        case ORBIT_ECCENTRICITY_ANGULAR:
            return "eccentricity_angular";
        case ORBIT_ECCENTRICITY_FREQUENCY:
            return "eccentricity_frequency";

        case ORBIT_PEAK_GW_FREQUENCY_WEN_NEWTONIAN:
            return "peak_gw_frequency_wen_newtonian";
        case ORBIT_PEAK_GW_FREQUENCY_WEN_RADIAL:
            return "peak_gw_frequency_wen_radial";

        case ORBIT_ENERGY_BINARY:
            return "energy_binary";
        case ORBIT_ENERGY_RATE_BINARY:
            return "binary_energy_rate";
        case ORBIT_ANGULAR_MOMENTUM_BINARY:
            return "angular_momentum_binary";
        case ORBIT_ANGULAR_MOMENTUM_RATE_BINARY:
            return "binary_angular_momentum_rate";

        case ORBIT_LAPLACE_RUNGE_LENZ_VECTOR:
            return "laplace_runge_lenz_vector";
        case ORBIT_ANGULAR_MOMENTUM_VECTOR:
            return "angular_momentum_vector";
        case ORBIT_ECCENTRICITY_VECTOR:
            return "eccentricity_vector";

        case ORBIT_FREQUENCY_RADIAL:
            return "omega_radial";
        case ORBIT_FREQUENCY_AZIMUTHAL:
            return "omega_azimuthal";
        case ORBIT_PERIOD_RADIAL:
            return "period_radial";
        case ORBIT_PERIOD_AZIMUTHAL:
            return "period_azimuthal";
        case ORBIT_PERICENTER_ADVANCE:
            return "pericenter_advance";
        case ORBIT_PRECESSION_RATE:
            return "precession_rate";
    }
    return "unknown";
}


static int parse_orbit_quantity(const char *token, OrbitQuantity *quantity)
{
    if (strcmp(token, "separation") == 0
        || strcmp(token, "binary_separation") == 0
        || strcmp(token, "r") == 0) {
        *quantity = ORBIT_SEPARATION;
        return 1;
    }
    if (strcmp(token, "radial_velocity") == 0
        || strcmp(token, "separation_rate") == 0
        || strcmp(token, "r_dot") == 0) {
        *quantity = ORBIT_RADIAL_VELOCITY;
        return 1;
    }
    if (strcmp(token, "orbital_frequency_instantaneous") == 0
        || strcmp(token, "instantaneous_orbital_frequency") == 0
        || strcmp(token, "omega_instantaneous") == 0
        || strcmp(token, "Omega_instantaneous") == 0) {
        *quantity = ORBIT_FREQUENCY_INSTANTANEOUS;
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
    if (strcmp(token, "semilatus_rectum") == 0
        || strcmp(token, "radial_semilatus_rectum") == 0
        || strcmp(token, "semilatus_rectum_radial") == 0) {
        *quantity = ORBIT_SEMILATUS_RECTUM;
        return 1;
    }

    if (strcmp(token, "semimajor_axis_newtonian") == 0
        || strcmp(token, "semi_major_axis_newtonian") == 0
        || strcmp(token, "newtonian_semimajor_axis") == 0
        || strcmp(token, "a_newtonian") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_NEWTONIAN;
        return 1;
    }
    if (strcmp(token, "semimajor_axis_radial") == 0
        || strcmp(token, "radial_semimajor_axis") == 0
        || strcmp(token, "semimajor_axis_radial_pn") == 0
        || strcmp(token, "a_r") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_RADIAL;
        return 1;
    }
    if (strcmp(token, "semimajor_axis_radial_analytic") == 0
        || strcmp(token, "radial_semimajor_axis_analytic") == 0
        || strcmp(token, "semimajor_axis_radial_pn_expanded") == 0
        || strcmp(token, "a_r_analytic") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_RADIAL_ANALYTIC;
        return 1;
    }
    if (strcmp(token, "semimajor_axis_energy") == 0
        || strcmp(token, "energy_semimajor_axis") == 0
        || strcmp(token, "a_E") == 0
        || strcmp(token, "a_energy") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_ENERGY;
        return 1;
    }
    if (strcmp(token, "semimajor_axis_mean_motion") == 0
        || strcmp(token, "mean_motion_semimajor_axis") == 0
        || strcmp(token, "semimajor_axis_n") == 0
        || strcmp(token, "a_n") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_MEAN_MOTION;
        return 1;
    }
    if (strcmp(token, "semimajor_axis_azimuthal") == 0
        || strcmp(token, "azimuthal_semimajor_axis") == 0
        || strcmp(token, "semimajor_axis_phi") == 0
        || strcmp(token, "a_phi") == 0) {
        *quantity = ORBIT_SEMIMAJOR_AXIS_AZIMUTHAL;
        return 1;
    }

    if (strcmp(token, "eccentricity_newtonian") == 0
        || strcmp(token, "newtonian_eccentricity") == 0
        || strcmp(token, "e_newtonian") == 0) {
        *quantity = ORBIT_ECCENTRICITY_NEWTONIAN;
        return 1;
    }
    if (strcmp(token, "eccentricity_radial") == 0
        || strcmp(token, "radial_eccentricity") == 0
        || strcmp(token, "eccentricity_radial_pn") == 0
        || strcmp(token, "e_r") == 0) {
        *quantity = ORBIT_ECCENTRICITY_RADIAL;
        return 1;
    }
    if (strcmp(token, "eccentricity_radial_analytic") == 0
        || strcmp(token, "radial_eccentricity_analytic") == 0
        || strcmp(token, "eccentricity_radial_pn_expanded") == 0
        || strcmp(token, "e_r_analytic") == 0) {
        *quantity = ORBIT_ECCENTRICITY_RADIAL_ANALYTIC;
        return 1;
    }
    if (strcmp(token, "eccentricity_time") == 0
        || strcmp(token, "time_eccentricity") == 0
        || strcmp(token, "eccentricity_time_pn") == 0
        || strcmp(token, "e_t") == 0) {
        *quantity = ORBIT_ECCENTRICITY_TIME;
        return 1;
    }
    if (strcmp(token, "eccentricity_angular") == 0
        || strcmp(token, "angular_eccentricity") == 0
        || strcmp(token, "eccentricity_angular_pn") == 0
        || strcmp(token, "eccentricity_phi") == 0
        || strcmp(token, "e_phi") == 0) {
        *quantity = ORBIT_ECCENTRICITY_ANGULAR;
        return 1;
    }
    if (strcmp(token, "eccentricity_frequency") == 0
        || strcmp(token, "frequency_eccentricity") == 0
        || strcmp(token, "eccentricity_omega") == 0
        || strcmp(token, "e_Omega") == 0
        || strcmp(token, "e_omega") == 0) {
        *quantity = ORBIT_ECCENTRICITY_FREQUENCY;
        return 1;
    }

    if (strcmp(token, "peak_gw_frequency_wen_newtonian") == 0
        || strcmp(token, "gw_peak_frequency_wen_newtonian") == 0
        || strcmp(token, "f_peak_wen_newtonian") == 0) {
        *quantity = ORBIT_PEAK_GW_FREQUENCY_WEN_NEWTONIAN;
        return 1;
    }
    if (strcmp(token, "peak_gw_frequency_wen_radial") == 0
        || strcmp(token, "gw_peak_frequency_wen_radial") == 0
        || strcmp(token, "f_peak_wen_radial") == 0) {
        *quantity = ORBIT_PEAK_GW_FREQUENCY_WEN_RADIAL;
        return 1;
    }

    if (strcmp(token, "energy_binary") == 0
        || strcmp(token, "binary_energy") == 0
        || strcmp(token, "isolated_pair_energy") == 0
        || strcmp(token, "energy") == 0) {
        *quantity = ORBIT_ENERGY_BINARY;
        return 1;
    }
    if (strcmp(token, "binary_energy_rate") == 0
        || strcmp(token, "energy_rate_binary") == 0
        || strcmp(token, "energy_rate") == 0
        || strcmp(token, "dE_dt") == 0) {
        *quantity = ORBIT_ENERGY_RATE_BINARY;
        return 1;
    }
    if (strcmp(token, "angular_momentum_binary") == 0
        || strcmp(token, "binary_angular_momentum") == 0
        || strcmp(token, "isolated_pair_angular_momentum") == 0
        || strcmp(token, "angular_momentum") == 0
        || strcmp(token, "J") == 0) {
        *quantity = ORBIT_ANGULAR_MOMENTUM_BINARY;
        return 1;
    }
    if (strcmp(token, "binary_angular_momentum_rate") == 0
        || strcmp(token, "angular_momentum_rate_binary") == 0
        || strcmp(token, "angular_momentum_rate") == 0
        || strcmp(token, "dJ_dt") == 0) {
        *quantity = ORBIT_ANGULAR_MOMENTUM_RATE_BINARY;
        return 1;
    }

    if (strcmp(token, "laplace_runge_lenz_vector") == 0
        || strcmp(token, "lrl_vector") == 0
        || strcmp(token, "LRL_vector") == 0
        || strcmp(token, "A_vector") == 0) {
        *quantity = ORBIT_LAPLACE_RUNGE_LENZ_VECTOR;
        return 1;
    }
    if (strcmp(token, "angular_momentum_vector") == 0
        || strcmp(token, "binary_angular_momentum_vector") == 0
        || strcmp(token, "L_vector") == 0) {
        *quantity = ORBIT_ANGULAR_MOMENTUM_VECTOR;
        return 1;
    }
    if (strcmp(token, "eccentricity_vector") == 0
        || strcmp(token, "e_vector") == 0) {
        *quantity = ORBIT_ECCENTRICITY_VECTOR;
        return 1;
    }

    if (strcmp(token, "omega_radial") == 0
        || strcmp(token, "radial_frequency") == 0
        || strcmp(token, "omega_r") == 0
        || strcmp(token, "Omega_r") == 0) {
        *quantity = ORBIT_FREQUENCY_RADIAL;
        return 1;
    }
    if (strcmp(token, "omega_azimuthal") == 0
        || strcmp(token, "azimuthal_frequency") == 0
        || strcmp(token, "omega_phi") == 0
        || strcmp(token, "Omega_phi") == 0) {
        *quantity = ORBIT_FREQUENCY_AZIMUTHAL;
        return 1;
    }
    if (strcmp(token, "period_radial") == 0
        || strcmp(token, "radial_period") == 0
        || strcmp(token, "P_r") == 0) {
        *quantity = ORBIT_PERIOD_RADIAL;
        return 1;
    }
    if (strcmp(token, "period_azimuthal") == 0
        || strcmp(token, "azimuthal_period") == 0
        || strcmp(token, "P_phi") == 0) {
        *quantity = ORBIT_PERIOD_AZIMUTHAL;
        return 1;
    }
    if (strcmp(token, "pericenter_advance") == 0
        || strcmp(token, "periastron_advance") == 0
        || strcmp(token, "delta_varpi") == 0) {
        *quantity = ORBIT_PERICENTER_ADVANCE;
        return 1;
    }
    if (strcmp(token, "precession_rate") == 0
        || strcmp(token, "pericenter_precession_rate") == 0
        || strcmp(token, "periastron_precession_rate") == 0
        || strcmp(token, "Omega_prec") == 0) {
        *quantity = ORBIT_PRECESSION_RATE;
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
            char message[768];
            snprintf(message, sizeof(message), "Unknown quantity \"%s\" for %s", token, name);
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
    for (size_t i = 0; i < orbit->num_quantities; ++i) {
        const OrbitQuantity quantity = orbit->quantities[i];
        const char *name = orbit_quantity_name(quantity);
        if (quantity == ORBIT_LAPLACE_RUNGE_LENZ_VECTOR
            || quantity == ORBIT_ANGULAR_MOMENTUM_VECTOR
            || quantity == ORBIT_ECCENTRICITY_VECTOR) {
            for (int axis = 0; axis < 3; ++axis)
                fprintf(orbit->file, "%s_%c\t", name, "xyz"[axis]);
        }
        else
            fprintf(orbit->file, "%s\t", name);
    }
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
        else if (strcmp(token, "angular_momentum") == 0)
            quantities |= OUTPUT_ANGULAR_MOMENTUM;
        else if (strcmp(token, "angular_momentum_vector") == 0)
            quantities |= OUTPUT_ANGULAR_MOMENTUM_VECTOR;
        else if (strcmp(token, "orbital_angular_momentum") == 0)
            quantities |= OUTPUT_ORBITAL_ANGULAR_MOMENTUM;
        else if (strcmp(token, "orbital_angular_momentum_vector") == 0)
            quantities |= OUTPUT_ORBITAL_ANGULAR_MOMENTUM_VECTOR;
        else if (strcmp(token, "spin_angular_momentum") == 0)
            quantities |= OUTPUT_SPIN_ANGULAR_MOMENTUM;
        else if (strcmp(token, "spin_angular_momentum_vector") == 0)
            quantities |= OUTPUT_SPIN_ANGULAR_MOMENTUM_VECTOR;
        else if (strcmp(token, "angular_momentum_com") == 0)
            quantities |= OUTPUT_ANGULAR_MOMENTUM_COM;
        else if (strcmp(token, "angular_momentum_vector_com") == 0)
            quantities |= OUTPUT_ANGULAR_MOMENTUM_VECTOR_COM;
        else if (strcmp(token, "orbital_angular_momentum_com") == 0)
            quantities |= OUTPUT_ORBITAL_ANGULAR_MOMENTUM_COM;
        else if (strcmp(token, "orbital_angular_momentum_vector_com") == 0)
            quantities |= OUTPUT_ORBITAL_ANGULAR_MOMENTUM_VECTOR_COM;
        else if (strcmp(token, "merger") == 0)
            quantities |= OUTPUT_MERGER;
        else if (strncmp(token, "orbit_", strlen("orbit_")) == 0)
            add_orbit_output(output, token, outdir, ode_params);
        else {
            char message[512];
            snprintf(message, sizeof(message), "Unknown output quantity \"%s\"", token);
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
    if (quantities & OUTPUT_ANGULAR_MOMENTUM)
        output->file_angular_momentum =
            open_output_file(outdir, "output_angular_momentum.dat");
    if (quantities & OUTPUT_ANGULAR_MOMENTUM_VECTOR)
        output->file_angular_momentum_vector =
            open_output_file(outdir, "output_angular_momentum_vector.dat");
    if (quantities & OUTPUT_ORBITAL_ANGULAR_MOMENTUM)
        output->file_orbital_angular_momentum =
            open_output_file(outdir, "output_orbital_angular_momentum.dat");
    if (quantities & OUTPUT_ORBITAL_ANGULAR_MOMENTUM_VECTOR)
        output->file_orbital_angular_momentum_vector =
            open_output_file(outdir, "output_orbital_angular_momentum_vector.dat");
    if (quantities & OUTPUT_SPIN_ANGULAR_MOMENTUM)
        output->file_spin_angular_momentum =
            open_output_file(outdir, "output_spin_angular_momentum.dat");
    if (quantities & OUTPUT_SPIN_ANGULAR_MOMENTUM_VECTOR)
        output->file_spin_angular_momentum_vector =
            open_output_file(outdir, "output_spin_angular_momentum_vector.dat");
    if (quantities & OUTPUT_ANGULAR_MOMENTUM_COM)
        output->file_angular_momentum_com =
            open_output_file(outdir, "output_angular_momentum_com.dat");
    if (quantities & OUTPUT_ANGULAR_MOMENTUM_VECTOR_COM)
        output->file_angular_momentum_vector_com =
            open_output_file(outdir, "output_angular_momentum_vector_com.dat");
    if (quantities & OUTPUT_ORBITAL_ANGULAR_MOMENTUM_COM)
        output->file_orbital_angular_momentum_com =
            open_output_file(outdir, "output_orbital_angular_momentum_com.dat");
    if (quantities & OUTPUT_ORBITAL_ANGULAR_MOMENTUM_VECTOR_COM)
        output->file_orbital_angular_momentum_vector_com =
            open_output_file(outdir, "output_orbital_angular_momentum_vector_com.dat");
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

    if (output->file_angular_momentum)
        fprintf(output->file_angular_momentum, "t\tangular_momentum");

    if (output->file_angular_momentum_vector)
        fprintf(output->file_angular_momentum_vector,
            "t\tangular_momentum_vector_x\tangular_momentum_vector_y\t"
            "angular_momentum_vector_z");

    if (output->file_orbital_angular_momentum)
        fprintf(output->file_orbital_angular_momentum, "t\torbital_angular_momentum");

    if (output->file_orbital_angular_momentum_vector)
        fprintf(output->file_orbital_angular_momentum_vector,
            "t\torbital_angular_momentum_vector_x\torbital_angular_momentum_vector_y\t"
            "orbital_angular_momentum_vector_z");

    if (output->file_spin_angular_momentum)
        fprintf(output->file_spin_angular_momentum, "t\tspin_angular_momentum");

    if (output->file_spin_angular_momentum_vector)
        fprintf(output->file_spin_angular_momentum_vector,
            "t\tspin_angular_momentum_vector_x\tspin_angular_momentum_vector_y\t"
            "spin_angular_momentum_vector_z");

    if (output->file_angular_momentum_com)
        fprintf(output->file_angular_momentum_com, "t\tangular_momentum_com");

    if (output->file_angular_momentum_vector_com)
        fprintf(output->file_angular_momentum_vector_com,
            "t\tangular_momentum_vector_com_x\tangular_momentum_vector_com_y\t"
            "angular_momentum_vector_com_z");

    if (output->file_orbital_angular_momentum_com)
        fprintf(output->file_orbital_angular_momentum_com,
            "t\torbital_angular_momentum_com");

    if (output->file_orbital_angular_momentum_vector_com)
        fprintf(output->file_orbital_angular_momentum_vector_com,
            "t\torbital_angular_momentum_vector_com_x\t"
            "orbital_angular_momentum_vector_com_y\torbital_angular_momentum_vector_com_z");

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
        || quantity == ORBIT_ECCENTRICITY_NEWTONIAN
        || quantity == ORBIT_PEAK_GW_FREQUENCY_WEN_NEWTONIAN
        || quantity == ORBIT_RADIAL_VELOCITY
        || quantity == ORBIT_FREQUENCY_INSTANTANEOUS;
}


static int orbit_quantity_needs_relative_state(OrbitQuantity quantity)
{
    return quantity == ORBIT_SEPARATION || orbit_quantity_needs_velocity(quantity);
}


static int orbit_quantity_needs_rhs(OrbitQuantity quantity)
{
    return quantity == ORBIT_ENERGY_RATE_BINARY
        || quantity == ORBIT_ANGULAR_MOMENTUM_RATE_BINARY;
}


static int orbit_quantity_needs_radial_turning_points(OrbitQuantity quantity)
{
    return quantity == ORBIT_PERICENTER_RADIAL
        || quantity == ORBIT_APOCENTER_RADIAL
        || quantity == ORBIT_SEMIMAJOR_AXIS_RADIAL
        || quantity == ORBIT_SEMILATUS_RECTUM
        || quantity == ORBIT_ECCENTRICITY_RADIAL
        || quantity == ORBIT_ECCENTRICITY_TIME
        || quantity == ORBIT_ECCENTRICITY_ANGULAR
        || quantity == ORBIT_ECCENTRICITY_FREQUENCY
        || quantity == ORBIT_PEAK_GW_FREQUENCY_WEN_RADIAL;
}


static int orbit_quantity_needs_frequency_state(OrbitQuantity quantity)
{
    return quantity == ORBIT_SEMIMAJOR_AXIS_MEAN_MOTION
        || quantity == ORBIT_SEMIMAJOR_AXIS_AZIMUTHAL
        || quantity == ORBIT_FREQUENCY_RADIAL
        || quantity == ORBIT_FREQUENCY_AZIMUTHAL
        || quantity == ORBIT_PERIOD_RADIAL
        || quantity == ORBIT_PERIOD_AZIMUTHAL
        || quantity == ORBIT_PERICENTER_ADVANCE
        || quantity == ORBIT_PRECESSION_RATE;
}


static int orbit_quantity_needs_vector_state(OrbitQuantity quantity)
{
    return quantity == ORBIT_LAPLACE_RUNGE_LENZ_VECTOR
        || quantity == ORBIT_ANGULAR_MOMENTUM_VECTOR
        || quantity == ORBIT_ECCENTRICITY_VECTOR;
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


static int any_active_orbit_needs_rhs(const OutputContext *output)
{
    for (size_t i = 0; i < output->num_orbits; ++i) {
        const OrbitOutput *orbit = &output->orbits[i];
        if (orbit->ended)
            continue;
        for (size_t j = 0; j < orbit->num_quantities; ++j) {
            if (orbit_quantity_needs_rhs(orbit->quantities[j]))
                return 1;
        }
    }
    return 0;
}


static void write_orbit_timestep(OrbitOutput *orbit, const struct ode_params *ode_params,
    const double *w, const double *velocities, const double *rhs, double t)
{
    int need_relative_state = 0;
    int need_pn_state = 0;
    int need_radial_turning_points = 0;
    int need_frequency_state = 0;
    int need_vector_state = 0;
    int need_rate_state = 0;
    for (size_t i = 0; i < orbit->num_quantities; ++i) {
        if (orbit_quantity_needs_relative_state(orbit->quantities[i]))
            need_relative_state = 1;
        else {
            need_pn_state = 1;
            if (orbit_quantity_needs_radial_turning_points(orbit->quantities[i]))
                need_radial_turning_points = 1;
            if (orbit_quantity_needs_frequency_state(orbit->quantities[i]))
                need_frequency_state = 1;
            if (orbit_quantity_needs_vector_state(orbit->quantities[i]))
                need_vector_state = 1;
            if (orbit_quantity_needs_rhs(orbit->quantities[i]))
                need_rate_state = 1;
        }
    }

    NewtonianRelativeOrbitState newtonian_state = {0};
    PNOrbitState pn_state = {0};
    PNOrbitFrequencyState frequency_state = {0};
    PNOrbitVectorState vector_state = {0};
    PNOrbitRateState rate_state = {0};

    if (!orbit->ended && need_relative_state)
        compute_newtonian_relative_orbit_state(w, velocities, ode_params, orbit->body_a,
            orbit->body_b, &newtonian_state);
    if (!orbit->ended && need_pn_state
        && compute_pn_orbit_state(w, ode_params, orbit->body_a, orbit->body_b, &pn_state)) {
        if (need_radial_turning_points)
            compute_pn_radial_turning_points(&pn_state);
        if (need_frequency_state)
            compute_pn_orbit_frequency_state(&pn_state, &frequency_state);
        if (need_vector_state)
            compute_pn_orbit_vector_state(&pn_state, &vector_state);
        if (need_rate_state)
            compute_pn_orbit_rate_state(rhs, ode_params, orbit->body_a, orbit->body_b,
                &pn_state, &rate_state);
    }

    fprintf(orbit->file, "\n%.20e\t", t);
    for (size_t i = 0; i < orbit->num_quantities; ++i) {
        double value = NAN;
        switch (orbit->quantities[i]) {
            case ORBIT_SEPARATION:
                value = separation_binary(&newtonian_state);
                break;
            case ORBIT_RADIAL_VELOCITY:
                value = radial_velocity_binary(&newtonian_state);
                break;
            case ORBIT_FREQUENCY_INSTANTANEOUS:
                value = orbital_frequency_instantaneous_binary(&newtonian_state);
                break;

            case ORBIT_PERICENTER_RADIAL:
                value = pericenter_radial_pn(&pn_state);
                break;
            case ORBIT_APOCENTER_RADIAL:
                value = apocenter_radial_pn(&pn_state);
                break;
            case ORBIT_SEMILATUS_RECTUM:
                value = semilatus_rectum_pn(&pn_state);
                break;

            case ORBIT_SEMIMAJOR_AXIS_NEWTONIAN:
                value = semimajor_axis_newtonian(&newtonian_state);
                break;
            case ORBIT_SEMIMAJOR_AXIS_RADIAL:
                value = semimajor_axis_radial_pn(&pn_state);
                break;
            case ORBIT_SEMIMAJOR_AXIS_RADIAL_ANALYTIC:
                value = semimajor_axis_radial_analytic_pn(&pn_state);
                break;
            case ORBIT_SEMIMAJOR_AXIS_ENERGY:
                value = semimajor_axis_energy_pn(&pn_state);
                break;
            case ORBIT_SEMIMAJOR_AXIS_MEAN_MOTION:
                value = semimajor_axis_mean_motion_pn(&frequency_state);
                break;
            case ORBIT_SEMIMAJOR_AXIS_AZIMUTHAL:
                value = semimajor_axis_azimuthal_pn(&frequency_state);
                break;

            case ORBIT_ECCENTRICITY_NEWTONIAN:
                value = eccentricity_newtonian(&newtonian_state);
                break;
            case ORBIT_ECCENTRICITY_RADIAL:
                value = eccentricity_radial_pn(&pn_state);
                break;
            case ORBIT_ECCENTRICITY_RADIAL_ANALYTIC:
                value = eccentricity_radial_analytic_pn(&pn_state);
                break;
            case ORBIT_ECCENTRICITY_TIME:
                value = eccentricity_time_pn(&pn_state);
                break;
            case ORBIT_ECCENTRICITY_ANGULAR:
                value = eccentricity_angular_pn(&pn_state);
                break;
            case ORBIT_ECCENTRICITY_FREQUENCY:
                value = eccentricity_frequency_pn(&pn_state);
                break;

            case ORBIT_PEAK_GW_FREQUENCY_WEN_NEWTONIAN:
                value = peak_gravitational_wave_frequency_wen_newtonian(&newtonian_state);
                break;
            case ORBIT_PEAK_GW_FREQUENCY_WEN_RADIAL:
                value = peak_gravitational_wave_frequency_wen_radial_pn(&pn_state);
                break;

            case ORBIT_ENERGY_BINARY:
                value = energy_binary_pn(&pn_state);
                break;
            case ORBIT_ENERGY_RATE_BINARY:
                value = energy_rate_binary_pn(&rate_state);
                break;
            case ORBIT_ANGULAR_MOMENTUM_BINARY:
                value = angular_momentum_binary_pn(&pn_state);
                break;
            case ORBIT_ANGULAR_MOMENTUM_RATE_BINARY:
                value = angular_momentum_rate_binary_pn(&rate_state);
                break;

            case ORBIT_LAPLACE_RUNGE_LENZ_VECTOR:
            case ORBIT_ANGULAR_MOMENTUM_VECTOR:
            case ORBIT_ECCENTRICITY_VECTOR: {
                double vector[3];
                if (orbit->quantities[i] == ORBIT_LAPLACE_RUNGE_LENZ_VECTOR)
                    laplace_runge_lenz_vector_binary_pn(&vector_state, vector);
                else if (orbit->quantities[i] == ORBIT_ANGULAR_MOMENTUM_VECTOR)
                    angular_momentum_vector_binary_pn(&vector_state, vector);
                else
                    eccentricity_vector_binary_pn(&vector_state, vector);
                for (int axis = 0; axis < 3; ++axis)
                    fprintf(orbit->file, "%.20e\t", vector[axis]);
                continue;
            }

            case ORBIT_FREQUENCY_RADIAL:
                value = frequency_radial_pn(&frequency_state);
                break;
            case ORBIT_FREQUENCY_AZIMUTHAL:
                value = frequency_azimuthal_pn(&frequency_state);
                break;
            case ORBIT_PERIOD_RADIAL:
                value = period_radial_pn(&frequency_state);
                break;
            case ORBIT_PERIOD_AZIMUTHAL:
                value = period_azimuthal_pn(&frequency_state);
                break;
            case ORBIT_PERICENTER_ADVANCE:
                value = pericenter_advance_pn(&frequency_state);
                break;
            case ORBIT_PRECESSION_RATE:
                value = precession_rate_pn(&frequency_state);
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
    const int need_rhs = any_active_orbit_needs_rhs(output);
    double velocities[array_half];
    double rhs[2*array_half];
    const double *velocity_values = NULL;
    const double *rhs_values = NULL;

    if (need_rhs) {
        rhs_pn_nbody(t, w, ode_params, rhs);
        rhs_values = rhs;
        if (need_velocities)
            velocity_values = rhs;
    }
    else if (need_velocities) {
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

    if (output->file_angular_momentum || output->file_angular_momentum_vector
        || output->file_orbital_angular_momentum
        || output->file_orbital_angular_momentum_vector
        || output->file_spin_angular_momentum
        || output->file_spin_angular_momentum_vector
        || output->file_angular_momentum_com
        || output->file_angular_momentum_vector_com
        || output->file_orbital_angular_momentum_com
        || output->file_orbital_angular_momentum_vector_com) {
        SystemAngularMomentumState angular_momentum_state = {0};
        compute_system_angular_momentum_state(w, ode_params, &angular_momentum_state);

        if (output->file_angular_momentum)
            fprintf(output->file_angular_momentum, "\n%.20e\t%.20e",
                t, total_angular_momentum(&angular_momentum_state));
        if (output->file_orbital_angular_momentum)
            fprintf(output->file_orbital_angular_momentum, "\n%.20e\t%.20e",
                t, orbital_angular_momentum(&angular_momentum_state));
        if (output->file_spin_angular_momentum)
            fprintf(output->file_spin_angular_momentum, "\n%.20e\t%.20e",
                t, spin_angular_momentum(&angular_momentum_state));
        if (output->file_angular_momentum_com)
            fprintf(output->file_angular_momentum_com, "\n%.20e\t%.20e",
                t, total_angular_momentum_com(&angular_momentum_state));
        if (output->file_orbital_angular_momentum_com)
            fprintf(output->file_orbital_angular_momentum_com, "\n%.20e\t%.20e",
                t, orbital_angular_momentum_com(&angular_momentum_state));

        double vector[3];
        if (output->file_angular_momentum_vector) {
            total_angular_momentum_vector(&angular_momentum_state, vector);
            fprintf(output->file_angular_momentum_vector,
                "\n%.20e\t%.20e\t%.20e\t%.20e", t,
                vector[0], vector[1], vector[2]);
        }
        if (output->file_orbital_angular_momentum_vector) {
            orbital_angular_momentum_vector(&angular_momentum_state, vector);
            fprintf(output->file_orbital_angular_momentum_vector,
                "\n%.20e\t%.20e\t%.20e\t%.20e", t,
                vector[0], vector[1], vector[2]);
        }
        if (output->file_spin_angular_momentum_vector) {
            spin_angular_momentum_vector(&angular_momentum_state, vector);
            fprintf(output->file_spin_angular_momentum_vector,
                "\n%.20e\t%.20e\t%.20e\t%.20e", t,
                vector[0], vector[1], vector[2]);
        }
        if (output->file_angular_momentum_vector_com) {
            total_angular_momentum_vector_com(&angular_momentum_state, vector);
            fprintf(output->file_angular_momentum_vector_com,
                "\n%.20e\t%.20e\t%.20e\t%.20e", t,
                vector[0], vector[1], vector[2]);
        }
        if (output->file_orbital_angular_momentum_vector_com) {
            orbital_angular_momentum_vector_com(&angular_momentum_state, vector);
            fprintf(output->file_orbital_angular_momentum_vector_com,
                "\n%.20e\t%.20e\t%.20e\t%.20e", t,
                vector[0], vector[1], vector[2]);
        }
    }

    // Newtonian definitions share one velocity evaluation. Relativistic definitions reuse the
    // pair invariants, and the root-based subset shares one turning-point solve.
    for (size_t i = 0; i < output->num_orbits; ++i)
        write_orbit_timestep(&output->orbits[i], ode_params, w, velocity_values, rhs_values, t);
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
    if (output->file_angular_momentum) fclose(output->file_angular_momentum);
    if (output->file_angular_momentum_vector) fclose(output->file_angular_momentum_vector);
    if (output->file_orbital_angular_momentum) fclose(output->file_orbital_angular_momentum);
    if (output->file_orbital_angular_momentum_vector)
        fclose(output->file_orbital_angular_momentum_vector);
    if (output->file_spin_angular_momentum) fclose(output->file_spin_angular_momentum);
    if (output->file_spin_angular_momentum_vector)
        fclose(output->file_spin_angular_momentum_vector);
    if (output->file_angular_momentum_com) fclose(output->file_angular_momentum_com);
    if (output->file_angular_momentum_vector_com)
        fclose(output->file_angular_momentum_vector_com);
    if (output->file_orbital_angular_momentum_com)
        fclose(output->file_orbital_angular_momentum_com);
    if (output->file_orbital_angular_momentum_vector_com)
        fclose(output->file_orbital_angular_momentum_vector_com);
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
