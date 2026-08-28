#ifndef PHYSICAL_QUANTITIES_H
#define PHYSICAL_QUANTITIES_H

struct ode_params;

double total_energy_conservative(double* w, struct ode_params* ode_params);

typedef struct NewtonianRelativeOrbitState {
    int valid;
    int num_dim;
    double total_mass;
    double separation;
    double relative_speed_squared;
    double radial_velocity_product;
    double relative_position[3];
    double relative_velocity[3];
} NewtonianRelativeOrbitState;


typedef struct PNRadialOrbitState {
    int valid;
    double total_mass;
    double symmetric_mass_ratio;
    double reduced_energy;
    double reduced_angular_momentum;
    double pericenter;
    double apocenter;
} PNRadialOrbitState;


int compute_newtonian_relative_orbit_state(const double *w, const double *velocities,
    const struct ode_params *ode_params, int body_a, int body_b,
    NewtonianRelativeOrbitState *state);

double semimajor_axis_newtonian(const NewtonianRelativeOrbitState *state);

double eccentricity_newtonian(const NewtonianRelativeOrbitState *state);

int compute_pn_radial_orbit_state(const double *w, const struct ode_params *ode_params,
    int body_a, int body_b, PNRadialOrbitState *state);

double pericenter_radial_pn(const PNRadialOrbitState *state);

double apocenter_radial_pn(const PNRadialOrbitState *state);

double semimajor_axis_radial_pn(const PNRadialOrbitState *state);

double eccentricity_radial_pn(const PNRadialOrbitState *state);

#endif
