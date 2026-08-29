#ifndef PHYSICAL_QUANTITIES_H
#define PHYSICAL_QUANTITIES_H

struct ode_params;

// State preparation
typedef struct SystemAngularMomentumState {
    int valid;
    double orbital_vector[3];
    double orbital_com_vector[3];
    double spin_vector[3];
    double total_vector[3];
    double total_com_vector[3];
    double orbital_magnitude;
    double orbital_com_magnitude;
    double spin_magnitude;
    double total_magnitude;
    double total_com_magnitude;
} SystemAngularMomentumState;

typedef struct NewtonianRelativeOrbitState {
    int valid;
    int velocity_valid;
    int num_dim;
    double total_mass;
    double separation;
    double relative_speed_squared;
    double radial_velocity_product;
    double relative_position[3];
    double relative_velocity[3];
} NewtonianRelativeOrbitState;


typedef struct PNOrbitState {
    int valid;
    int radial_turning_points_valid;
    int use_1pn;
    int use_2pn;
    double total_mass;
    double symmetric_mass_ratio;
    double reduced_energy;
    double reduced_angular_momentum;
    double reduced_separation;
    double relative_position[3];
    double relative_momentum[3];
    double angular_momentum_vector[3];
    double pericenter;
    double apocenter;
} PNOrbitState;


typedef struct PNOrbitVectorState {
    int valid;
    double angular_momentum_vector[3];
    double laplace_runge_lenz_vector[3];
    double eccentricity_vector[3];
} PNOrbitVectorState;


typedef struct PNOrbitFrequencyState {
    int valid;
    double total_mass;
    double radial_frequency;
    double azimuthal_frequency;
    double periastron_advance_factor;
} PNOrbitFrequencyState;


typedef struct PNOrbitRateState {
    int valid;
    double energy_rate;
    double angular_momentum_rate;
} PNOrbitRateState;

int compute_system_angular_momentum_state(const double *w,
    const struct ode_params *ode_params, SystemAngularMomentumState *state);

int compute_newtonian_relative_orbit_state(const double *w, const double *velocities,
    const struct ode_params *ode_params, int body_a, int body_b,
    NewtonianRelativeOrbitState *state);

int compute_pn_orbit_state(const double *w, const struct ode_params *ode_params,
    int body_a, int body_b, PNOrbitState *state);

int compute_pn_radial_turning_points(PNOrbitState *state);

int compute_pn_orbit_frequency_state(const PNOrbitState *orbit_state,
    PNOrbitFrequencyState *frequency_state);

int compute_pn_orbit_vector_state(const PNOrbitState *orbit_state,
    PNOrbitVectorState *vector_state);

int compute_pn_orbit_rate_state(const double *rhs, const struct ode_params *ode_params,
    int body_a, int body_b, const PNOrbitState *orbit_state, PNOrbitRateState *rate_state);


// Total energy and angular momentum
double total_energy_conservative(double* w, struct ode_params* ode_params);

double orbital_angular_momentum(const SystemAngularMomentumState *state);
double orbital_angular_momentum_com(const SystemAngularMomentumState *state);
double spin_angular_momentum(const SystemAngularMomentumState *state);
double total_angular_momentum(const SystemAngularMomentumState *state);
double total_angular_momentum_com(const SystemAngularMomentumState *state);

int orbital_angular_momentum_vector(
    const SystemAngularMomentumState *state, double vector[3]);
int orbital_angular_momentum_vector_com(
    const SystemAngularMomentumState *state, double vector[3]);
int spin_angular_momentum_vector(
    const SystemAngularMomentumState *state, double vector[3]);
int total_angular_momentum_vector(
    const SystemAngularMomentumState *state, double vector[3]);
int total_angular_momentum_vector_com(
    const SystemAngularMomentumState *state, double vector[3]);

// Direct pair kinematics
double separation_binary(const NewtonianRelativeOrbitState *state);

double radial_velocity_binary(const NewtonianRelativeOrbitState *state);

double orbital_frequency_instantaneous_binary(const NewtonianRelativeOrbitState *state);


// Radial turning-point geometry
double pericenter_radial_pn(const PNOrbitState *state);

double apocenter_radial_pn(const PNOrbitState *state);

double semilatus_rectum_pn(const PNOrbitState *state);


// Semimajor-axis definitions
double semimajor_axis_newtonian(const NewtonianRelativeOrbitState *state);

double semimajor_axis_radial_pn(const PNOrbitState *state);

double semimajor_axis_radial_analytic_pn(const PNOrbitState *state);

double semimajor_axis_energy_pn(const PNOrbitState *state);

double semimajor_axis_mean_motion_pn(const PNOrbitFrequencyState *state);

double semimajor_axis_azimuthal_pn(const PNOrbitFrequencyState *state);


// Eccentricity definitions
double eccentricity_newtonian(const NewtonianRelativeOrbitState *state);

double eccentricity_radial_pn(const PNOrbitState *state);

double eccentricity_radial_analytic_pn(const PNOrbitState *state);

double eccentricity_time_pn(const PNOrbitState *state);

double eccentricity_angular_pn(const PNOrbitState *state);

double eccentricity_frequency_pn(const PNOrbitState *state);


// Characteristic gravitational-wave frequencies
double peak_gravitational_wave_frequency_wen_newtonian(
    const NewtonianRelativeOrbitState *state);

double peak_gravitational_wave_frequency_wen_radial_pn(const PNOrbitState *state);


// Pair energy
double energy_binary_pn(const PNOrbitState *state);

double energy_rate_binary_pn(const PNOrbitRateState *state);


// Pair angular momentum
double angular_momentum_binary_pn(const PNOrbitState *state);

double angular_momentum_rate_binary_pn(const PNOrbitRateState *state);

int angular_momentum_vector_binary_pn(const PNOrbitVectorState *state, double vector[3]);


// Other canonical pair vectors
int laplace_runge_lenz_vector_binary_pn(const PNOrbitVectorState *state, double vector[3]);

int eccentricity_vector_binary_pn(const PNOrbitVectorState *state, double vector[3]);


// Orbit-averaged frequencies, periods and precession
double frequency_radial_pn(const PNOrbitFrequencyState *state);

double frequency_azimuthal_pn(const PNOrbitFrequencyState *state);

double period_radial_pn(const PNOrbitFrequencyState *state);

double period_azimuthal_pn(const PNOrbitFrequencyState *state);

double periastron_advance_factor_pn(const PNOrbitFrequencyState *state);

double pericenter_advance_pn(const PNOrbitFrequencyState *state);

double precession_rate_pn(const PNOrbitFrequencyState *state);

#endif
