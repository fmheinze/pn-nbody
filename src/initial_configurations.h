#ifndef INITIAL_CONFIGURATIONS_H
#define INITIAL_CONFIGURATIONS_H

struct ode_params;


// Binary struct
struct binary_params {
    double a;           // Semi-major axis
    double b;           // Semi-minor axis
    double e;           // Eccentricity
    double r_a;         // Apoapsis distance
    double r_p;         // Periapsis distance
    double p;           // Semi-latus rectum
    double f0;          // True anomaly
    double i;           // Inclination
    double Omega;       // Longitude of the ascending node
    double omega;       // Argument of periapsis

    double h_hat[3];    // Orbital angular momentum direction
    double e_hat[3];    // Periapsis / eccentricity-vector direction
};


// Struct for relativistic monoenergetic cluster
struct rel_mono_model {
    int n;
    int cap;

    double *x;       // dimensionless areal radius x = r sqrt(4 pi rho_c)
    double *mu;      // dimensionless mass mu = m sqrt(4 pi rho_c)
    double *z;       // z = E0 exp(-nu) / m0 = local Lorentz factor
    double *xiso;    // dimensionless isotropic radius
    double *cdf;     // proper rest-mass CDF

    double zc;
    double x_surf;
    double mu_surf;
    double mu0_surf;
    double compactness;   // R / M = x_surf / mu_surf
};


void position_binary(double com_pos[3], double h_input[3], double e_input[3], double w0[12]);

void ic_binary(struct ode_params* ode_params, struct binary_params* binary_params,
    double m1, double m2, double* w0);

void ic_hierarchical_triple(struct ode_params* ode_params,
    struct binary_params* inner_binary_params, struct binary_params* outer_binary_params,
    double* w0);

void ic_binary_single_scattering(struct ode_params* ode_params,
    struct binary_params* binary_params, double d0, double p0_rel, double b, double* w0);

void ic_binary_single_scattering_circ(double d0, double p0_rel, double b, double binary_phi0,
    double binary_r0, double binary_pt0, double binary_pr0, double* orientation, double* w0);

void ic_binary_binary_scattering(struct ode_params* ode_params,
    struct binary_params* binary1_params, struct binary_params* binary2_params, double d0,
    double p0_rel, double b, double* w0);

void ic_binary_binary_scattering_circ(double d0, double p0_rel, double b, double binary1_phi0,
    double binary1_r0, double binary1_pt0, double binary1_pr0, double* orientation_1,
    double binary_phi0_2, double binary2_r0, double binary2_pt0, double binary2_pr0,
    double* orientation_2, double* w0);

void ic_figure_eight_orbit(struct ode_params* params, double width, double* w0);

void ic_newtonian_plummer_cluster(struct ode_params* params, double compactness,
    double virial_ratio, double rmax_factor, double min_sep_factor,
    unsigned long long seed, double* w0);

void ic_relativistic_monoenergetic_cluster(struct ode_params* params,
    double target_compactness, double central_z, int solve_central_z,
    int ngrid, unsigned long long seed, double min_sep_factor,
    int remove_com, double* w0);

#endif
