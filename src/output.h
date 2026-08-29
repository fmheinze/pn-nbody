#ifndef OUTPUT_H
#define OUTPUT_H

#include <stddef.h>
#include <stdio.h>

struct ode_params;
struct OrbitOutput;


typedef struct OutputContext {
    FILE *file_mass;
    FILE *file_pos;
    FILE *file_mom;
    FILE *file_vel;
    FILE *file_spin;
    FILE *file_energy;
    FILE *file_angular_momentum;
    FILE *file_angular_momentum_vector;
    FILE *file_orbital_angular_momentum;
    FILE *file_orbital_angular_momentum_vector;
    FILE *file_spin_angular_momentum;
    FILE *file_spin_angular_momentum_vector;
    FILE *file_angular_momentum_com;
    FILE *file_angular_momentum_vector_com;
    FILE *file_orbital_angular_momentum_com;
    FILE *file_orbital_angular_momentum_vector_com;
    FILE *file_merger;

    struct OrbitOutput *orbits;
    size_t num_orbits;
} OutputContext;


void output_init(OutputContext *output, struct ode_params *ode_params);

void output_write_timestep(OutputContext *output, struct ode_params *ode_params, double *w,
    double t);

void output_follow_merger(OutputContext *output, int parent_a, int parent_b, int remnant);

void output_close(OutputContext *output);

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
);

#endif
