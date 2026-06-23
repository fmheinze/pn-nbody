#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>

struct ode_params;

void output_init(FILE** file_mass, FILE** file_pos, FILE** file_mom, FILE** file_spin, 
    FILE** file_energy, FILE** file_merger, struct ode_params* ode_params) ;
void output_write_timestep(FILE* file_pos, FILE* file_mom, FILE* file_spin, FILE* file_energy, 
    struct ode_params* ode_params, double* w, double t) ;
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
