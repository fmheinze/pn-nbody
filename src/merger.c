/**
 * @file merger.c
 * @brief Routines for merging bodies in the simulations.
 *
 * Routines for merging bodies in the simulations, such as finding merger pairs, updating the
 * ode_params and the state vector after merger, as well as prescriptions for the properties
 * of the merger remnant.
 */

#include "eom.h"
#include "output.h"
#include "stdio.h"
#include <math.h>
#include <float.h>

#define MERGE_FACTOR 2.0


/**
 * @brief Finds the closest merging pair.
 *
 * Finds the closest active pair satisfying:
 * 
 * r_ij < MERGE_FACTOR * (m_i + m_j)
 *
 * @param[in]   params     Parameter struct containing general information about the system
 * @param[in]   w          State vector w = [positions, momenta]
 * @param[out]  i_merge    Index of the first merger object
 * @param[out]  j_merge    Index of the second merger object
 * @return 1 if such a pair is found, 0 otherwise. The pair is written to *i_merge and *j_merge.
 */
int find_merger_pair(struct ode_params *params, double *w, int *i_merge, int *j_merge)
{
    int found = 0;
    double best_dist2 = DBL_MAX;

    const int n_bodies = params->num_bodies_initial;
    const int num_dim = params->num_dim;

    for (int i = 0; i < n_bodies; i++) {
        if (!params->active[i])
            continue;

        for (int j = i + 1; j < n_bodies; j++) {
            if (!params->active[j])
                continue;

            const double mi = params->masses[i];
            const double mj = params->masses[j];
            const double M_tot = mi + mj;

            const double r_merge = MERGE_FACTOR * M_tot;
            const double r_merge2 = r_merge * r_merge;

            double dist2 = 0.0;

            for (int n = 0; n < num_dim; n++) {
                const double dx = w[num_dim * i + n] - w[num_dim * j + n];
                dist2 += dx * dx;
            }

            if (dist2 < r_merge2 && dist2 < best_dist2) {
                best_dist2 = dist2;
                *i_merge = i;
                *j_merge = j;
                found = 1;
            }
        }
    }
    return found;
}


/**
 * @brief Merges active bodies i and j.
 *
 * Current rough prescription:
 *   M_f = m_i + m_j
 *   x_f = mass-weighted COM position
 *   p_f = p_i + p_j
 *
 * Inactive slot j is kept finite internally:
 *   x_j = x_f
 *   p_j = 0
 * 
 * @param   params        Parameter struct containing general information about the system
 * @param   w             State vector w = [positions, momenta]
 * @param   i             Index of the first merger object
 * @param   j             Index of the second merger object
 * @param   t             Current time
 * @param   file_merger   Pointer to the merger output file
 */
void merge_pair(
    struct ode_params *params,
    double *w,
    int i,
    int j,
    double t,
    FILE *file_merger
) {
    const int num_dim = params->num_dim;
    const int array_half = num_dim * params->num_bodies_initial;

    if (!params->active[i] || !params->active[j])
        return;

    if (i == j)
        return;

    const double mi = params->masses[i];
    const double mj = params->masses[j];
    const double m_remnant = mi + mj;

    if (m_remnant <= 0.0)
        return;

    const long long id_i = params->body_id[i];
    const long long id_j = params->body_id[j];
    const long long id_remnant = params->next_body_id++;

    const int gen_i = params->generation[i];
    const int gen_j = params->generation[j];
    const int gen_remnant = 1 + (gen_i > gen_j ? gen_i : gen_j);

    double x_i_old[num_dim];
    double x_j_old[num_dim];
    double x_rem[num_dim];

    double p_i_old[num_dim];
    double p_j_old[num_dim];
    double p_rem[num_dim];

    double r2 = 0.0;

    /*
     * Save old parent states and compute remnant state.
     */
    for (int n = 0; n < num_dim; n++) {
        const int xi = num_dim * i + n;
        const int xj = num_dim * j + n;

        const int pi = array_half + num_dim * i + n;
        const int pj = array_half + num_dim * j + n;

        x_i_old[n] = w[xi];
        x_j_old[n] = w[xj];

        p_i_old[n] = w[pi];
        p_j_old[n] = w[pj];

        const double dx = x_i_old[n] - x_j_old[n];
        r2 += dx * dx;

        x_rem[n] = (mi * x_i_old[n] + mj * x_j_old[n]) / m_remnant;
        p_rem[n] = p_i_old[n] + p_j_old[n];
    }

    const double r_ij = sqrt(r2);

    /*
     * Write merger event before overwriting bookkeeping.
     */
    if (file_merger != NULL) {
        output_write_merger_event(
            file_merger,
            t,
            params,
            i,
            j,
            i,
            id_i,
            id_j,
            id_remnant,
            gen_i,
            gen_j,
            gen_remnant,
            mi,
            mj,
            m_remnant,
            x_i_old,
            x_j_old,
            x_rem,
            p_i_old,
            p_j_old,
            p_rem,
            r_ij
        );
    }

    for (int n = 0; n < num_dim; n++) {
        const int xi = num_dim * i + n;
        const int xj = num_dim * j + n;

        const int pi = array_half + num_dim * i + n;
        const int pj = array_half + num_dim * j + n;

        w[xi] = x_rem[n];
        w[pi] = p_rem[n];

        w[xj] = x_rem[n];
        w[pj] = 0.0;
    }

    params->masses[i] = m_remnant;
    params->masses[j] = 0.0;

    params->active[i] = 1;
    params->active[j] = 0;
    params->num_active--;

    params->body_id[i] = id_remnant;
    params->body_id[j] = -1;

    params->generation[i] = gen_remnant;
    params->generation[j] = -1;
}


/**
 * @brief Repeatedly searches for merger pairs and merges them.
 *
 * This handles triple/quadruple situations by merging one pair,
 * then re-running the search on the updated system.
 * 
 * @param   params     Parameter struct containing general information about the system
 * @param   w          State vector w = [positions, momenta]
 * @return 1 if at least one merger happened, 0 otherwise.
 */
int test_and_merge_bodies(struct ode_params *params, double *w, double t, FILE *file_merger)
{
    int merged_any = 0;

    while (1) {
        int i_merge = -1;
        int j_merge = -1;

        const int found = find_merger_pair(params, w, &i_merge, &j_merge);

        if (!found)
            break;

        merge_pair(params, w, i_merge, j_merge, t, file_merger);
        merged_any = 1;
    }

    return merged_any;
}