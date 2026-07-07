#ifndef ACTIVE_LIST_H
#define ACTIVE_LIST_H

#include "eom.h"

/**
 * @brief Compact list of currently active body indices.
 *
 * The code often loops over num_bodies_initial and skips inactive entries.
 * This helper builds a dense list of active body slots once per RHS/cache build,
 * so expensive PN loops can later iterate only over active bodies.
 */
typedef struct {
    int num_bodies_initial;
    int num_active;
    int *ids;
} ActiveList;

/**
 * @brief Build the active-body list from ode_params->active.
 *
 * The function counts active flags directly instead of trusting
 * ode_params->num_active. This makes the helper robust during future refactors.
 */
void active_list_init(ActiveList *list, const struct ode_params *ode_params);

/**
 * @brief Release memory owned by the active-body list.
 */
void active_list_free(ActiveList *list);

#endif
