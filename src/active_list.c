/**
 * @file active_list.c
 * @brief Compact list of active bodies for faster N-body loops.
 *
 * List of active bodies to avoid loops over inactive bodies. Storage is sized once for all body
 * slots. Refreshing the list only rewrites the active indices, so repeated RHS evaluations do not
 * allocate memory.
 */

#include <stdlib.h>
#include "active_list.h"
#include "utils.h"


/**
 * @brief Allocate storage and build the active-body list.
 *
 * The function counts active flags directly instead of trusting ode_params->num_active.
 */
void active_list_init(ActiveList *list, const struct ode_params *ode_params)
{
    if (list == NULL)
        errorexit("active_list_init received NULL list");
    if (ode_params == NULL)
        errorexit("active_list_init received NULL ode_params");
    if (ode_params->num_bodies_initial < 0)
        errorexit("num_bodies_initial must be non-negative");
    if (ode_params->num_bodies_initial > 0 && ode_params->active == NULL)
        errorexit("ode_params->active is NULL");

    list->num_bodies_initial = ode_params->num_bodies_initial;
    list->num_active = 0;
    list->ids = NULL;

    if (list->num_bodies_initial > 0) {
        list->ids = (int *)malloc((size_t)list->num_bodies_initial * sizeof(*list->ids));
        if (list->ids == NULL)
            errorexit("Memory allocation failed for ActiveList ids");
    }

    active_list_refresh(list, ode_params);
}


/**
 * @brief Refresh active indices without allocating or freeing memory.
 */
void active_list_refresh(ActiveList *list, const struct ode_params *ode_params)
{
    if (list == NULL)
        errorexit("active_list_refresh received NULL list");
    if (ode_params == NULL)
        errorexit("active_list_refresh received NULL ode_params");
    if (list->num_bodies_initial != ode_params->num_bodies_initial)
        errorexit("ActiveList size does not match ode_params");
    if (list->num_bodies_initial > 0 && list->ids == NULL)
        errorexit("ActiveList storage is not initialized");
    if (list->num_bodies_initial > 0 && ode_params->active == NULL)
        errorexit("ode_params->active is NULL");

    list->num_active = 0;
    for (int a = 0; a < list->num_bodies_initial; ++a) {
        if (ode_params->active[a])
            list->ids[list->num_active++] = a;
    }
}


/**
 * @brief Release memory owned by the active-body list.
 */
void active_list_free(ActiveList *list)
{
    if (list == NULL)
        return;

    free(list->ids);
    list->ids = NULL;
    list->num_active = 0;
    list->num_bodies_initial = 0;
}
