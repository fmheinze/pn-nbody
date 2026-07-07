/**
 * @file active_list.c
 * @brief Compact list of active bodies for faster N-body loops.
 */

#include <stdlib.h>
#include "active_list.h"
#include "utils.h"


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

    for (int a = 0; a < ode_params->num_bodies_initial; ++a) {
        if (ode_params->active[a])
            list->num_active++;
    }

    if (list->num_active == 0)
        return;

    list->ids = (int *)malloc((size_t)list->num_active * sizeof(int));
    if (list->ids == NULL)
        errorexit("Memory allocation failed for ActiveList ids");

    int k = 0;
    for (int a = 0; a < ode_params->num_bodies_initial; ++a) {
        if (ode_params->active[a])
            list->ids[k++] = a;
    }
}


void active_list_free(ActiveList *list)
{
    if (list == NULL)
        return;

    free(list->ids);
    list->ids = NULL;
    list->num_active = 0;
    list->num_bodies_initial = 0;
}
