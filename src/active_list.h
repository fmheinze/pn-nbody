#ifndef ACTIVE_LIST_H
#define ACTIVE_LIST_H

#include "eom.h"


typedef struct {
    int num_bodies_initial;
    int num_active;
    int *ids;
} ActiveList;


void active_list_init(ActiveList *list, const struct ode_params *ode_params);
void active_list_refresh(ActiveList *list, const struct ode_params *ode_params);
void active_list_free(ActiveList *list);

#endif
