#ifndef MERGER_H
#define MERGER_H

struct ode_params;
struct OutputContext;

int test_and_merge_bodies(struct ode_params *params, double *w, double t,
    struct OutputContext *output);

#endif
