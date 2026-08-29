#ifndef PN_BINARY_H
#define PN_BINARY_H


double pn_binary_reduced_hamiltonian(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn);

double pn_binary_reduced_hamiltonian_dx(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn);

double pn_binary_reduced_hamiltonian_dpr(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn);

double pn_binary_reduced_hamiltonian_dj(double x, double pr_hat, double j, double nu,
    int use_1pn, int use_2pn);

double pn_binary_turning_hamiltonian(double x, double j, double nu,
    int use_1pn, int use_2pn);

double pn_binary_solve_circular_j(double x, double nu, int use_1pn, int use_2pn);

double pn_binary_solve_eccentric_j(double xp, double xa, double nu,
    int use_1pn, int use_2pn);

double pn_binary_solve_pr_hat_abs(double x, double j, double reduced_energy, double nu,
    int use_1pn, int use_2pn);

double pn_binary_dh_dpr_coeff_at_zero(double x, double j, double nu,
    int use_1pn, int use_2pn);

int pn_binary_solve_turning_points(double reduced_energy, double j, double nu,
    int use_1pn, int use_2pn, double current_x,
    double *x_pericenter, double *x_apocenter);

#endif
