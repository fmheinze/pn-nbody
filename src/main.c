/**
 * @file main.c
 * @brief Main program entry point and high-level control flow.
 *
 * High-level control flow for reading the command line, initializing parameters and initial
 * conditions and performing the numerical ODE integration.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <limits.h>
#include "integration.h"
#include "utils.h"
#include "pn_eom.h"
#include "parameters.h"
#include "initialize.h"
#include "pn_eom_hamiltonians.h"


#ifdef __APPLE__
  #include <mach-o/dyld.h>
#endif


/**
 * @brief Parses command-line arguments and creates output directory.
 * 
 * Parses command-line arguments, creates the project output directory, and initializes the 
 * parameter database by adding the first parameters (parameter file name and output directory)
 * 
 * @param[in]   argc    Number of command-line arguments
 * @param[in]   argv    Array of command-line argument strings
 */
void read_command_line(int argc, char** argv)
{
    if (argc != 2)
        errorexit("Usage: pn_nbody /path/to/parfile.par");

    // Build parfile path (append .par if missing)
    char* parfile = calloc(strlen(argv[1]) + 5, 1);
    if (!parfile) errorexit("Out of memory");
    strcpy(parfile, argv[1]);
    if (!strstr(parfile, ".par")) strcat(parfile, ".par");

    // Extract base name (without .par)
    char parfile_copy[PATH_MAX];
    snprintf(parfile_copy, sizeof(parfile_copy), "%s", parfile);
    char *base_ptr = basename(parfile_copy);
    char *dot = strstr(base_ptr, ".par");
    if (dot) *dot = '\0';

    // Get executable directory and then its parent
    char exe_dir[PATH_MAX];
    if (get_executable_dir(exe_dir) != 0) {
        fprintf(stderr, "Error: could not resolve executable directory: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    char parent_dir[PATH_MAX];
    snprintf(parent_dir, sizeof(parent_dir), "%s/..", exe_dir);

    char parent_real[PATH_MAX];
    if (!realpath(parent_dir, parent_real)) {
        fprintf(stderr, "Error: realpath(%s) failed: %s\n", parent_dir, strerror(errno));
        exit(EXIT_FAILURE);
    }

    // Ensure <parent>/output exists, then create <parent>/output/<base>
    char output_root[PATH_MAX];
    snprintf(output_root, sizeof(output_root), "%s/output", parent_real);
    mkdir_or_die(output_root, 0755);

    char outdir[PATH_MAX];
    snprintf(outdir, sizeof(outdir), "%s/%s", output_root, base_ptr);
    mkdir_or_die(outdir, 0755);

    // Pass parameters onward
    add_parameter("outdir", outdir, "output directory");
    add_parameter("parfile", parfile, "name of parameter file given in command line");

    free(parfile);
}


/**
 * @brief Program entry point.
 * 
 * Parses command-line arguments, initializes the simulation, executes the numerical ODE 
 * integration, and performs the final cleanup.
 * 
 * @param[in]   argc    Number of command-line arguments
 * @param[in]   argv    Array of command-line argument strings
 * @return EXIT_SUCCESS on successful completion,
 *         EXIT_FALIURE if an error occurs during execution.
 */
int main(int argc, char** argv) 
{
    print_divider();
    printf("Welcome to pn-nbody\n");
    print_divider();

    // Read an initialize parameters
    printf("Specified parameters:\n");
    read_command_line(argc, argv);
    parse_parameter_file(get_parameter_string("parfile"));
    initialize_parameters();
    struct ode_params params = initialize_ode_params();
    print_divider();

    // Initialize state vector with specified initial conditions
    double* w = initialize_state_vector(&params);

    // Run simulation
    printf("Running simulation...\n"); 
    if (get_parameter_int("impulse_method")) {
        ode_integrator_impulse(w, rhs_pn_nbody, compute_dUTT4_dx, &params);
    }
    else {
        ode_integrator(w, rhs_pn_nbody, &params);
    }
    printf("\n"); 
    print_divider();

    // Finalize
    printf("The simulation has finished! Thanks for using pn-nbody!\n");
    free_vector(w);
    free_ode_params(&params);
    return EXIT_SUCCESS;
}
