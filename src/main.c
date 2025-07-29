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


void read_command_line(int argc, char** argv)
{
    if (argc != 2)
        errorexit("Usage: pn_nbody /path/to/parfile.par");

    /* Determine the name of the parameter file and the output directory,
    the output folder has the same name as the parfile (without .par) */

    // Build parfile path
    char* parfile = (char*) calloc(strlen(argv[1]) + 40, sizeof(char));
    strcpy(parfile, argv[1]);
    if (!strstr(parfile, ".par")) strcat(parfile, ".par");

    // Extract base name (without .par)
    char* base = basename(parfile);
    char* dot = strstr(base, ".par");
    if (dot) *dot = '\0';

    // Get executable path
    char base_path[] = {"/Users/fheinze/Desktop/pn-nbody"};

    // Create output directory
    char outdir[PATH_MAX];
    snprintf(outdir, sizeof(outdir), "%s/output/%s", base_path, base);
    if (mkdir(outdir, 0755) != 0)
        if (errno != EEXIST)
            errorexit("Could not create output directory");
    
    /* Add outdir and parfile as parameters, which should be the first parameters
    that lead to the creation of the parameter database */
    add_parameter("outdir", outdir, "output directory");
    add_parameter("parfile", parfile, "name of parameter file given in command line");
    
    free(parfile);
}


int main(int argc, char** argv) {
    print_divider();
    printf("Welcome to pn_nbody\n");
    print_divider();

    // Parameters
    printf("Initial parameters:\n");
    read_command_line(argc, argv);
    parse_parameter_file(get_parameter_string("parfile"));
    initialize_parameters();
    print_divider();

    // Initial values
    struct ode_params params = initialize_ode_params();
    double* w = initialize_state_vector(&params);

    // Run simulation
    printf("Starting simulation...\n"); 
    ode_integrator(w, rhs_pn_threebody, &params);
    printf("\n"); 
    print_divider();

    // Finalize
    printf("The simulation has finished! Thanks for using pn_nbody!\n");
    free_vector(w);
    free_ode_params(&params);
    return 0;
}
