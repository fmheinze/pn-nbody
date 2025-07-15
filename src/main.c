#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
    the output folder has the same name and location as the parfile (without .par) */
    char* parfile = (char*) calloc(sizeof(char), strlen(argv[1])+40);
    char* outdir  = (char*) calloc(sizeof(char), strlen(argv[1])+40);
    strcpy(parfile, argv[1]);
    if (!strstr(parfile, ".par")) strcat(parfile, ".par");
    strcpy(outdir, parfile);
    *strstr(outdir, ".par") = '\0';
    
    /* Add outdir and parfile as parameters, which should be the first parameters
    that lead to the creation of the parameter database */
    add_parameter("outdir", outdir, "output directory");
    add_parameter("parfile", parfile, "name of parameter file given in command line");
    
    free(parfile);
    free(outdir);
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
