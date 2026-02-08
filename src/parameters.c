/**
 * @file parameters.c
 * @brief Functions for managing the parameter database.
 *
 * Functions for managing the parameter database, e.g. parsing parameter files, creating,
 * querying and modifying parameters.
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "parameters.h"


// Parameter database
tParameter* pdb;

// Number of parameters
int npdb, npdbmax = 1000;


// ------------------------------------------------------------------------------------------------
// Helper functions for the parameter database
// ------------------------------------------------------------------------------------------------

// Find parameter in the parameter database
static tParameter* find_parameter(const char* name, const int fatal) {
    if (!name) 
        errorexit("find_parameter: called without parameter name");

    for (int i = 0; i < npdb; i++)
        if (!strcmp(pdb[i].name, name))
        return &pdb[i];

    if (fatal) {
        printf("Could not find parameter \"%s\"\n", name);
        errorexit("this one is required!");
    }
    return 0;
}


// Create a new parameter in parameter database, merge if already there
static void make_parameter(const char *name, const char *value, const char *description) {
    static int firstcall = 1;
    tParameter* p;

    // If this function is called for the first time, allocate memory for parameter database
    if (firstcall) {
        firstcall = 0;
        pdb = (tParameter *) calloc(npdbmax, sizeof(tParameter));
        if (!pdb) errorexit("make_parameter: out of memory");
        npdb = 0;
    }

    p = find_parameter(name, 0);

    // If the parameter cannot be found in the dat base, create it
    if (!p) {
        p = &pdb[npdb++];
        p->name = (char *) calloc(strlen(name)+1, sizeof(char));
        p->value = (char *) calloc(strlen(value)+1, sizeof(char));
        strcpy(p->name, name);
        strcpy(p->value, value);
    } 
    // If it is already there, update the description
    else {
        free(p->description);
    }
    p->description = (char *) calloc(strlen(description)+1, sizeof(char));
    strcpy(p->description, description);

    if (npdb >= npdbmax) {
        printf("The maximum number of %d parameters has been exceeded!", npdbmax);
        errorexit("There is no more space for new parameters");
    }
}


// Set parameter in the parameter database
static void set_parameter(const char* name, const char* value) {
    tParameter* p = find_parameter(name, 1);
    if (!p)
        errorexit("set_parameter: parameter not found or could not be created.");

    free(p->value);
    p->value = strdup(value);
    if (!p->value)
        errorexit("set_parameter: out of memory while duplicating value string.");
}


// ------------------------------------------------------------------------------------------------
// Functions for external calls
// ------------------------------------------------------------------------------------------------

/**
 * @brief Parses a given parameter file.
 * 
 * Parses a given parameter file and writes the parameter values to the parameter database.
 * If a line starts with a #, it is seen as a comment. Multiple white spaces and " are collapsed
 * to a single space.
 * 
 * @param[in]   parfile     Path to the parameter file
 */
void parse_parameter_file(const char *parfile)
{
    FILE* fp;
    int c, i, j;
    int nbuffer, lpar, lval;
    char *buffer, *par, *val;

    // Read file into memory, add one space at end and beginning
    fp = fopen(parfile, "r");
    if (!fp) {
        printf("Could not open parameter file \"%s\"\n", parfile);
        errorexit("");
    }
    buffer = 0;
    for (i = nbuffer = 0;; i++) {
        if (i >= nbuffer-2) {
            if (nbuffer > 1000000) 
                errorexit("Sanity forbids parameter files bigger than 1MB");
            buffer = (char *) realloc(buffer, sizeof(char)*(nbuffer += 1000));
            if (!buffer) errorexit("Out of memory while reading parameter file.");
        }
        if (i == 0) buffer[i++] = ' ';
        if ((c = fgetc(fp)) == EOF) break;
        buffer[i] = c;
    }
    fclose(fp);
    buffer[i++] = ' ';
    buffer[i] = '\0';
    nbuffer = strlen(buffer);

    // Replace comments and quotes with white spaces
    for (i = 0; i < nbuffer; i++) { 
        if (buffer[i] == '#') 
            while (i < nbuffer && buffer[i] != '\n')
                buffer[i++] = ' ';
        if (buffer[i] == '"')
            buffer[i] = ' ';
    }

    // Collapse all white space to a single space
    for (i = j = 1; i < nbuffer; i++) {
        if (!isspace(buffer[i])) buffer[j++] = buffer[i];
        else if (!isspace(buffer[j-1])) buffer[j++] = ' ';
    }
    buffer[j] = '\0';
    nbuffer = strlen(buffer);

    // Remove spaces around "="
    for (i = j = 1; i < nbuffer; i++) {
        if (buffer[i] != ' ' || (buffer[i-1] != '=' && buffer[i+1] != '=')) 
            buffer[j++] = buffer[i];
    }
    buffer[j-1] = '\0';
    nbuffer = strlen(buffer);

    // Split parameter names and values by inserting zeros
    for (i = 1; i < nbuffer; i++) {
        if (buffer[i] == '=') {
            buffer[i] = '\0';
            for (j = i-1; j > 0 && buffer[j-1] != ' '; j--)
            ;
            buffer[j-1] = '\0';
        }
    }
    
    // Loop over all parameter/value pairs and make/set the parameter values
    for (i = 1; i < nbuffer; i += lpar + lval + 2) {
        par = buffer+i;
        lpar = strlen(par);
        val = par + lpar + 1;
        lval = strlen(val);

        if (!find_parameter(par, 0))
            make_parameter(par, val, "not in library of initialized parameters");
        else
            set_parameter(par, val);
    }  
    free(buffer);
}


// Free parameter database memory
void free_parameters() {
    if (!pdb) return;

    for (int i = 0; i < npdb; i++) {
        free(pdb[i].name);
        free(pdb[i].value);
        free(pdb[i].description);
    }

    free(pdb);
    pdb = NULL;
    npdb = 0;
}


/* --- Creation functions --- */

void add_parameter(const char* name, const char* value, const char* description) {
    make_parameter(name, value, description);
    printf("%-30s  =  %s\n", name, get_parameter_string(name));
}


void add_parameter_i(const char* name, const int i, const char* value, const char* description) {
    char new[100];
    sprintf(new, "%s%d", name, i);
    make_parameter(new, value, description);
    printf("%-30s  =  %s\n", new, get_parameter_string(new));
}


/* --- Assignment functions --- */

void set_parameter_string(const char* name, const char *value) {
    set_parameter(name, value);
}


void set_parameter_int(const char* name, const int i) {
    char value[100];
    sprintf(value, "%d", i);
    set_parameter(name, value);
}


void set_parameter_double(const char* name, const double d) {
    char value[100];
    sprintf(value, "%.20e", d);
    set_parameter(name, value);
}


void set_double_array(const char *name, int n, const double *a) {
    if (n <= 0 || !a) return;

    // Estimate total string size
    int estimated_len = n * 25;
    char *value = malloc(estimated_len);
    if (!value) {
        perror("malloc");
        return;
    }

    value[0] = '\0';
    char buf[32];

    for (int i = 0; i < n; ++i) {
        snprintf(buf, sizeof(buf), "%.16e", a[i]);
        strcat(value, buf);
        if (i < n - 1) strcat(value, " ");
    }

    tParameter *p = find_parameter(name, 1);
    if (!p) {
        free(value);
        errorexit("set_double_array: parameter allocation failed.");
    }
    free(p->value);
    p->value = value;
}


int set_if_unset_double(double *dst, double val) {
    if (!is_set_double(*dst) && isfinite(val)) {
        *dst = val;
        return 1;
    }
    return 0;
}


/* --- Query functions --- */

int is_set_double(double x) { 
    return x >= 0.0; 
}


char *get_parameter_string(const char* name) {
    tParameter* p = find_parameter(name, 1);
    return p->value; 
}


int get_parameter_int(const char* name) {
    tParameter* p = find_parameter(name, 1);
    return atoi(p->value);
}


double get_parameter_double(const char* name) {
    tParameter* p = find_parameter(name, 1);
    return atof(p->value);
}


double get_parameter_double_i(const char* name, int i) {
    char new[100];
    sprintf(new, "%s%d", name, i);
    return get_parameter_double(new);
}


double* get_parameter_double_array(const char* name) {
    char* param_string = get_parameter_string(name);

    // Make a copy to safely tokenize
    char* copy = strdup(param_string);
    if (!copy) {
        perror("strdup");
        return NULL;
    }

    // First pass: count tokens
    int count = 0;
    char* token = strtok(copy, " ");
    while (token) {
        count++;
        token = strtok(NULL, " ");
    }
    free(copy);

    // Allocate array with enough memory
    double* array = malloc(sizeof(double) * count);
    if (!array) {
        perror("malloc");
        return NULL;
    }

    // Second pass: parse doubles
    copy = strdup(param_string);
    if (!copy) {
        free(array);
        perror("strdup");
        return NULL;
    }
    token = strtok(copy, " ");
    int i = 0;
    while (token && i < count) {
        char* endptr;
        array[i] = strtod(token, &endptr);
        if (*endptr != '\0') {
            fprintf(stderr, "Invalid double: \"%s\" in parameter \"%s\"\n", token, name);
            free(copy);
            free(array);
            return NULL;
        }
        i++;
        token = strtok(NULL, " ");
    }

    free(copy);
    return array;
}


double* get_parameter_double_array_i(const char* name, const int i) {
    char new[100];
    sprintf(new, "%s%d", name, i);
    return get_parameter_double_array(new);
}


int get_parameter_array_count(const char* name) {
    char* param_string = get_parameter_string(name);

    // Make a copy to safely tokenize
    char* copy = strdup(param_string);
    if (!copy) {
        perror("strdup");
        return -1;
    }

    // Count tokens
    int count = 0;
    char* token = strtok(copy, " ");
    while (token) {
        count++;
        token = strtok(NULL, " ");
    }
    free(copy);
    return count;
}

double get_parameter_double_array_entry(const char* name, int index) {
    int count = get_parameter_array_count(name);
    double *array = get_parameter_double_array(name);
    if (!array)
        errorexit("get_double_entry: failed to retrieve array.");

    if (index < 0 || index >= count) {
        fprintf(stderr, "Index %d out of bounds (size %d) for parameter \"%s\"\n", index, count, name);
        free(array);
        errorexit("get_double_entry: index out of bounds");
    }

    double result = array[index];
    free(array);
    return result;
}


double get_binary_parameter_double_i(const char *par, int i) {
    char key[64];
    if (i == 0)
        snprintf(key, sizeof(key), "binary_%s", par);
    else
        snprintf(key, sizeof(key), "binary%d_%s", i, par);
    return get_parameter_double(key);
}
