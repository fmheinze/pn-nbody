#ifndef PARAMETERS_H
#define PARAMETERS_H

// Parameter struct
typedef struct {
    char* name;
    char* value;
    char* description;
} tParameter;

void parse_parameter_file(const char *parfile);
void make_parameter(const char *name, const char *value, const char *description);
void free_parameters();
void set_parameter(const char *name, const char *value);
tParameter* find_parameter(const char *name, const int fatal);
void print_parameters(void);

void add_parameter(const char *name, const char *value, const char *description);
void add_parameter_i(const char *name, const int i, const char *value, const char *description);
void AppPar(const char *name, const char *value);
void set_parameter_string(const char *name, const char *value);
void set_parameter_int(const char *name, const int i);
void set_parameter_double(const char *name, const double d);
void set_double_array(const char *name, int n, const double *a);
char *get_parameter_string(const char *name);
int get_parameter_int(const char *name);
double get_parameter_double(const char *name);
double get_parameter_double_i(const char* name, int i);
double* get_parameter_double_array(const char* name);
double* get_parameter_double_array_i(const char* name, const int i);
int get_parameter_array_count(const char* name);
double get_parameter_double_array_entry(const char* name, int index);

#endif