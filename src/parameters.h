#ifndef PARAMETERS_H
#define PARAMETERS_H

// Parameter struct
typedef struct {
    char* name;
    char* value;
    char* description;
} tParameter;

void parse_parameter_file(const char *parfile);

void add_parameter(const char *name, const char *value, const char *description);
void add_parameter_i(const char *name, const int i, const char *value, const char *description);

void set_parameter_string(const char *name, const char *value);
void set_parameter_int(const char *name, const int i);
void set_parameter_double(const char *name, const double d);
void set_double_array(const char *name, int n, const double *a);
int set_if_unset_double(double *dst, double val);

int is_set_double(double x);
char *get_parameter_string(const char *name);
int get_parameter_int(const char *name);
double get_parameter_double(const char *name);
double get_parameter_double_i(const char* name, int i);
double* get_parameter_double_array(const char* name);
double* get_parameter_double_array_i(const char* name, const int i);
int get_parameter_array_count(const char* name);
double get_parameter_double_array_entry(const char* name, int index);
double get_binary_parameter_double_i(const char *par, int i);

#endif
