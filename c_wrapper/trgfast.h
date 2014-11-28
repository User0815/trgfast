#include<stdbool.h>


typedef struct _trg_options trg_options; 

typedef double (*trg_background_function)(double, double);

typedef struct _trg_object trg_object;




struct _trg_options {
    double k_linear;
    double k_th;
    int Nq;
    double tolerance;
    double rkthresh;
    int timesteps; 
    bool linear;
    int growth_length; 
    int verbosity;
    int dump_dir;
    bool want_ps;
    bool want_raw;
    bool want_im;
    bool want_diff;
    bool want_a;
    bool want_kdq;
    double scale_a;
    double curv_limit;
    double slope_limit;
    bool symmetry;
    int extra;
    double offset;

    int eta_length;
    bool scale_dependent_background;
}; 

struct _trg_object {
    trg_options *options;
    double *ps_in, **ps_out;
    int k_length_in, *k_length_out;
    trg_background_function O21, O22;    
    double *eta;
    double *z;
    int z_length;
};

trg_object* trg_init(
            double *z, int z_length,
            double* ps_in, int k_length,
            trg_background_function O21,
            trg_background_function O22
        );

void trg_free(trg_object* object);

double* trg_concatenate(int size, const double *k, const double *ps);

void trg_extract_spectrum(trg_object *trg, int zi, int species, double **k, double **ps);

void trg_evolution(trg_object* trg);

// Declaration of the Fortran function
extern void time_evolution_c(
                double* eta, int* eta_len, 
                double* ps_in, int* ps_len, 
                double* O_eta, int* O_eta_len,
                double* O_k, int* O_k_len,
                double* OmegaBulk, int* OBulk_len, 
                double* ps_out, double* growth_out, double* opts
            );
