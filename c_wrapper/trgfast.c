// @file trgfast.c
// A C wrapper that calls the Fortran routine.
// @author Adrian Vollmer
// @date 25.11.2014


#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>

#include"trgfast.h"


#define len_options 23

const trg_options trg_default_options = { .001, 0.01, 200, .01, .001, 30,
    false, 100, 2, 0, false, false, false, false, false, false, 1., 10.,
    20., true, 2, 1e-4 };



// Some helpful functions

trg_object* trg_init(
            double *z, int z_length,
            double *ps_in, int k_length,
            trg_background_function O21,
            trg_background_function O22
        ){

    trg_object *trg = malloc(sizeof(trg_object));
    trg_options *options = malloc(sizeof(trg_options));
    trg->options = options;
    memcpy(trg->options, &trg_default_options, sizeof(trg_options));

    
    trg->O21 = O21;
    trg->O22 = O22;

    trg->ps_in = ps_in;
    trg->k_length_in = k_length;
    trg->ps_out = NULL;
    trg->z = z;
    trg->z_length = z_length;

    trg->options->scale_dependent_background = false;
    trg->options->eta_length=50;

    return trg;
}

void trg_free(trg_object* trg){
    int i;
    free(trg->options);
    free(trg->ps_in);
    for (i=0; i<trg->z_length-1; i++){
        free(trg->ps_out[i]);
    }
    free(trg->ps_out);
    free(trg);
}

double* trg_concatenate(int size, const double *k, const double *ps){

    double *result = malloc(sizeof(double));
    result = malloc(2*size*sizeof(double));

    memcpy(result, k, size*sizeof(double));
    memcpy(result+size, ps, size*sizeof(double));

    return result;
}


// Some private functions

void prepare_background_functions(
            trg_object *trg,
            double** O_eta,  int *O_eta_len,
            double** O_k,    int *O_k_len,
            double** OmegaBulk, int *OBulk_len
        ){

    int i, j;
    double kmin = 0.001, kmax=1000;

    *O_eta_len = trg->options->eta_length;
    *O_k_len = 1;
    if (trg->options->scale_dependent_background) 
        *O_k_len = 100;
    *OBulk_len = 2**O_k_len**O_eta_len;

    *O_k = malloc(*O_k_len * sizeof(double));
    
    if (*O_k_len > 1){
        for (i=0; i<*O_k_len; i++){
            (*O_k)[i] = exp(log(kmin)+log(kmax-kmin)*i/(*O_k_len-1));
        }
    }

    *O_eta = malloc(*O_eta_len * sizeof(double));
    *OmegaBulk = malloc(*OBulk_len * sizeof(double));
    for (i=0; i<*O_eta_len; i++){
        (*O_eta)[i] = i*log(1+trg->z[0])/(*O_eta_len-1);
        for (j=0; j<*O_k_len; j++){
            (*OmegaBulk)[i+*O_eta_len*2*j] =
                    trg->O21(exp((*O_eta)[i])/(1+trg->z[0]), (*O_k)[j]);
            (*OmegaBulk)[i+*O_eta_len*(2*j+1)] =
                    trg->O22(exp((*O_eta)[i])/(1+trg->z[0]), (*O_k)[j]);
        }
    }
}

void prepare_ps_out(trg_object *trg, double *ps_out, int size){
    int i, counter = 0;
    trg->k_length_out = malloc((trg->z_length-1)*sizeof(int));
    trg->ps_out = malloc((trg->z_length-1)*sizeof(double*));
    for (i=0; i<trg->z_length-1; i++) {
        trg->k_length_out[i] = size/(trg->z_length-1);
        trg->ps_out[i] = malloc(trg->k_length_out[i]*sizeof(double));
        memcpy(trg->ps_out[i], ps_out + counter, trg->k_length_out[i] * sizeof(double));
        counter += trg->k_length_out[i];
        trg->k_length_out[i] /= 4;
    }
}


void prepare_options_array(trg_object* trg, double(* options)[len_options]){
    (*options)[0] = -1; // no default options
    (*options)[1] = trg->options->k_linear;
    (*options)[2] = trg->options->k_th;
    (*options)[3] = trg->options->Nq;
    (*options)[4] = trg->options->tolerance;
    (*options)[5] = trg->options->rkthresh;
    (*options)[6] = trg->options->timesteps; 
    (*options)[7] = trg->options->linear>0;
    (*options)[8] = trg->options->growth_length; 
    (*options)[9] = trg->options->verbosity;
    (*options)[10] = trg->options->dump_dir;
    (*options)[11] = trg->options->want_ps>0;
    (*options)[12] = trg->options->want_raw>0;
    (*options)[13] = trg->options->want_im>0;
    (*options)[14] = trg->options->want_diff>0;
    (*options)[15] = trg->options->want_a>0;
    (*options)[16] = trg->options->want_kdq>0;
    (*options)[17] = trg->options->scale_a;
    (*options)[18] = trg->options->curv_limit;
    (*options)[19] = trg->options->slope_limit;
    (*options)[20] = trg->options->symmetry;
    (*options)[21] = trg->options->extra;
    (*options)[22] = trg->options->offset;
}


void trg_extract_spectrum(trg_object *trg, int zi, int species, double **k, double **ps){ 

        //TODO if zi>z_length: error
    int k_len = trg->k_length_out[zi];
    *k = malloc(k_len*sizeof(double));
    *ps = malloc(k_len*sizeof(double));
    memcpy(*k, trg->ps_out[zi], k_len*sizeof(double));
    memcpy(*ps, (double*)(trg->ps_out[zi])+species*k_len, k_len*sizeof(double));
}

// Finally, the interface

void trg_evolution(trg_object* trg){

        // Declare arrays and variables.
    double options[len_options], growth_out[1000];
    double *ps_in = NULL,  *O_eta = NULL, *O_k = NULL, *OmegaBulk = NULL;
    int eta_len, ps_len, O_eta_len, O_k_len, OBulk_len, i;


    ps_in = trg->ps_in;
    ps_len = 2*trg->k_length_in;

    
    prepare_background_functions(trg, &O_eta, &O_eta_len, &O_k,
            &O_k_len, &OmegaBulk, &OBulk_len);

    prepare_options_array(trg, &options);

        //Convert redshifts to eta
    eta_len=trg->z_length;
    double *eta = malloc(eta_len*sizeof(double));
    memcpy(eta, trg->z, eta_len*sizeof(double));
    for (i=0; i<eta_len; i++) eta[i] = log((1+trg->z[0])/(1+trg->z[i]));

    double *ps_out = malloc(4*trg->k_length_in*(eta_len-1)*sizeof(double));

        // Call Fortran function.
    time_evolution_c(eta, &eta_len, 
            ps_in, &ps_len,
            O_eta, &O_eta_len, 
            O_k, &O_k_len,
            OmegaBulk, &OBulk_len, 
            ps_out, 
            growth_out, 
            options);

    prepare_ps_out(trg, ps_out, 4*trg->k_length_in*(eta_len-1)); 

        // Clean up
    free(eta);
    free(O_eta);
    free(O_k);
    free(OmegaBulk);
    free(ps_out);
}

