
#include"trgfast.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Define the background functions for LCDM
double OmegaM0 = .25, OmegaL0 = .75, h = .73;

double Omega21(double a, double k){
    return -1.5 * OmegaM0/(pow(a,3)*OmegaL0 + OmegaM0);
}

double Omega22(double a, double k){
    return 2. - (-2*pow(a,3)*OmegaL0 + OmegaM0)/(2.*(pow(a,3)*OmegaL0 + OmegaM0));
}


// Load a file with two columns into two arrays
int import(const char *filename, double **k, double **ps) {
    FILE *myfile;
    char line[1024];
    int i, file_length;
    myfile = fopen(filename, "r");
    if (NULL==myfile) { printf("File not found: %s\n", filename); return -1; };
    file_length=0;
    while(fscanf(myfile, "%[^\n]\n", line) != EOF) {
        if (line[0]!='#') {
            file_length++;
        }
    }
    rewind(myfile);
    *k = malloc(file_length * sizeof(double));
    *ps = malloc(file_length * sizeof(double));
    i = 0;
    while(fscanf(myfile, "%[^\n]\n", line) != EOF) {
        if (line[0]!='#') {
            sscanf(line, "%le %le\n", &(*k)[i], &(*ps)[i]);
            i++;
        }
    }
    if (NULL!=myfile) fclose(myfile);
    return file_length;
}


// Writes a 1D array into a 4-column file
void export(const char *filename, double *k, double *ps, int array_len) {
    FILE *myfile;
    int i;
    myfile = fopen(filename, "w");
    for (i=0; i<array_len; i++) {
        if (ps[i] > 0) 
            fprintf(myfile, "%le %le\n", k[i], ps[i]);
    }
    fclose(myfile);
}


int main(){

        // Load power spectrum from file
    double *k = NULL, *ps = NULL;
    int k_length = import("pk.dat", &k, &ps);

    double z[3] = {100,3,0}; 

        // Initialize trg object 
    double *ps_in = trg_concatenate(k_length, k, ps);
    trg_object* trg = trg_init(z, 3, ps_in, k_length, &Omega21, &Omega22);

        // We will reuse k and ps later, so free them here:
    free(k);
    free(ps);

        // Adjust options to your liking (or leave them to their default
        // value)
    trg->options->timesteps = 40;
    trg->options->linear = false;

        // Run time evolution
    trg_evolution(trg);

        // Do something with the result:

        // Extract one power spectrum, e.g. P_11 at the last redshift
        // (the 2nd argument runs from 0 to z_length-2 (here 1), and the 3rd
        // argument corresponds to 1->P11, 2->P12, 3->P22)
    trg_extract_spectrum(trg, 1, 1, &k, &ps);
        
        // Save it to a file
    export("ps_nl.dat", k, ps, k_length);

        // Clean up
    free(k);
    free(ps);
    trg_free(trg);

    return 0;
}
