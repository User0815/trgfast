// @file driver.c
// A C program that calls the Fortran routine.
// @author Adrian Vollmer
// @date 20.07.2012


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Declaration of the Fortran functions.
extern void time_evolution_c(
                double* eta, int* eta_len, 
                double* ps_in, int* ps_len, 
                double* O_eta, int* O_eta_len,
                double* O_k, int* O_k_len,
                double* OmegaBulk, int* OBulk_len, 
                double* ps_out, double* growth_out, double* opts);

const int len_options = 23;

// Define the matrix Omega for LCDM
double OmegaM0 = .25, OmegaL0 = .75, h = .73;

double Omega21(double a, double k){
    return -1.5 * OmegaM0/(pow(a,3)*OmegaL0 + OmegaM0);
}

double Omega22(double a, double k){
    return 2. - (-2*pow(a,3)*OmegaL0 + OmegaM0)/(2.*(pow(a,3)*OmegaL0 + OmegaM0));
}

// Load a 2-column table of double numbers into a 1D array. 
// Empty lines and lines starting with a '#' are being ignored.
int import(const char *filename, double *array, int array_len) {
    FILE *myfile;
    char line[1024];
    int i, file_length;
    myfile = fopen(filename, "r");
    if (NULL==myfile) { printf("File not found: %s\n", filename); return -1; };
    file_length=0;
    while(fscanf(myfile, "%[^\n]\n", line) != EOF) 
        if (line[0]!='#') file_length++;
    rewind(myfile);
    i=0;
    while(fscanf(myfile, "%[^\n]\n", line) != EOF) {
        if (line[0]!='#') {
            sscanf(line, "%le %le\n", &array[i], &array[i+file_length]);
            // sscanf(line, "%le %le\n", &array[2*i], &array[2*i+1]);
            i+=1;
            if (i+file_length == array_len) break;
        }
    }
    if (NULL!=myfile) fclose(myfile);
    return i;
}


// Writes a 1D array into a 4-column file
void export(const char *filename, double *array, int array_len) {
    FILE *myfile;
    int i;
    myfile = fopen(filename, "w");
    for (i=0; i<array_len/4; i++) {
        if (array[i+array_len/4] > 0) 
            fprintf(myfile, "%le %le %le %le\n",
                array[i], array[i+array_len/4], 
                          array[i+2*array_len/4], 
                          array[i+3*array_len/4]);
    }
    fclose(myfile);
}


int main() {
        // Declare arrays and variables.
    double *ps_in, *eta, *O_eta, *O_k, *OmegaBulk;
    int eta_len, ps_len, O_eta_len, O_k_len, OBulk_len;
    double options[len_options], zini;
    int i, j;

    zini = 100;
    eta_len = 2;
    eta = malloc(eta_len * sizeof(O_eta));
    eta[1] = log(1+zini); // final eta
    eta[0] = 0;  // initial eta

        // Read linear matter power spectrum at z=0 from file.
    ps_in = malloc(2 * 2024 * sizeof(ps_in));
    ps_len = 2*import("pk.dat", ps_in, 2*2024);

        // These arrays need to be a fixed size
    double ps_out[2*ps_len*(eta_len-1)], growth_out[1000];

        // Prepare background functions
    O_k_len = 1;
    O_k = malloc(O_k_len * sizeof(O_k));
        // This is only needed if the background depends on the scale
    // for (i=0; i<O_k_len; i++){
    //     O_k[i] = exp(log(kmin)+log(kmax-kmin)*i/(O_k_len-1));
    // }
    
    O_eta_len = 100;
    O_eta = malloc(O_eta_len * sizeof(O_eta));
    OBulk_len = 2*O_eta_len*O_k_len;
    OmegaBulk = malloc(OBulk_len * sizeof(OmegaBulk));
    for (i=0; i<O_eta_len; i++){
        O_eta[i] = i*log(1+zini)/(O_eta_len-1);
        for (j=0; j<O_k_len; j++){
            OmegaBulk[i+O_eta_len*2*j] = Omega21(exp(O_eta[i])/(1+zini), O_k[j]);
            OmegaBulk[i+O_eta_len*(2*j+1)] = Omega22(exp(O_eta[i])/(1+zini), O_k[j]);
        }
    }

        // Fill the arrays with the parameters
    options[0] = 1; // just use the default

        // Call Fortran function.
    time_evolution_c(eta, &eta_len, 
            ps_in, &ps_len,
            O_eta, &O_eta_len, 
            O_k, &O_k_len,
            OmegaBulk, &OBulk_len, 
            (double*)ps_out, 
            (double*)growth_out, 
            options);

        // Write result to a file.
    export("c_ps_nl.dat", ps_out, 2*ps_len*(eta_len-1));

        // Clean up
    free(ps_in);
    free(eta);
    free(O_eta);
    free(O_k);
    free(OmegaBulk);
    return 0;
}
