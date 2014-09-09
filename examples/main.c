/// @file main.c
/// A C program that calls the Fortran routine.
// @author Adrian Vollmer
// @date 20.07.2012


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Declaration of the Fortran functions.
extern void trgfast_init_c(
            double *ps,
            double *par_real,
            int *par_int,
            double O21(double*,double*),
            double O22(double*, double*),
            char *dir
        );

extern void trgfast_get_ps_c(
            double *z,
            double *curlyH,
            int *length,
            double *ps_out
        );

extern void trgfast_free_c();

const int len_int = 12;
const int len_real = 6;



// Define the matrix Omega for LCDM
double OmegaM0 = .25, OmegaL0 = .75, h = .73;

double Omega21(double *a, double *k){
    return -1.5 * OmegaM0/(pow(*a,3)*OmegaL0 + OmegaM0);
}

double Omega22(double *a, double *k){
    return 2. - (-2*pow(*a,3)*OmegaL0 + OmegaM0)/(2.*(pow(*a,3)*OmegaL0 + OmegaM0));
}

double curlyH(double a){
    return a*(h/3e3)*sqrt(OmegaL0 + OmegaM0/pow(a,3));
}


// Allocate an array of pointers which each point at an array of doubles.
double** Make2DDoubleArray(int arraySizeX, int arraySizeY) {
    double** theArray;
    int i;
    theArray =  malloc(arraySizeX*sizeof(double*));
    for (i = 0; i < arraySizeX; i++)
        theArray[i] =  malloc(arraySizeY*sizeof(double));   
    return theArray;
} 


// Load a 2-column table of float numbers into a 2D array with at least
// max_n rows. Empty lines and lines starting with a '#' are being ignored.
int import(const char *filename, double **table, int max_n) {
    FILE *myfile;
    char line[1024];
    int i;
    myfile = fopen(filename, "r");
    if (NULL==myfile) { printf("File not found: %s\n", filename); return -1; };
    i=0;
    while(fscanf(myfile, "%[^\n]\n", line) != EOF) {
        if (line[0]!='#') {
            sscanf(line, "%le %le\n", &table[i][0], &table[i][1]);
            i++;
            if (i==max_n) break;
        }
    }
    if (NULL!=myfile) fclose(myfile);
    return i;
}


void dump(const char *filename, double **table, int max_n) {
    FILE *myfile;
    int i;
    myfile = fopen(filename, "w");
    for (i=0; i<max_n; i++) {
        fprintf(myfile, "%le %le %le %le\n",
                table[i][0], table[i][1], table[i][2], table[i][3]);
    }
    fclose(myfile);
}


int main() {
        // Declare arrays and variables.
    double** ps = Make2DDoubleArray(1024,2);
    double par_real[len_real];
    double cH, zfin;
    int par_int[len_int];
    int n, i;

        // Read linear matter power spectrum at z=zini from file.
    n = import("pk.dat", ps, 1024);

        // Bundle all three spectra into one array.
    double psin[n][2];
    double psout[n][4];
    for (i=0; i<n; i++) {
        psin[i][0] = ps[i][0];
        psin[i][1] = ps[i][1];
    }
    
        // Fill the arrays with the parameters
    par_real[1] = 100d0;
    par_real[2] = 0d0;
    par_real[3] = -1.d0; //.999684// f(z)
    par_real[4] = .966d0; // ns
    par_real[5] = 17d0; // k_th
    par_real[6] = 1d0; // scale_A
    par_int[1] = size(ps0,2);
    par_int[2] = 1; // full nonlinear
    par_int[3] = 1; // verbosity
    par_int[4] = 50; // time evolution step size
    par_int[5] = 200; // sampling points for Kdq
    par_int[6] = 1;  // extrapolation method
    par_int[7] = 0;  // dump Kdq?
    par_int[8] = 0;  // dump A?
    par_int[9] = 1;  // dump PS?
    par_int[10] = 0; // dump diff?
    par_int[11] = 0; // dump growth?
    par_int[12] = 0; // initial conditions

        // Call Fortran functions.
        // Compute non-linear power spectrum
    trgfast_init_c((double*)psin, par_real, par_int, Omega21, Omega22, "../output/");
        // Request the power spectrum at an arbitrary redshift (must pass
        // pointers)
    zfin = 0.;
    cH = curlyH(1/(1+zfin));
    trgfast_get_ps_c(&zfin, &cH, &n, (double*)psout);

        // Do whatever with the result, e.g. write to a file.
        // First adjust dimensions of ps
    for (i=0; i<1024; i++) free(ps[i]);
    free(ps);

    ps = Make2DDoubleArray(n,4);

    for (i=0; i<n; i++){
        ps[i][0] = psout[i][0];
        ps[i][1] = psout[i][1];
        ps[i][2] = psout[i][2];
        ps[i][3] = psout[i][3];
    }
    dump("../output/c_ps_nl.dat", ps, n);

        // Here you could request the power spectrum at other z's

        // Clean up
    trgfast_free_c();

    for (i=0; i<n; i++) free(ps[i]);
    free(ps);
    return 0;
}
