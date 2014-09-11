#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"mathlink.h"

/* symbols from Fortan are all lower case */

extern void init_a_c(double* ps, int* k_len, double* opts);
extern double a_c(double* k, double* res);
extern void clean_up_c();
extern void clean_up_ode_c();
extern void f_ode_c(double* eta, double* X, double* Xprime);
extern void init_ode_c( double* O_eta, int* O_eta_len, 
                        double* O_k, int* O_k_len, 
                        double* Omega, int* Omega_len, 
                        double* k, int* k_len,
                        double* opts);
extern void time_evolution_c(
                double* eta, int* eta_len, 
                double* ps_in, int* ps_len, 
                double* O_eta, int* O_eta_len,
                double* O_k, int* O_k_len,
                double* OmegaBulk, int* OBulk_len, 
                double* ps_out, double* growth_out, double* opts);


void time_evo(double* eta, long eta_len, 
              double* ps_in, long ps_len, 
              double* O_eta, long O_eta_len, 
              double* O_k, long O_k_len,
              double* OmegaBulk, long OBulk_len, 
              double* opts, long opts_len
             ){
    double ps_out[2*ps_len*(eta_len-1)];
    double growth_out[1000];

    time_evolution_c(eta, (int*)&eta_len, ps_in, (int*)&ps_len,
            O_eta, (int*)&O_eta_len, O_k, (int*)&O_k_len,
            OmegaBulk, (int*)&OBulk_len, (double*)ps_out, 
            (double*)growth_out, opts);
        
    MLPutFunction(stdlink, "List", 2);
    MLPutReal64List(stdlink, (double*)ps_out, 2*ps_len*(eta_len-1));
    MLPutReal64List(stdlink, (double*)growth_out, 1000);
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

void clean_up_A(){
    clean_up_c();
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

void init_A(double* ps, long ps_len, double* opts, long opts_len){
    init_a_c(ps, (int*)&ps_len, opts);
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

void getA(double k){
    double res[14];
    a_c(&k, (double*)res);
    MLPutReal64List(stdlink, (double*)res, 14);
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

void init_ode(double* O_eta, long O_eta_len,
              double* O_k, long O_k_len,
              double* OmegaBulk, long OBulk_len, 
              double* k, long k_len, 
              double* opts, long opts_len){
    init_ode_c(k, (int*)&k_len, O_eta, (int*)&O_eta_len, O_k, (int*)&O_k_len, 
            OmegaBulk, (int*)&OBulk_len, opts);
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

void f_ode(double eta, double* X, long X_len){
    double Xprime[X_len];
    f_ode_c(&eta, X, (double*)Xprime);
    MLPutReal64List(stdlink, (double*)Xprime, X_len);
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

void clean_up_ode(){
    clean_up_ode_c();
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}


int main(int argc, char* argv[]) {
    return MLMain(argc, argv);
}
