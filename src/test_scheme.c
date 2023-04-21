#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <openssl/rand.h>
#include <gmp.h>
#include "ds_benchmark.h"
#include <time.h>
#include <mpfr.h>
#include <openssl/evp.h>
#include <math.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod_poly.h"
#include "flint/fmpz_mod.h"
#include "aux.h"
#include "sampler.h"
#include "ZigguratO.h"
#include "scheme.h"
#include "permutations.h"
#include "nizkp.h"
#include "param.h"



mpz_t Q;
mpz_t BOUND;

fmpz_t fQ;
fmpz_t ALPHA;
fmpz_t OMEGA;
fmpz_t IALPHA;
fmpz_t IOMEGA;
fmpz_t IDOS;

fmpz_mod_ctx_t ctx;
fmpz_mod_poly_t * modulus;
fmpz_mod_poly_t * modulus_inv;

mpf_t SIGMA;
int LAMB;
int SEEDOL;
int SEEDM;
size_t N;
size_t KAPPA;
size_t K;
int D;
int deltaOL;
int deltaM;

mpf_t YBar[NUMBEROFRECTANGLES + 1];
mpz_t XBar[NUMBEROFRECTANGLES + 1];

int ROW_NUMBER;

int BQ;
int bQ;
unsigned char auxSampler;

const EVP_MD *MD_COM_AUX;
int md_len_com;


void test_scheme_verifier(int lambdas[], int degrees[], int Ds[], mpz_t Q){

    fmpz_mod_poly_t ** a;
    fmpz_init_array_fft(&a,K);
    fmpz_mod_poly_t ** b;
    fmpz_init_array_fft(&b,K);
    fmpz_mod_poly_t ** c;
    fmpz_init_array_fft(&c,K);

    mpz_t * e = (mpz_t *) malloc( sizeof(mpz_t) * N * K);

    mpz_p_init_array(e,N*K);

    fmpz_mod_poly_t * m;
    fmpz_init_fft(&m);
    fmpz_sampler_zq_poly_FFT(m);

    fmpz_mod_poly_t * r;
    fmpz_init_fft(&r);

    char buffer[FILENAME_MAX];
    gmp_snprintf(buffer, sizeof(buffer), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER],FILEVERIFICATION, Q);

    key_gen(a,b,NULL,NULL);
    commit(m,a,b,c,e,r);
    int res;
    TIME_OPERATION_ITERATIONS(res = verify(c, m, r, e, a, b), "verifier", ITERATIONS_SCHEME,buffer)

    test_message("opening verification",res);

    fmpz_clear_fft(m);
    fmpz_clear_fft(r);

    fmpz_clear_array_fft(a,K);
    fmpz_clear_array_fft(b,K);
    fmpz_clear_array_fft(c,K);

    mpz_p_free_array(e, N*K);
    free(e);
}


void test_scheme_keygen(int lambdas[], int degrees[], int Ds[], mpz_t Q){

    fmpz_mod_poly_t ** a;
    fmpz_init_array_fft(&a,K);
    fmpz_mod_poly_t ** b;
    fmpz_init_array_fft(&b,K);
    fmpz_mod_poly_t * a_inv;
    fmpz_init_fft(&a_inv);
    int * pos = (int *) malloc( sizeof(int) * D );

    char buffer[FILENAME_MAX];
    gmp_snprintf(buffer, sizeof(buffer), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILEKEYGEN, Q);

    TIME_OPERATION_ITERATIONS(key_gen(a,b,a_inv,pos), "keygen", ITERATIONS_SCHEME,buffer)

    test_message("keygen without errors",1);

    fmpz_clear_array_fft(a,K);
    fmpz_clear_array_fft(b,K);

    fmpz_clear_fft(a_inv);
    free(pos);

    printf("\n");
}


void test_scheme_commitment(int lambdas[], int degrees[], int Ds[], mpz_t Q){

    fmpz_mod_poly_t ** a;
    fmpz_init_array_fft(&a,K);
    fmpz_mod_poly_t ** b;
    fmpz_init_array_fft(&b,K);
    fmpz_mod_poly_t ** c;
    fmpz_init_array_fft(&c,K);

    mpz_t * e = (mpz_t *) malloc( sizeof(mpz_t) * N * K);

    mpz_p_init_array(e,N*K);

    fmpz_mod_poly_t * m;
    fmpz_init_fft(&m);
    fmpz_sampler_zq_poly_FFT(m);

    fmpz_mod_poly_t * r;
    fmpz_init_fft(&r);

    char buffer[FILENAME_MAX];
    gmp_snprintf(buffer, sizeof(buffer), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILECOMMITMENT, Q);

    key_gen(a,b,NULL,NULL);
    TIME_OPERATION_ITERATIONS(commit(m,a,b,c,e,r), "commitment", ITERATIONS_SCHEME,buffer)

    test_message("commitment without errors",1);

    fmpz_clear_fft(m);
    fmpz_clear_fft(r);

    fmpz_clear_array_fft(a,K);
    fmpz_clear_array_fft(b,K);
    fmpz_clear_array_fft(c,K);

    mpz_p_free_array(e, N*K);
    free(e);

    printf("\n");
}


void test_knowledge_valid_opening(int lambdas[], int degrees[], int Ds[], mpz_t Q){

        fmpz_mod_poly_t ** a;
        fmpz_init_array_fft(&a,K);
        fmpz_mod_poly_t ** b;
        fmpz_init_array_fft(&b,K);

        fmpz_mod_poly_t ** c;
        fmpz_init_array_fft(&c,K);

        mpz_t * e = (mpz_t *) malloc( sizeof(mpz_t) * N * K);

        mpz_p_init_array(e,N*K);

        fmpz_mod_poly_t * m;
        fmpz_init_fft(&m);
        fmpz_sampler_zq_poly_FFT(m);

        fmpz_mod_poly_t * r;
        fmpz_init_fft(&r);

        fmpz_mod_poly_t * a_inv;
        fmpz_init_fft(&a_inv);
        int * pos = (int *) malloc( sizeof(int) * D );

        key_gen(a,b,a_inv,pos);
        commit(m,a,b,c,e,r);

        unsigned char ** seed = (unsigned char **) malloc(deltaOL * sizeof(unsigned char *));
        for (int rep = 0; rep < deltaOL; rep++){
            seed[rep] = malloc(LAMB*sizeof(unsigned char));
        }

        mpz_t *** g = (mpz_t ***) malloc(deltaOL *  sizeof(mpz_t **));
        for (int rep = 0; rep < deltaOL; rep++){
            g[rep] = (mpz_t **) malloc( (KAPPA + 1) * sizeof(mpz_t *) );
            for(int i = 0; i < KAPPA + 1; i++){
                g[rep][i] = (mpz_t *)malloc(2*N*K * sizeof(mpz_t));
                mpz_p_init_array(g[rep][i],2*N*K);
            }
        }


        fmpz_mod_poly_t *** y = (fmpz_mod_poly_t ***) malloc(deltaOL *  sizeof(fmpz_mod_poly_t **));
        for (int rep = 0; rep < deltaOL; rep++){
            fmpz_init_array_fft(&y[rep],K);
        }


        unsigned char ** md_value_com1 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        unsigned char ** md_value_com2 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        unsigned char ** d1 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        unsigned char ** d2 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        for (int rep = 0; rep < deltaOL; rep++){
            md_value_com1[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            md_value_com2[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            d1[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
            d2[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
        }

        fmpz_mod_poly_t ** s;
        fmpz_init_array_fft(&s,deltaOL);

        unsigned long *** e_tilde = (unsigned long ***) malloc(deltaOL*sizeof(unsigned long **));
        for (int rep = 0; rep < deltaOL; rep++){
            e_tilde[rep] = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
            for(int i = 0; i < KAPPA + 1; i++){
                e_tilde[rep][i] = (unsigned long *) malloc(N*K/4);
            }
        }

        char buffer_prover[FILENAME_MAX];
        gmp_snprintf(buffer_prover, sizeof(buffer_prover), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILEPROVEROPENING, Q);
        char buffer_verifier[FILENAME_MAX];
        gmp_snprintf(buffer_verifier, sizeof(buffer_verifier), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILEVERIFIEROPENING, Q);


        TIME_OPERATION_ITERATIONS(prover_opening(a, b, c, m, r, e, deltaOL, md_value_com1, md_value_com2, g, seed, y, s, d1, e_tilde, d2), "prover opening", ITERATIONS_SCHEME, buffer_prover);

        int res;
        TIME_OPERATION_ITERATIONS(res = verifier_opening(a, b, c, deltaOL, md_value_com1, md_value_com2, g, seed, y, s, d1, e_tilde, d2, pos), "verifier opening", ITERATIONS_SCHEME, buffer_verifier);


        test_message("nizkp opening verifier",res);

        for (int rep = 0; rep < deltaOL; rep++){
            free(seed[rep]);
            free(md_value_com1[rep]);
            free(md_value_com2[rep]);
            free(d1[rep]);
            free(d2[rep]);
        }

        free(seed);
        free(md_value_com1);
        free(md_value_com2);
        free(d1);
        free(d2);

        fmpz_clear_array_fft(s,deltaOL);
        for (int rep = 0; rep < deltaOL; rep++){
            fmpz_clear_array_fft(y[rep],K);
        }
        free(y);

        for (int rep = 0; rep < deltaOL; rep++){
            mpz_p_free_matrix(g[rep], KAPPA + 1, 2*N*K);
        }
        free(g);

        for (int rep = 0; rep < deltaOL; rep++){
            for(int i = 0; i < KAPPA + 1; i++){
                free(e_tilde[rep][i]);
            }
            free(e_tilde[rep]);
        }
        free(e_tilde);

        fmpz_clear_fft(m);
        fmpz_clear_fft(r);

        fmpz_clear_array_fft(a,K);
        fmpz_clear_array_fft(b,K);
        fmpz_clear_array_fft(c,K);

        fmpz_clear_fft(a_inv);
        free(pos);

        mpz_p_free_array(e, N*K);
        free(e);
}

void test_knowledge_linear_relation(int lambdas[], int degrees[], int Ds[], mpz_t Q){

        fmpz_mod_poly_t ** a, ** b;
        fmpz_init_array_fft(&a,K);
        fmpz_init_array_fft(&b,K);

        fmpz_mod_poly_t ** c1,  ** c2,  ** c3;
        fmpz_init_array_fft(&c1,K);
        fmpz_init_array_fft(&c2,K);
        fmpz_init_array_fft(&c3,K);

        mpz_t * e1 = (mpz_t *) malloc( sizeof(mpz_t) * N * K);
        mpz_t * e2 = (mpz_t *) malloc( sizeof(mpz_t) * N * K);
        mpz_t * e3 = (mpz_t *) malloc( sizeof(mpz_t) * N * K);

        mpz_p_init_array(e1, N*K);
        mpz_p_init_array(e2, N*K);
        mpz_p_init_array(e3, N*K);

        fmpz_mod_poly_t * m1;
        fmpz_init_fft(&m1);
        fmpz_sampler_zq_poly_FFT(m1);

        fmpz_mod_poly_t * m2;
        fmpz_init_fft(&m2);
        fmpz_sampler_zq_poly_FFT(m2);

        fmpz_mod_poly_t * m3;
        fmpz_init_fft(&m3);

        fmpz_mod_poly_t * lambda1;
        fmpz_init_fft(&lambda1);
        fmpz_sampler_zq_poly_FFT(lambda1);

        fmpz_mod_poly_t * lambda2;
        fmpz_init_fft(&lambda2);
        fmpz_sampler_zq_poly_FFT(lambda2);

        fmpz_mod_poly_t * lambda1_aux;
        fmpz_init_fft(&lambda1_aux);

        fmpz_mod_poly_t * lambda2_aux;
        fmpz_init_fft(&lambda2_aux);

        fmpz_mod_poly_mult_FFT(lambda1_aux, m1, lambda1);
        fmpz_mod_poly_mult_FFT(lambda2_aux, m2, lambda2);

        fmpz_mod_poly_add_FFT(m3, lambda1_aux, lambda2_aux);

        fmpz_mod_poly_t * r1, * r2, * r3;
        fmpz_init_fft(&r1);
        fmpz_init_fft(&r2);
        fmpz_init_fft(&r3);

        fmpz_mod_poly_t * a_inv;
        fmpz_init_fft(&a_inv);
        int * pos = (int *) malloc( sizeof(int) * D );

        key_gen(a,b,a_inv,pos);

        commit(m1,a,b,c1,e1,r1);
        commit(m2,a,b,c2,e2,r2);
        commit(m3,a,b,c3,e3,r3);

        unsigned char ** seed1 = (unsigned char **) malloc(deltaOL * sizeof(unsigned char *));
        unsigned char ** seed2 = (unsigned char **) malloc(deltaOL * sizeof(unsigned char *));
        unsigned char ** seed3 = (unsigned char **) malloc(deltaOL * sizeof(unsigned char *));
        for (int rep = 0; rep < deltaOL; rep++){
            seed1[rep] = malloc(LAMB*sizeof(unsigned char));
            seed2[rep] = malloc(LAMB*sizeof(unsigned char));
            seed3[rep] = malloc(LAMB*sizeof(unsigned char));
        }

        mpz_t *** g1 = (mpz_t ***) malloc(deltaOL *  sizeof(mpz_t **));
        mpz_t *** g2 = (mpz_t ***) malloc(deltaOL *  sizeof(mpz_t **));
        mpz_t *** g3 = (mpz_t ***) malloc(deltaOL *  sizeof(mpz_t **));
        for (int rep = 0; rep < deltaOL; rep++){
            g1[rep] = (mpz_t **) malloc( (KAPPA + 1) * sizeof(mpz_t *) );
            g2[rep] = (mpz_t **) malloc( (KAPPA + 1) * sizeof(mpz_t *) );
            g3[rep] = (mpz_t **) malloc( (KAPPA + 1) * sizeof(mpz_t *) );
            for(int i = 0; i < KAPPA + 1; i++){
                g1[rep][i] = (mpz_t *)malloc(2*N*K * sizeof(mpz_t));
                g2[rep][i] = (mpz_t *)malloc(2*N*K * sizeof(mpz_t));
                g3[rep][i] = (mpz_t *)malloc(2*N*K * sizeof(mpz_t));
                mpz_p_init_array(g1[rep][i],2*N*K);
                mpz_p_init_array(g2[rep][i],2*N*K);
                mpz_p_init_array(g3[rep][i],2*N*K);
            }
        }


        fmpz_mod_poly_t *** y1 = (fmpz_mod_poly_t ***) malloc(deltaOL *  sizeof(fmpz_mod_poly_t **));
        fmpz_mod_poly_t *** y2 = (fmpz_mod_poly_t ***) malloc(deltaOL *  sizeof(fmpz_mod_poly_t **));
        fmpz_mod_poly_t *** y3 = (fmpz_mod_poly_t ***) malloc(deltaOL *  sizeof(fmpz_mod_poly_t **));
        for (int rep = 0; rep < deltaOL; rep++){
            fmpz_init_array_fft(&y1[rep],K);
            fmpz_init_array_fft(&y2[rep],K);
            fmpz_init_array_fft(&y3[rep],K);
        }


        unsigned char ** md_value_com1 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        unsigned char ** md_value_com2 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        unsigned char ** d1 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        unsigned char ** d2 = (unsigned char **) malloc(deltaOL*sizeof(unsigned char*));
        for (int rep = 0; rep < deltaOL; rep++){
            md_value_com1[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            md_value_com2[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            d1[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
            d2[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
        }

        fmpz_mod_poly_t ** s1, ** s2, ** s3;
        fmpz_init_array_fft(&s1,deltaOL);
        fmpz_init_array_fft(&s2,deltaOL);
        fmpz_init_array_fft(&s3,deltaOL);

        unsigned long *** e_tilde1 = (unsigned long ***) malloc(deltaOL*sizeof(unsigned long **));
        unsigned long *** e_tilde2 = (unsigned long ***) malloc(deltaOL*sizeof(unsigned long **));
        unsigned long *** e_tilde3 = (unsigned long ***) malloc(deltaOL*sizeof(unsigned long **));
        for (int rep = 0; rep < deltaOL; rep++){
            e_tilde1[rep] = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
            e_tilde2[rep] = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
            e_tilde3[rep] = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
            for(int i = 0; i < KAPPA + 1; i++){
                e_tilde1[rep][i] = (unsigned long *) malloc(N*K/4);
                e_tilde2[rep][i] = (unsigned long *) malloc(N*K/4);
                e_tilde3[rep][i] = (unsigned long *) malloc(N*K/4);
            }
        }

        char buffer_prover[FILENAME_MAX];
        gmp_snprintf(buffer_prover, sizeof(buffer_prover), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILEPROVERLINEAR, Q);
        char buffer_verifier[FILENAME_MAX];
        gmp_snprintf(buffer_verifier, sizeof(buffer_verifier), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILEVERIFIERLINEAR, Q);



        TIME_OPERATION_ITERATIONS(prover_linear(a, b, c1, c2, c3, m1, m2, m3, r1, r2, r3, e1, e2, e3, lambda1, lambda2, deltaOL, md_value_com1, md_value_com2, g1, g2, g3, seed1, seed2, seed3, y1, y2, y3, s1, s2, s3, d1, e_tilde1, e_tilde2, e_tilde3, d2), "prover linear", ITERATIONS_SCHEME, buffer_prover);


        int res;
        TIME_OPERATION_ITERATIONS(res = verifier_linear(a, b, c1, c2, c3, lambda1, lambda2, deltaOL, md_value_com1, md_value_com2, g1, g2, g3, seed1, seed2, seed3, y1, y2, y3, s1, s2, s3, d1, e_tilde1, e_tilde2, e_tilde3, d2, pos), "verifier linear", ITERATIONS_SCHEME, buffer_verifier);


        test_message("nizkp linear verifier",res);

        fmpz_clear_fft(lambda1);
        fmpz_clear_fft(lambda2);
        fmpz_clear_fft(lambda1_aux);
        fmpz_clear_fft(lambda2_aux);

        for (int rep = 0; rep < deltaOL; rep++){
            free(seed1[rep]);
            free(seed2[rep]);
            free(seed3[rep]);
            free(md_value_com1[rep]);
            free(md_value_com2[rep]);
            free(d1[rep]);
            free(d2[rep]);
        }

        free(seed1);
        free(seed2);
        free(seed3);
        free(md_value_com1);
        free(md_value_com2);
        free(d1);
        free(d2);

        fmpz_clear_array_fft(s1,deltaOL);
        fmpz_clear_array_fft(s2,deltaOL);
        fmpz_clear_array_fft(s3,deltaOL);
        for (int rep = 0; rep < deltaOL; rep++){
            fmpz_clear_array_fft(y1[rep],K);
            fmpz_clear_array_fft(y2[rep],K);
            fmpz_clear_array_fft(y3[rep],K);
        }
        free(y1);
        free(y2);
        free(y3);

        for (int rep = 0; rep < deltaOL; rep++){
            mpz_p_free_matrix(g1[rep], KAPPA + 1, 2*N*K);
            mpz_p_free_matrix(g2[rep], KAPPA + 1, 2*N*K);
            mpz_p_free_matrix(g3[rep], KAPPA + 1, 2*N*K);
        }
        free(g1);
        free(g2);
        free(g3);



        for (int rep = 0; rep < deltaOL; rep++){
            for(int i = 0; i < KAPPA + 1; i++){
                free(e_tilde1[rep][i]);
                free(e_tilde2[rep][i]);
                free(e_tilde3[rep][i]);
            }
            free(e_tilde1[rep]);
            free(e_tilde2[rep]);
            free(e_tilde3[rep]);
        }

        free(e_tilde1);
        free(e_tilde2);
        free(e_tilde3);

        fmpz_clear_fft(m1);
        fmpz_clear_fft(m2);
        fmpz_clear_fft(m3);
        fmpz_clear_fft(r1);
        fmpz_clear_fft(r2);
        fmpz_clear_fft(r3);

        fmpz_clear_array_fft(a,K);
        fmpz_clear_array_fft(b,K);
        fmpz_clear_array_fft(c1,K);
        fmpz_clear_array_fft(c2,K);
        fmpz_clear_array_fft(c3,K);

        fmpz_clear_fft(a_inv);
        free(pos);

        mpz_p_free_array(e1, N*K);
        mpz_p_free_array(e2, N*K);
        mpz_p_free_array(e3, N*K);
        free(e1);
        free(e2);
        free(e3);
}


void test_knowledge_multiplicative_relation(int lambdas[], int degrees[], int Ds[], mpz_t Q){

        fmpz_mod_poly_t ** a;
        fmpz_init_array_fft(&a,K);
        fmpz_mod_poly_t ** b;
        fmpz_init_array_fft(&b,K);

        fmpz_mod_poly_t ** c1;
        fmpz_init_array_fft(&c1,K);
        fmpz_mod_poly_t ** c2;
        fmpz_init_array_fft(&c2,K);
        fmpz_mod_poly_t ** c3;
        fmpz_init_array_fft(&c3,K);

        mpz_t * e1 = (mpz_t *) malloc( sizeof(mpz_t) * N * K);
        mpz_t * e2 = (mpz_t *) malloc( sizeof(mpz_t) * N * K);
        mpz_t * e3 = (mpz_t *) malloc( sizeof(mpz_t) * N * K);

        mpz_p_init_array(e1, N*K);
        mpz_p_init_array(e2, N*K);
        mpz_p_init_array(e3, N*K);

        fmpz_mod_poly_t * m1;
        fmpz_init_fft(&m1);
        fmpz_sampler_zq_poly_FFT(m1);

        fmpz_mod_poly_t * m2;
        fmpz_init_fft(&m2);
        fmpz_sampler_zq_poly_FFT(m2);

        fmpz_mod_poly_t * m3;
        fmpz_init_fft(&m3);

        fmpz_mod_poly_mult_FFT(m3,m1,m2);

        fmpz_mod_poly_t * r1, * r2, * r3;
        fmpz_init_fft(&r1);
        fmpz_init_fft(&r2);
        fmpz_init_fft(&r3);

        fmpz_mod_poly_t * a_inv;
        fmpz_init_fft(&a_inv);
        int * pos = (int *) malloc( sizeof(int) * D );

        key_gen(a,b,a_inv,pos);

        commit(m1,a,b,c1,e1,r1);
        commit(m2,a,b,c2,e2,r2);
        commit(m3,a,b,c3,e3,r3);

        unsigned char ** seed1 = (unsigned char **) malloc(deltaM * sizeof(unsigned char *));
        unsigned char ** seed2 = (unsigned char **) malloc(deltaM * sizeof(unsigned char *));
        unsigned char ** seed3 = (unsigned char **) malloc(deltaM * sizeof(unsigned char *));
        for (int rep = 0; rep < deltaM; rep++){
            seed1[rep] = malloc(LAMB*sizeof(unsigned char));
            seed2[rep] = malloc(LAMB*sizeof(unsigned char));
            seed3[rep] = malloc(LAMB*sizeof(unsigned char));
        }

        mpz_t *** g1 = (mpz_t ***) malloc(deltaM *  sizeof(mpz_t **));
        mpz_t *** g2 = (mpz_t ***) malloc(deltaM *  sizeof(mpz_t **));
        mpz_t *** g3 = (mpz_t ***) malloc(deltaM *  sizeof(mpz_t **));
        for (int rep = 0; rep < deltaM; rep++){
            g1[rep] = (mpz_t **) malloc( (KAPPA + 1) * sizeof(mpz_t *) );
            g2[rep] = (mpz_t **) malloc( (KAPPA + 1) * sizeof(mpz_t *) );
            g3[rep] = (mpz_t **) malloc( (KAPPA + 1) * sizeof(mpz_t *) );
            for(int i = 0; i < KAPPA + 1; i++){
                g1[rep][i] = (mpz_t *)malloc(2*N*K * sizeof(mpz_t));
                g2[rep][i] = (mpz_t *)malloc(2*N*K * sizeof(mpz_t));
                g3[rep][i] = (mpz_t *)malloc(2*N*K * sizeof(mpz_t));
                mpz_p_init_array(g1[rep][i],2*N*K);
                mpz_p_init_array(g2[rep][i],2*N*K);
                mpz_p_init_array(g3[rep][i],2*N*K);
            }
        }

        fmpz_mod_poly_t *** y1 = (fmpz_mod_poly_t ***) malloc(deltaM *  sizeof(fmpz_mod_poly_t **));
        fmpz_mod_poly_t *** y2 = (fmpz_mod_poly_t ***) malloc(deltaM *  sizeof(fmpz_mod_poly_t **));
        fmpz_mod_poly_t *** y3 = (fmpz_mod_poly_t ***) malloc(deltaM *  sizeof(fmpz_mod_poly_t **));
        for (int rep = 0; rep < deltaM; rep++){
            fmpz_init_array_fft(&y1[rep],K);
            fmpz_init_array_fft(&y2[rep],K);
            fmpz_init_array_fft(&y3[rep],K);
        }

        unsigned char ** md_value_com1 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** md_value_com2 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** md_value_com3 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** md_value_com4 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** md_value_com5 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** d1 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** d2 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** d3 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** d4 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        unsigned char ** d5 = (unsigned char **) malloc(deltaM*sizeof(unsigned char*));
        for (int rep = 0; rep < deltaM; rep++){
            md_value_com1[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            md_value_com2[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            md_value_com3[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            md_value_com4[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            md_value_com5[rep] = (unsigned char *) malloc(md_len_com*sizeof(unsigned char));
            d1[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
            d2[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
            d3[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
            d4[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
            d5[rep] = (unsigned char *) malloc(LAMB*sizeof(unsigned char));
        }

        fmpz_mod_poly_t ** s1, ** s2, ** s3;
        fmpz_init_array_fft(&s1,deltaM);
        fmpz_init_array_fft(&s2,deltaM);
        fmpz_init_array_fft(&s3,deltaM);

        fmpz_mod_poly_t ** t_x, ** t_plus;
        fmpz_init_array_fft(&t_x,deltaM);
        fmpz_init_array_fft(&t_plus,deltaM);

        unsigned long ***  e_tilde1 = (unsigned long ***) malloc(deltaM*sizeof(unsigned long **));
        unsigned long ***  e_tilde2 = (unsigned long ***) malloc(deltaM*sizeof(unsigned long **));
        unsigned long ***  e_tilde3 = (unsigned long ***) malloc(deltaM*sizeof(unsigned long **));
        for (int rep = 0; rep < deltaM; rep++){
             e_tilde1[rep] = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
             e_tilde2[rep] = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
             e_tilde3[rep] = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
            for(int i = 0; i < KAPPA + 1; i++){
                 e_tilde1[rep][i] = (unsigned long *) malloc(N*K/4);
                 e_tilde2[rep][i] = (unsigned long *) malloc(N*K/4);
                 e_tilde3[rep][i] = (unsigned long *) malloc(N*K/4);
            }
        }

        fmpz_mod_poly_t ** mu3, ** mu_x, ** mu_plus;
        fmpz_init_array_fft(&mu3,deltaM);
        fmpz_init_array_fft(&mu_x,deltaM);
        fmpz_init_array_fft(&mu_plus,deltaM);

        char buffer_prover[FILENAME_MAX];
        gmp_snprintf(buffer_prover, sizeof(buffer_prover), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILEPROVERMULTIPLICATIVE, Q);
        char buffer_verifier[FILENAME_MAX];
        gmp_snprintf(buffer_verifier, sizeof(buffer_verifier), "benchmarks/%d-%d-%d/%s_%Zd.csv",  lambdas[ROW_NUMBER], degrees[ROW_NUMBER], Ds[ROW_NUMBER], FILEVERIFIERMULTIPLICATIVE, Q);


        TIME_OPERATION_ITERATIONS(prover_multiplicative(a, b, c1, c2, c3, m1, m2, m3, r1, r2, r3, e1, e2, e3, deltaM, md_value_com1, md_value_com2, md_value_com3, md_value_com4, g1, g2, g3, md_value_com5, seed1, seed2, seed3, y1, y2, y3, s1, s2, s3, t_x, t_plus, d1, d4, d5, e_tilde1, e_tilde2, e_tilde3, mu3, mu_x, mu_plus, d2, d3), "prover multiplicative", ITERATIONS_SCHEME, buffer_prover);

        int res;
        TIME_OPERATION_ITERATIONS(res = verifier_multiplicative(a, b, c1, c2, c3, deltaM, md_value_com1, md_value_com2, md_value_com3, md_value_com4, md_value_com5, g1, g2, g3, seed1, seed2, seed3, y1, y2, y3, s1, s2, s3, t_x, t_plus, d1, d4, d5, e_tilde1, e_tilde2, e_tilde3, mu3, mu_x, mu_plus, d2, d3, a_inv, pos), "verifier multiplicative", ITERATIONS_SCHEME, buffer_verifier);


        test_message("nizkp multiplicative verifier",res);

        for (int rep = 0; rep < deltaM; rep++){
            free(seed1[rep]);
            free(seed2[rep]);
            free(seed3[rep]);
            free(md_value_com1[rep]);
            free(md_value_com2[rep]);
            free(md_value_com3[rep]);
            free(md_value_com4[rep]);
            free(md_value_com5[rep]);
            free(d1[rep]);
            free(d2[rep]);
            free(d3[rep]);
            free(d4[rep]);
            free(d5[rep]);
        }

        free(seed1);
        free(seed2);
        free(seed3);
        free(md_value_com1);
        free(md_value_com2);
        free(md_value_com3);
        free(md_value_com4);
        free(md_value_com5);
        free(d1);
        free(d2);
        free(d3);
        free(d4);
        free(d5);

        fmpz_clear_array_fft(s1,deltaM);
        fmpz_clear_array_fft(s2,deltaM);
        fmpz_clear_array_fft(s3,deltaM);

        fmpz_clear_array_fft(mu3,deltaM);
        fmpz_clear_array_fft(t_x,deltaM);
        fmpz_clear_array_fft(t_plus,deltaM);
        fmpz_clear_array_fft(mu_x,deltaM);
        fmpz_clear_array_fft(mu_plus,deltaM);

        for (int rep = 0; rep < deltaM; rep++){
            fmpz_clear_array_fft(y1[rep],K);
            fmpz_clear_array_fft(y2[rep],K);
            fmpz_clear_array_fft(y3[rep],K);
        }
        free(y1);
        free(y2);
        free(y3);

        for (int rep = 0; rep < deltaM; rep++){
            mpz_p_free_matrix(g1[rep], KAPPA + 1, 2*N*K);
            mpz_p_free_matrix(g2[rep], KAPPA + 1, 2*N*K);
            mpz_p_free_matrix(g3[rep], KAPPA + 1, 2*N*K);
        }
        free(g1);
        free(g2);
        free(g3);

        for (int rep = 0; rep < deltaM; rep++){
            for(int i = 0; i < KAPPA + 1; i++){
                free(e_tilde1[rep][i]);
                free(e_tilde2[rep][i]);
                free(e_tilde3[rep][i]);
            }
            free(e_tilde1[rep]);
            free(e_tilde2[rep]);
            free(e_tilde3[rep]);
        }

        free(e_tilde1);
        free(e_tilde2);
        free(e_tilde3);

        fmpz_clear_fft(m1);
        fmpz_clear_fft(m2);
        fmpz_clear_fft(m3);
        fmpz_clear_fft(r1);
        fmpz_clear_fft(r2);
        fmpz_clear_fft(r3);

        fmpz_clear_array_fft(a,K);
        fmpz_clear_array_fft(b,K);
        fmpz_clear_array_fft(c1,K);
        fmpz_clear_array_fft(c2,K);
        fmpz_clear_array_fft(c3,K);

        fmpz_clear_fft(a_inv);
        free(pos);

        mpz_p_free_array(e1, N*K);
        mpz_p_free_array(e2, N*K);
        mpz_p_free_array(e3, N*K);
        free(e1);
        free(e2);
        free(e3);
}


void read_csv(char * file, int lambdas[], int degrees[], mpz_t moduli[], mpf_t sigmas[], mpz_t Bs[], int Ks[], int Ds[], int deltaOLs[], int deltaMs[], int number_params){
    FILE * fp;
    fp = fopen (file,"r");
    char comma;
    int length_string = 0;
    char * s;

    // skip header line and count the number of parameters in each row
    int count_params = 1;
    do {
        comma = fgetc(fp);
        if (comma == ',') {count_params++;};
    } while (comma != '\n');

    for(int j = 0; j < number_params; j++){
        for(int i = j*count_params; i < count_params + j*count_params; i++){
            s = (char *)calloc(100,sizeof(char));
            comma = fgetc(fp);

            while(isdigit(comma) || comma == '.'){
                s[length_string] = comma;
                comma = fgetc(fp);
                length_string++;
            }


            s[length_string]='\0';

            if(strlen(s) != 0){
                switch (i%count_params){
                    case 0:
                        lambdas[j] = atoi(s);
                        break;
                    case 1:
                        degrees[j] = atoi(s);
                        break;
                    case 2:
                        mpz_set_str(moduli[j],s,10);
                        break;
                    case 3:
                        Ds[j] = atoi(s);
                        break;
                    case 4:
                        Ks[j] = atoi(s);
                        break;
                    case 5:
                        mpf_set_str(sigmas[j],s,10);
                        break;
                    case 6:
                        mpz_set_str(Bs[j],s,10);
                        break;
                    case 7:
                        deltaOLs[j] = atoi(s);
                        break;
                    case 8:
                        deltaMs[j] = atoi(s);
                        break;
                    default:
                        break;
                }
            } else {
                i--;
            }
            length_string = 0;
            free(s);
        }
    }
    fclose (fp);
}


int main(int argc, char **argv) {

    mpf_set_default_prec(512);
    printf("\n");

    if(argc < 3){
        printf("Incorrect input arguments \n");
        return 1;
    }

    printf("%s\n", argv[1]);
    printf("%s\n", argv[2]);


    ROW_NUMBER = atoi(argv[2]);

    //ASSIGNATION OF GLOBAL VARIABLES

    mpz_init(Q);
    mpf_init(SIGMA);

    fmpz_t fDOS;
    fmpz_init(fDOS);
    fmpz_set_ui(fDOS,2);

    int number_params = 0;
    FILE * fp;
    fp = fopen (FILEPARAMETERS,"r");
    char auxchar;
    do {
        auxchar = fgetc(fp);
        if (auxchar == '\n') {number_params++;};
    } while (auxchar != EOF);

    number_params--;

    fclose(fp);

    int lambdas[number_params];
    int degrees[number_params];
    mpz_t Bs[number_params];
    int Ks[number_params];
    int Ds[number_params];
    int deltaOLs[number_params];
    int deltaMs[number_params];
    mpz_t moduli[number_params];
    mpf_t sigmas[number_params];
    mpz_p_init_array(moduli,number_params);
    mpz_p_init_array(Bs,number_params);
    mpf_p_init_array(sigmas,number_params);



    read_csv(FILEPARAMETERS,lambdas,degrees,moduli,sigmas,Bs,Ks,Ds,deltaOLs,deltaMs,number_params);

    mpf_t X[NUMBEROFRECTANGLES+1];
    mpf_p_init_array(X,NUMBEROFRECTANGLES+1);
    mpf_p_init_array(YBar,NUMBEROFRECTANGLES+1);
    mpz_p_init_array(XBar,NUMBEROFRECTANGLES+1);


    printf("---- ROW = %d ------\n", ROW_NUMBER);


    mpf_set(SIGMA,sigmas[ROW_NUMBER]);
    mpz_set(Q,moduli[ROW_NUMBER]);
    samplerAuxiliaryParameters(&bQ,&BQ,&auxSampler);
    LAMB = ceil((double)lambdas[ROW_NUMBER]/8);
    N = degrees[ROW_NUMBER];
    mpz_set(BOUND,Bs[ROW_NUMBER]);
    KAPPA = mpz_sizeinbase(BOUND,2)-1;
    K = Ks[ROW_NUMBER];
    D = Ds[ROW_NUMBER];
    deltaOL = deltaOLs[ROW_NUMBER];
    deltaM = deltaMs[ROW_NUMBER];
    SEEDOL = ceil(((double)lambdas[ROW_NUMBER]+deltaOL*bQ)/8);
    SEEDM = ceil(((double)lambdas[ROW_NUMBER]+deltaM*bQ)/8);



    fmpz_set_mpz(fQ,Q);
    root_of_unity(ALPHA,fQ,2*D);
    root_of_unity(OMEGA,fQ,D);


    fmpz_invmod(IALPHA,ALPHA,fQ);
    fmpz_invmod(IOMEGA,OMEGA,fQ);
    fmpz_invmod(IDOS,fDOS,fQ);

    fmpz_mod_ctx_init(ctx, fQ);
    modulus = (fmpz_mod_poly_t *) malloc(D*sizeof(fmpz_mod_poly_t));
    fmpz_t aux_modulus;
    fmpz_init(aux_modulus);
    modulus_inv = (fmpz_mod_poly_t *) malloc(D*sizeof(fmpz_mod_poly_t));
    fmpz_mod_poly_t aux_modulus_inv;
    fmpz_mod_poly_init(aux_modulus_inv,ctx);
    fmpz_mod_poly_set_coeff_ui(aux_modulus_inv,N/D,1,ctx);
    for (int i = 0; i < D; i++) {
        fmpz_mod_poly_init(modulus[i],ctx);
        fmpz_powm_ui(aux_modulus,OMEGA,reverseBits(i,(int)floor(log2(D))),fQ);
        fmpz_mul(aux_modulus,aux_modulus,ALPHA);
        fmpz_neg(aux_modulus,aux_modulus);
        fmpz_mod(aux_modulus,aux_modulus,fQ);

        fmpz_mod_poly_set_coeff_ui(modulus[i], N/D, 1, ctx);
        fmpz_mod_poly_set_coeff_fmpz(modulus[i], 0, aux_modulus,ctx);

        fmpz_mod_poly_init(modulus_inv[i],ctx);
        fmpz_mod_poly_reverse(modulus_inv[i],modulus[i],N/D+1,ctx);
        int e = fmpz_mod_poly_invmod(modulus_inv[i],modulus_inv[i],aux_modulus_inv,ctx);
        if (e == 0) {
            red();
            printf("ERROR: computing finv for fmpz_mod_poly_mulmod_preinv\n");
            reset();
            exit(0);
        }
    }
    fmpz_clear(aux_modulus);
    fmpz_mod_poly_clear(aux_modulus_inv,ctx);


    printf("N = %zu\n", N);
    gmp_printf ("Q = %Zd\n", Q);
    printf("K = %zu\n", K);
    gmp_printf ("BOUND = %Zd\n", BOUND);
    printf("D = %d\n", D);
    printf("deltaOL = %d\n", deltaOL);
    printf("deltaM = %d\n", deltaM);

    gmp_printf ("sigma = %.4Ff\n", SIGMA);
    print_fmpz("ALPHA", ALPHA);
    print_fmpz("OMEGA", OMEGA);
    print_fmpz("IALPHA", IALPHA);
    print_fmpz("IOMEGA", IOMEGA);
    print_fmpz("IDOS", IDOS);

    printf("Bits of Q = %zu\n", mpz_sizeinbase (Q, 2));

    mpf_t t;
    mpf_init_set_ui(t,TAILCUT);

    compute_rectangles(NUMBEROFRECTANGLES, t, SIGMA, X, XBar, YBar);

    // mpz_p_tofile(XBar,NUMBEROFRECTANGLES+1,FILEZIGGURAT);
    // mpf_p_tofile(YBar,NUMBEROFRECTANGLES+1,FILEZIGGURAT);

    MD_CTX = EVP_MD_CTX_new();
    MD_SHAKE = EVP_get_digestbyname(XOF);
    MD_COM_AUX = EVP_get_digestbyname(MESSAGEDIGEST);

    #ifdef DEBUG
        if(!MD_COM_AUX) {
            red();
            printf("Unknown message digest");
            reset();
            exit(0);
        }
    #endif

    md_len_com = EVP_MD_get_size(MD_COM_AUX);

    if (strcmp(argv[1],"keygen") == 0){
        test_scheme_keygen(lambdas, degrees, Ds, Q);
    } else if (strcmp(argv[1],"commitment") == 0){
        test_scheme_commitment(lambdas, degrees, Ds, Q);
    }else if (strcmp(argv[1],"verifier") == 0){
        test_scheme_verifier(lambdas, degrees, Ds, Q);
    } else if (strcmp(argv[1],"opening") == 0){
        test_knowledge_valid_opening(lambdas, degrees, Ds, Q);
    } else if (strcmp(argv[1],"linear") == 0){
        test_knowledge_linear_relation(lambdas, degrees, Ds, Q);
    } else if (strcmp(argv[1],"multiplicative") == 0){
        test_knowledge_multiplicative_relation(lambdas, degrees, Ds, Q);
    }


    mpf_p_free_array(X,NUMBEROFRECTANGLES+1);
    mpf_p_free_array(YBar,NUMBEROFRECTANGLES+1);
    mpz_p_free_array(XBar,NUMBEROFRECTANGLES+1);

    mpz_p_free_array(moduli,number_params);
    mpz_p_free_array(Bs,number_params);
    mpf_p_free_array(sigmas,number_params);

    mpz_clear(BOUND);

    fmpz_clear(fQ);
    fmpz_clear(ALPHA);
    fmpz_clear(IALPHA);
    fmpz_clear(OMEGA);
    fmpz_clear(IOMEGA);
    fmpz_clear(IDOS);

    mpf_clear(SIGMA);
    mpf_clear(t);
    mpz_clear(Q);

    fmpz_clear_fft(modulus);
    fmpz_clear_fft(modulus_inv);

    EVP_MD_CTX_free(MD_CTX);

    return 1;
}
