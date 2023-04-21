#include <stdio.h>
#include <string.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod.h"
#include "flint/fmpz_mod_poly.h"
#include "param.h"
#include "aux.h"
#include "sampler.h"
#include "permutations.h"

int verifier1_aux(const fmpz_t alpha, mpz_t * gs[KAPPA + 1], const unsigned long * e_tilde[KAPPA + 1], EVP_MD_CTX * mdctx_com){

    #ifdef VERBOSE
        printf("\n **** VERIFIER (B = 1) **** \n \n");
    #endif

    mpz_t alpha2;
    mpz_init(alpha2);
    fmpz_get_mpz(alpha2, alpha);

    int sum = 0;
    int res;

    for(int j = 0; j < KAPPA + 1; j++){
        sum += vector_check(e_tilde[j]);
    }

    if(sum != KAPPA + 1){
        #ifdef DEBUG
            red();
            printf("Error in the verification (permutations of e' not correct)\n");
            reset();
        #endif
        res = 0;
    } else {
        #ifdef DEBUG
            green();
            printf("Permutations of e' correct\n");
            reset();
        #endif
        res = 1;
    }

    for(int j = 0; j < KAPPA + 1; j++){
        for(int i = 0; i < 2*N*K; i++){
            if (get_bit(e_tilde[j],i)){
                mpz_sub(gs[j][i],gs[j][i],alpha2);
            }
        }
        hash_vector(mdctx_com, gs[j], 2*N*K);
    }

    for(int j = 0; j < KAPPA + 1; j++){
        if(!EVP_DigestUpdate(mdctx_com, e_tilde[j], N*K/4)){
            #ifdef DEBUG
                printf("Error computing hash");
            #endif
            exit(0);
        }
    }

    mpz_clear(alpha2);

    return res;

}

int verifier0_aux(const unsigned char seed[], const fmpz_t alpha, mpz_t * gs[KAPPA + 1], const fmpz_mod_poly_t ** y, const fmpz_mod_poly_t * s, const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t ** C, EVP_MD_CTX *mdctx_com, fmpz_mod_poly_t ** gs_inv, int * pos) {
    unsigned char * buffer = (unsigned char *) malloc((size_t)(KAPPA + 1)*LAMB * sizeof(unsigned char));
    #pragma omp critical
    {
        EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
        EVP_DigestUpdate(MD_CTX, seed, LAMB);
        EVP_DigestFinalXOF(MD_CTX, buffer, (KAPPA + 1)*LAMB);
        EVP_MD_CTX_reset(MD_CTX);
    }

    EVP_DigestUpdate(mdctx_com, seed, LAMB); //comit to seed


    for(int j = 0; j < K; j++) {
        hash_poly_FFT(mdctx_com, y[j]);
    }

    mpz_t ** gs_inverse = gs;


    for(int i = 0;i < KAPPA + 1; i++){
        permuteInv(gs_inverse[i],2*N*K,buffer+i*LAMB);
    }

    mpz_t * f_scaled_sum = (mpz_t *)malloc(N*K * sizeof(mpz_t));
    mpz_p_init_array(f_scaled_sum,N*K);
    for(int i = KAPPA; i >= 0; i--){
        for (int j = 0; j < N*K; j++){
            mpz_mul_2exp(f_scaled_sum[j],f_scaled_sum[j],1);
        }
        mpz_p_add(gs_inverse[i],f_scaled_sum,N*K,N*K,f_scaled_sum);
    }

    fmpz_mod_poly_t ** bs;
    fmpz_init_array_fft(&bs,K);

    for(int i = 0; i < K; i++){
        fmpz_mod_poly_mult_FFT(bs[i],B[i],s);
    }

    fmpz_mod_poly_t ** c_bar_scaled;
    fmpz_init_array_fft(&c_bar_scaled,K);

    fmpz_mod_poly_t * aux;
    fmpz_init_fft(&aux);

    mpz_t powkappa;
    mpz_init_set_ui(powkappa, 1);
    mpz_mul_2exp(powkappa,powkappa,KAPPA);


    for(int i = 0; i < N; i++){
        fmpz_mod_poly_set_coeff_mpz(aux[i/(N/D)], i%(N/D), powkappa, ctx);
    }

    fmpz_p_fcr(aux,0,D,ALPHA,OMEGA,D,fQ);

    for(int i = 0; i < K; i++){
        fmpz_mod_poly_add_FFT(c_bar_scaled[i], C[i], aux);
        fmpz_mod_poly_scalar_mul_fmpz_FFT(c_bar_scaled[i], c_bar_scaled[i], alpha);
        fmpz_mod_poly_sub_FFT(c_bar_scaled[i], c_bar_scaled[i], bs[i]);
    }



    for(int i = 0; i < K; i++){
        for(int j = 0; j < N; j++){
            fmpz_mod_poly_set_coeff_mpz(gs_inv[i][j/(N/D)], j%(N/D), f_scaled_sum[j+i*N], ctx);
        }
        fmpz_p_fcr(gs_inv[i],0,D,ALPHA,OMEGA,D,fQ);
    }

    for(int i = 0; i < K; i++){
        fmpz_mod_poly_sub_FFT(gs_inv[i], y[i], gs_inv[i]);
        fmpz_mod_poly_add_FFT(gs_inv[i], gs_inv[i], c_bar_scaled[i]);
    }

    int isLatticePoint = lattice_membership(gs_inv,A,pos);

    free(buffer);
    mpz_p_free_array(f_scaled_sum,N*K);
    free(f_scaled_sum);
    fmpz_clear_array_fft(bs,K);
    fmpz_clear_array_fft(c_bar_scaled,K);
    fmpz_clear_fft(aux);

    mpz_clear(powkappa);

    if(isLatticePoint){
        #ifdef DEBUG
            green();
            printf("Test of lattice membership passed \n");
            reset();
        #endif
        return 1;
    } else {
        #ifdef DEBUG
            red();
            printf("Test of lattice membership not passed \n");
            reset();
        #endif
        return 0;
    }

}


int prover_initial_aux(const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t * r, const unsigned long ** e_primes, EVP_MD_CTX * mdctx_com1, EVP_MD_CTX * mdctx_com2, mpz_t ** f_permuted, unsigned char * seed, unsigned long ** e_tilde, fmpz_mod_poly_t * mu, fmpz_mod_poly_t * rho, fmpz_mod_poly_t ** y){

    /* ------------------------------------------------------------------------------------------------------ */
    /* --------------------------- (1) Sample kappa + 1 permutations from G_{2nk} --------------------------- */
    /* ------------------------------------------------------------------------------------------------------ */
    #ifdef VERBOSE
        clock_t begin = clock();
    #endif

    RAND_priv_bytes(seed, LAMB);
    unsigned char * buffer = (unsigned char *) malloc((size_t)(KAPPA + 1)*LAMB * sizeof(unsigned char));
    #pragma omp critical
    {
        EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
        EVP_DigestUpdate(MD_CTX, seed, LAMB);
        EVP_DigestFinalXOF(MD_CTX, buffer, (KAPPA + 1)*LAMB);
        EVP_MD_CTX_reset(MD_CTX);
    }

    EVP_DigestUpdate(mdctx_com1, seed, LAMB); //commit to seed

    #ifdef VERBOSE
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Sample kappa + 1 permutations from G_{2nk} = %f\n", time_spent);
    #endif

    /* ------------------------------------------------------------------------------------------------------ */
    /* ----------------------------- (2) Sample kappa + 1 vectors from Zq^{2nk} ----------------------------- */
    /* ------------------------------------------------------------------------------------------------------ */

    #ifdef VERBOSE
        begin = clock();
    #endif

    mpz_t ** f = f_permuted;


    #ifdef VERBOSE
        printf("KAPPA = %zu\n", KAPPA);
    #endif

    for(int i = 0; i < KAPPA + 1; i++){
      sampler_zq_list(f[i],2*N*K);
      #ifdef VERBOSE
          mpz_p_print(f[i],2*N*K,"");
      #endif
    }

    #ifdef VERBOSE
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Sample kappa + 1 vectors from Zq^{2nk} = %f\n", time_spent);
    #endif


    /* ------------------------------------------------------------------------------------------------------ */
    /* ----------------------- Computation of a * mu + b * rho + I' * sum(2^j * f_j) ------------------------ */
    /* ------------------------------------------------------------------------------------------------------ */

    /* (3) Sample mu and rho uniformly at random from R_q */

    fmpz_sampler_zq_poly_FFT(rho);


    #ifdef VERBOSE
        printf("\nComputation of I' * sum(2^j * f_j)\n");
    #endif

    /* Directly compute scaled sum */

    #ifdef VERBOSE
        begin = clock();
    #endif

    mpz_t * f_scaled_sum = (mpz_t *)malloc(N*K * sizeof(mpz_t));
    mpz_p_init_array(f_scaled_sum,N*K);


    for(int i = KAPPA; i >= 0; i--){
      for (int j = 0; j < N*K; j++){
        mpz_mul_2exp(f_scaled_sum[j],f_scaled_sum[j],1);
      }
      mpz_p_add(f[i],f_scaled_sum,N*K,N*K,f_scaled_sum);
    }
    #ifdef VERBOSE
        mpz_p_print(f_scaled_sum,N*K,"\nI' * sum(2^j * f_j)");
    #endif

    for(int i = 0; i < K; i++){
        for(int j = 0; j < N; j++){
            fmpz_mod_poly_set_coeff_mpz(y[i][j/(N/D)], j%(N/D), f_scaled_sum[j+i*N], ctx);
        }
        fmpz_p_fcr(y[i],0,D,ALPHA,OMEGA,D,fQ);
    }

    mpz_p_free_array(f_scaled_sum,N*K);
    free(f_scaled_sum);

    /* Compute A*mu + B*rho + I' * sum(2^j * f_j) */

    fmpz_mod_poly_t * aux;
    fmpz_init_fft(&aux);
    for(int i = 0; i < K; i++){
        fmpz_mod_poly_mult_FFT(aux, mu, A[i]);
        fmpz_mod_poly_add_FFT(y[i], y[i], aux);
        fmpz_mod_poly_mult_FFT(aux, rho, B[i]);
        fmpz_mod_poly_add_FFT(y[i], y[i], aux);
    }
    fmpz_clear_fft(aux);

    #ifdef VERBOSE
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Computation of a * mu + b * rho + I' * sum(2^j * f_j) = %f\n", time_spent);
    #endif

    #ifdef VERBOSE
        begin = clock();
    #endif

    for(int j = 0; j < K; j++) {
        hash_poly_FFT(mdctx_com1, y[j]);
    }

    #ifdef VERBOSE
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Hashing k vectors from Zq^{n} = %f\n", time_spent);
    #endif


    /* ------------------------------------------------------------------------------------------------------ */
    /* --------------------------------------- Compute and commit to pi(f_j) ---------------------------------- */
    /* ------------------------------------------------------------------------------------------------------ */

    #ifdef VERBOSE
        begin = clock();
    #endif


    for(int i = 0; i < KAPPA + 1; i++){
      permute(f_permuted[i],2*N*K,buffer+i*LAMB);
    }

    #ifdef VERBOSE
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Permute kappa + 1 vectors from Zq^{2nk} = %f\n", time_spent);
    #endif

    #ifdef VERBOSE
        begin = clock();
    #endif

    for(int i = 0; i < KAPPA + 1; i++){
        hash_vector(mdctx_com2, f_permuted[i], 2*N*K);
    }

    #ifdef VERBOSE
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Hashing kappa + 1 vectors from Zq^{2nk} = %f\n", time_spent);
    #endif


    /* ------------------------------------------------------------------------------------------------------ */
    /* ---------------------------- Computation of pi(e'_j) for j in {0,...,KAPPA} ------------------------------ */
    /* ------------------------------------------------------------------------------------------------------ */

    #ifdef VERBOSE
        printf("\n EXTENSION OF THE ERROR PERMUTED \n");
    #endif

    #ifdef VERBOSE
        begin = clock();
    #endif

    for(int i = 0; i < KAPPA + 1; i++){
        memcpy(e_tilde[i],e_primes[i],N*K/4);
    }

    for(int j = 0; j < KAPPA + 1; j++){
        permuteBit(e_tilde[j],2*N*K,buffer+j*LAMB);
        EVP_DigestUpdate(mdctx_com2, e_tilde[j], N*K/4);
    }

    #ifdef VERBOSE
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Computation and hashing of pi_j(e'_j)  = %f\n", time_spent);
    #endif

    free(buffer);

    return 1;
}

int prover_middle_aux(const fmpz_t alpha, const unsigned long ** e_tilde, mpz_t ** f_permuted){

    /* ------------------------------------------------------------------------------------------------------ */
    /* --------------------------------------- Compute pi(f_j + alpha e_j) ---------------------------------- */
    /* ------------------------------------------------------------------------------------------------------ */
    #ifdef VERBOSE
        clock_t begin = clock();
    #endif

    mpz_t alpha2;
    mpz_init(alpha2);
    fmpz_get_mpz(alpha2, alpha);


    for(int j = 0; j < KAPPA + 1; j++){


        for(int i = 0; i < 2*K*N; i++){
            if (get_bit(e_tilde[j],i)) {
                mpz_add(f_permuted[j][i],f_permuted[j][i],alpha2);
                mpz_mod(f_permuted[j][i],f_permuted[j][i],Q);
            }
        }


        #ifdef VERBOSE
            mpz_p_print(f_permuted[j],2*N*K,"\n");
        #endif
    }


    #ifdef VERBOSE
        printf("pi(f_j + alpha e_j) \n");
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Compute pi(f_j + alpha e_j) = %f\n", time_spent);
    #endif

    mpz_clear(alpha2);

    return 1;
}

int prover0_aux(const fmpz_mod_poly_t * r, fmpz_t alpha, fmpz_mod_poly_t * rho, fmpz_mod_poly_t * rho_r_scaled){

    /* ------------------------------------------------------------------------------------------------------ */
    /* --------------------------------------- Compute rho_r_scaled ----------------------------------------- */
    /* ------------------------------------------------------------------------------------------------------ */
    #ifdef VERBOSE
        clock_t begin = clock();
    #endif

    fmpz_mod_poly_scalar_mul_fmpz_FFT(rho_r_scaled, r, alpha);
    fmpz_mod_poly_add_FFT(rho_r_scaled, rho_r_scaled, rho);

    #ifdef VERBOSE
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Compute rho + alpha*r = %f\n", time_spent);
    #endif

    return 1;
}
