#include <stdio.h>
#include "aux.h"
#include "param.h"
#include "sampler.h"
#include "auxzkp.h"
#include <omp.h>

// constant that determines the ratio of randomness sampled each time we requiere sampling uniformly at random
const double expansion = 2.1;

/**********************************************************************
* Name:  alphas
*
* Description:
*    given numCom arrays of delta auxiliary commitments uses them to determine alphaSeed
*    via a XOF in order to use it to obtain len fmpz_t that are stored in alpha
*
* Arguments (read only):
*    int delta: repetitions to obtain soundness, is the size of the arrays of commitments
*    int len: number of alphas to obtain, it can be delta or 2*delta if we need alphas and betas
*    const fmpz_mod_poly_t ** statement: statement of the proof
*    int stLen: number of polynomials in the statement
*    const unsigned char *** aux_coms: numCom arrays of delta auxiliary commitments that determine the alphas
*    int numCom: number of arrays of auxiliary commitments that determine the alphas
*
* Populates:
*    fmpz_t * alpha: array of len alphas deterministacally obtained pseudorandomly from aux_coms via alphaSeed
*    unsigned char * alphaSeed: seed deterministacally obtained pseudorandomly from aux_coms used to get the alphas
*
***********************************************************************/
void alphas(fmpz_t * alpha, unsigned char * alphaSeed, int SEEDS, int delta, int len, const fmpz_mod_poly_t ** statement, int stLen, const unsigned char *** aux_coms, int numCom) {
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    unsigned char round = 0;
    EVP_DigestUpdate(MD_CTX, &round, sizeof(unsigned char));
    for(int j = 0; j < stLen; j++) {
        hash_poly_FFT(MD_CTX, statement[j]);
    }
    for (int i=0; i<numCom; i++) {
        for (int rep = 0; rep < delta; rep++){
            EVP_DigestUpdate(MD_CTX, aux_coms[i][rep], md_len_com);
        }
    }
    EVP_DigestFinalXOF(MD_CTX, alphaSeed, SEEDS);
    EVP_MD_CTX_reset(MD_CTX);

    size_t alphaBufferLength = (size_t) delta*BQ*expansion; // this should be optimized
    unsigned char * alphaBuffer = (unsigned char *) malloc(alphaBufferLength * sizeof(unsigned char));
    size_t alphaBufferPosition = 0;
    unsigned char alphaCounter = 0;
    refreshBuffer(alphaBuffer,alphaBufferLength,alphaSeed,&alphaCounter);

    for (int rep = 0; rep < len; rep++){
        fmpz_init(alpha[rep]);
        fmpz_sampler_zq_B(alpha[rep], alphaBuffer, alphaBufferLength, &alphaBufferPosition, alphaSeed, &alphaCounter);
    }

    free(alphaBuffer);
}

/**********************************************************************
* Name:  prover_opening
*
* Description:
*    given a commitment key and a commitment opening computes a nizkpok of such opening
*
* Arguments (read only):
*    const fmpz_mod_poly_t ** A, B: commitment key, a pair of K polynomials
*    const fmpz_mod_poly_t ** C: commitment A * m + B * r + e, K polynomials
*    const fmpz_mod_poly_t * m: commited message, a polynomial
*    const fmpz_mod_poly_t * r: commitment randomness, a uniformly random polynomial
*    const mpz_t * e: commitment randomness error, coefficients of a vector of polynomials following a discrete gaussian bounded by B
*    int delta: repetitions to obtain soundness
*
* Populates:
* - Common elements:
*    unsigned char ** md_value_com1: array of delta auxiliary commitments com1 = Hash( seed || y || d1 )
*    unsigned char ** md_value_com2: array of delta auxiliary commitments com2 = Hash( fj permuted || e' permuted || d2 )
*    mpz_t *** g: array of delta vectors of coefficients of K polynomials pi_j(fj+alpha·e'j)
* - Answers for verifier 0:
*    unsigned char ** seed: array of delta seeds that determine the permutations pi_j
*    fmpz_mod_poly_t *** y: array of delta vectors of polynomials A * mu + B * rho + I' * sum(2^i * f_i)
*    fmpz_mod_poly_t ** s: array of delta polynomials rho + alpha * r
*    unsigned char ** d1: array of delta openings for md_value_com1
* - Answers for verifier 1:
*    unsigned long *** e_tilde: array of delta e', permuted bits of the errors
*    unsigned char ** d2: array of delta openings for md_value_com2
*
***********************************************************************/
int prover_opening(const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t ** C, const fmpz_mod_poly_t * m, const fmpz_mod_poly_t * r, const mpz_t * e, int delta, unsigned char ** md_value_com1, unsigned char ** md_value_com2, mpz_t *** g, unsigned char ** seed,  fmpz_mod_poly_t *** y, fmpz_mod_poly_t ** s, unsigned char ** d1, unsigned long *** e_tilde, unsigned char ** d2){

    EVP_MD_CTX *mdctx_com1, *mdctx_com2;

    fmpz_mod_poly_t * mu;

    fmpz_mod_poly_t ** rho;
    fmpz_init_array_fft(&rho,delta);

    /* Compute extensions e_prime_j, i.e, compute e_prime_j = (e_bar_j || w) such that e_prime_j
    has the same number of 0's and 1's */

    unsigned long ** e_primes = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
    for(int i = 0; i < KAPPA + 1; i++){
        e_primes[i] = (unsigned long *) malloc(N*K/4);
    }
    expand(e, e_primes);

    #pragma omp parallel private(mdctx_com1,mdctx_com2,mu)
    {
        mdctx_com1 = EVP_MD_CTX_new();
        mdctx_com2 = EVP_MD_CTX_new();

        fmpz_init_fft(&mu);

        #pragma omp for
        for (int rep = 0; rep < delta; rep++){
            EVP_DigestInit_ex2(mdctx_com1, MD_COM_AUX, NULL);
            EVP_DigestInit_ex2(mdctx_com2, MD_COM_AUX, NULL);

            fmpz_sampler_zq_poly_FFT(mu);

            prover_initial_aux(A, B, r, e_primes, mdctx_com1, mdctx_com2, g[rep], seed[rep], e_tilde[rep], mu, rho[rep], y[rep]);

            RAND_priv_bytes(d1[rep],LAMB);
            EVP_DigestUpdate(mdctx_com1, d1[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com1, md_value_com1[rep], NULL);

            #ifdef VERBOSE
                printf("%3d - ",rep);
                print_hex("Prover (com1)", md_value_com1[rep], md_len_com);
                printf("%3d - ",rep);
                print_hex("Prover (d1)", d1[rep], LAMB);
            #endif


            RAND_priv_bytes(d2[rep],LAMB);
            EVP_DigestUpdate(mdctx_com2, d2[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com2, md_value_com2[rep], NULL);

            #ifdef VERBOSE
                printf("%3d - ",rep);
                print_hex("Prover (com2)", md_value_com2[rep], md_len_com);
                printf("%3d - ",rep);
                print_hex("Prover (d2)", d2[rep], LAMB);
            #endif

            EVP_MD_CTX_reset(mdctx_com1);
            EVP_MD_CTX_reset(mdctx_com2);
        }
        EVP_MD_CTX_free(mdctx_com1);
        EVP_MD_CTX_free(mdctx_com2);

        fmpz_clear_fft(mu);
    }

    for(int i = 0; i < KAPPA + 1; i++){
        free(e_primes[i]);
    }
    free(e_primes);




    // Obtain first challenge from initial commitments via a XOF
    fmpz_t * alpha = (fmpz_t *) malloc(delta*sizeof(fmpz_t));
    unsigned char alphaSeed[SEEDOL];
    const unsigned char ** aux_coms[2] = {(const unsigned char **) md_value_com1, (const unsigned char **) md_value_com2};
    const fmpz_mod_poly_t * statement[3*K];
    for (int i = 0; i < K; i++){
        statement[i] = A[i];
        statement[i+K] = B[i];
        statement[i+2*K] = C[i];
    }
    alphas(alpha, alphaSeed, SEEDOL, delta, delta, statement, 3*K, aux_coms, 2);

    #ifdef VERBOSE
        print_hex("Prover (alphaSeed)", alphaSeed, LAMB);
    #endif

    #pragma omp parallel for
    for (int rep = 0; rep < delta; rep++){
        prover_middle_aux(alpha[rep], e_tilde[rep], g[rep]);
    }

    // Obtain second challenge from previous messages via a XOF
    unsigned char round = 1;
    unsigned char bSeed[SEEDOL];
    unsigned char * secondChallenge = malloc((int) ceil(delta/8.0)*sizeof(unsigned char));
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, &round, sizeof(unsigned char));

    for(int j = 0; j < 3*K; j++) {
        hash_poly_FFT(MD_CTX, statement[j]);
    }

    EVP_DigestUpdate(MD_CTX, alphaSeed, SEEDOL);

    for (int rep = 0; rep < delta; rep++){
        for(int j = 0; j < KAPPA + 1; j++){
            hash_vector(MD_CTX, g[rep][j], 2*N*K);
        }

    }
    EVP_DigestFinalXOF(MD_CTX, bSeed, SEEDOL);
    EVP_MD_CTX_reset(MD_CTX);
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, bSeed, SEEDOL);
    EVP_DigestFinalXOF(MD_CTX, secondChallenge, (int) ceil(delta/8.0));
    EVP_MD_CTX_reset(MD_CTX);


    #ifdef VERBOSE
        print_hex("Prover (secondChallenge)", secondChallenge, (int) ceil(delta/8.0));
    #endif

    #pragma omp parallel for
    for (int rep = 0; rep < delta; rep++){
        if (!(secondChallenge[rep/8] & (1 << rep%8))) {
            prover0_aux(r, alpha[rep], rho[rep], s[rep]);
        }
    }

    #pragma omp parallel for
    for (int rep = 0; rep < delta; rep++){
        fmpz_clear(alpha[rep]);
    }
    free(alpha);
    fmpz_clear_array_fft(rho,delta);
    free(secondChallenge);

    return 1;
}

/**********************************************************************
* Name:  verifier_opening
*
* Description:
*    given a commitment key, a commitment and a nizkpok of a valid opening verifies the latter
*
* Arguments (read only, EXCEPT elements of g that are overwritten during the verification):
*    const fmpz_mod_poly_t ** A, B: commitment key, a pair of K polynomials
*    const fmpz_mod_poly_t ** C: commitment A * m + B * r + e, K polynomials
*    int delta: repetitions to obtain soundness
* - Common elements:
*    unsigned char ** com1: delta auxiliary commirments com1 = Hash( seed || y || d1 )
*    unsigned char ** com2: delta auxiliary commirments com2 = Hash( fj permuted || e' permuted || d2 )
*    mpz_t *** g: array of delta vectors of coefficients of K polynomials pi_j(fj+alpha·e'j)
* - Answers for verifier 0:
*    unsigned char ** seed: array of delta seeds that determine the permutations pi_j
*    fmpz_mod_poly_t *** y: array of delta vectors of polynomials A * mu + B * rho + I' * sum(2^i * f_i)
*    fmpz_mod_poly_t ** s: array of delta polynomials rho + alpha * r
*    unsigned char ** d1: array of delta openings for com1
* - Answers for verifier 1:
*    unsigned long *** e_tilde: array of delta e', permuted bits of the errors
*    unsigned char ** d2: array of delta openings for com2
* - Auxiliary elements:
*    const int * pos: position of elements of A invertible modulo each of the factors of x^n+1
*
***********************************************************************/
int verifier_opening(const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t ** C, int delta, const unsigned char ** com1,const unsigned char ** com2, mpz_t *** g, const unsigned char ** seed, const fmpz_mod_poly_t *** y, const fmpz_mod_poly_t ** s, const unsigned char ** d1, const unsigned long *** e_tilde, const unsigned char ** d2, const int * pos){

    // Obtain first challenge from initial commitments via a XOF
    fmpz_t * alpha = (fmpz_t *) malloc(delta*sizeof(fmpz_t));
    unsigned char alphaSeed[SEEDOL];
    const unsigned char ** aux_coms[2] = {com1, com2};
    const fmpz_mod_poly_t * statement[3*K];
    for (int i = 0; i < K; i++){
        statement[i] = A[i];
        statement[i+K] = B[i];
        statement[i+2*K] = C[i];
    }
    alphas(alpha, alphaSeed, SEEDOL, delta, delta, statement, 3*K, aux_coms, 2);

    #ifdef VERBOSE
        print_hex("Verifier (alphaSeed)", alphaSeed, LAMB);
    #endif

    // Obtain second challenge from previous messages via a XOF
    unsigned char round = 1;
    unsigned char bSeed[SEEDOL];
    unsigned char * secondChallenge = malloc((int) ceil(delta/8.0)*sizeof(unsigned char));
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, &round, sizeof(unsigned char));
    for(int j = 0; j < 3*K; j++) {
        hash_poly_FFT(MD_CTX, statement[j]);
    }
    EVP_DigestUpdate(MD_CTX, alphaSeed, SEEDOL);
    for (int rep = 0; rep < delta; rep++){
        for(int j = 0; j < KAPPA + 1; j++){
            hash_vector(MD_CTX, g[rep][j], 2*N*K);
        }
    }
    EVP_DigestFinalXOF(MD_CTX, bSeed, SEEDOL);
    EVP_MD_CTX_reset(MD_CTX);
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, bSeed, SEEDOL);
    EVP_DigestFinalXOF(MD_CTX, secondChallenge, (int) ceil(delta/8.0));
    EVP_MD_CTX_reset(MD_CTX);

    #ifdef VERBOSE
        print_hex("Verifier (secondChallenge)", secondChallenge, (int) ceil(delta/8.0));
    #endif

    int res = 1;

    EVP_MD_CTX *mdctx_com;
    unsigned char md_value_com[md_len_com];

    fmpz_mod_poly_t ** gs_inv;

    #pragma omp parallel private(mdctx_com,md_value_com,gs_inv) reduction(&:res)
    {
        mdctx_com = EVP_MD_CTX_new();
        fmpz_init_array_fft(&gs_inv,K);

        #pragma omp for
        for (int rep = 0; rep < delta; rep++){
            EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
            if (!(secondChallenge[rep/8] & (1 << rep%8))) {
                /* verifier 0 */
                res &= verifier0_aux(seed[rep], alpha[rep], g[rep], y[rep], s[rep], A, B, C, mdctx_com, gs_inv, pos);
                EVP_DigestUpdate(mdctx_com, d1[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                res &= verify_aux_com("com1", com1[rep], md_value_com, md_len_com);
                #ifdef VERBOSE
                    printf("%3d - ",rep);
                    print_hex("Verifier (com1)", md_value_com, md_len_com);
                    printf("%3d - ",rep);
                    print_hex("Verifier (d1)", d1[rep], LAMB);
                #endif
            } else {
                /* verifier 1 */
                res &= verifier1_aux(alpha[rep], g[rep], e_tilde[rep], mdctx_com);
                EVP_DigestUpdate(mdctx_com, d2[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                res &= verify_aux_com("com2", com2[rep], md_value_com, md_len_com);
                #ifdef VERBOSE
                    printf("%3d - ",rep);
                    print_hex("Verifier (com2)", md_value_com, md_len_com);
                    printf("%3d - ",rep);
                    print_hex("Verifier (d2)", d2[rep], LAMB);
                #endif
            }
            EVP_MD_CTX_reset(mdctx_com);
        }
        fmpz_clear_array_fft(gs_inv,K);
        EVP_MD_CTX_free(mdctx_com);
    }

    for (int rep = 0; rep < delta; rep++){
        fmpz_clear(alpha[rep]);
    }
    free(alpha);
    free(secondChallenge);

    return res;
}

/**********************************************************************
* Name:  prover_linear
*
* Description:
*    given a commitment key and 3 commitment openings satisfying a linear relation computes a nizkpok of such openings
*
* Arguments (read only):
*    const fmpz_mod_poly_t ** A, B: commitment key, a pair of K polynomials
*    const fmpz_mod_poly_t ** Ci: i-th commitment A * m + B * r + e, K polynomials
*    const fmpz_mod_poly_t * mi: i-th commited message, a polynomial
*    const fmpz_mod_poly_t * ri: i-th commitment randomness, a uniformly random polynomial
*    const mpz_t * ei: i-th commitment randomness error, coefficients of a vector of polynomials following a discrete gaussian bounded by B
*    const fmpz_mod_poly_t * lambda1, lambda2: linear coefficients, a pair of polynomials
*    int delta: repetitions to obtain soundness
*
* Populates:
* - Common elements:
*    unsigned char ** md_value_com1: array of delta auxiliary commitments com1 = Hash( {seedi} || {yi} || d1 )
*    unsigned char ** md_value_com2: array of delta auxiliary commitments com2 = Hash( {fij permuted} || {e'ij} permuted || d2 )
*    mpz_t *** gi: i-th array of delta vectors of coefficients of K polynomials pi_ij(fij+alpha·e'ij)
* - Answers for verifier 0:
*    unsigned char ** seedi: i-th array of delta seeds that determine the permutations pi_j
*    fmpz_mod_poly_t *** yi: i-th  array of delta vectors of polynomials A * mui + B * rhoi + I' * sum(2^l * f_ij)
*    fmpz_mod_poly_t ** si: i-th  array of delta polynomials rhoi + alpha * ri
*    unsigned char ** d1: array of delta openings for md_value_com1
* - Answers for verifier 1:
*    unsigned long *** e_tildei: i-th array of delta e'i, permuted bits of the errors
*    unsigned char ** d2: array of delta openings for md_value_com2
*
***********************************************************************/
int prover_linear(const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t ** C1,  const fmpz_mod_poly_t ** C2,  const fmpz_mod_poly_t ** C3, const fmpz_mod_poly_t * m1, const fmpz_mod_poly_t * m2, const fmpz_mod_poly_t * m3, const fmpz_mod_poly_t * r1, const fmpz_mod_poly_t * r2, const fmpz_mod_poly_t * r3, const mpz_t * e1, const mpz_t * e2, const mpz_t * e3, const fmpz_mod_poly_t * lambda1, const fmpz_mod_poly_t * lambda2, int delta, unsigned char ** md_value_com1, unsigned char ** md_value_com2, mpz_t *** g1, mpz_t *** g2, mpz_t *** g3, unsigned char ** seed1, unsigned char ** seed2, unsigned char ** seed3, fmpz_mod_poly_t *** y1, fmpz_mod_poly_t *** y2, fmpz_mod_poly_t *** y3, fmpz_mod_poly_t ** s1, fmpz_mod_poly_t ** s2, fmpz_mod_poly_t ** s3, unsigned char ** d1, unsigned long *** e_tilde1, unsigned long *** e_tilde2, unsigned long *** e_tilde3, unsigned char ** d2) {

    EVP_MD_CTX *mdctx_com1, *mdctx_com2;


    fmpz_mod_poly_t * mu1, * mu2, * mu3;
    fmpz_mod_poly_t ** rho1, ** rho2, ** rho3;
    fmpz_init_array_fft(&rho1,delta);
    fmpz_init_array_fft(&rho2,delta);
    fmpz_init_array_fft(&rho3,delta);

    fmpz_mod_poly_t * mu_lambda_aux;

    unsigned long ** e_primes1 = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
    unsigned long ** e_primes2 = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
    unsigned long ** e_primes3 = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
    for(int i = 0; i < KAPPA + 1; i++){
        e_primes1[i] = (unsigned long *) malloc(N*K/4);
        e_primes2[i] = (unsigned long *) malloc(N*K/4);
        e_primes3[i] = (unsigned long *) malloc(N*K/4);
    }
    expand(e1, e_primes1);
    expand(e2, e_primes2);
    expand(e3, e_primes3);

    #pragma omp parallel private(mdctx_com1,mdctx_com2,mu1,mu2,mu3,mu_lambda_aux)
    {
        mdctx_com1 = EVP_MD_CTX_new();
        mdctx_com2 = EVP_MD_CTX_new();

        fmpz_init_fft(&mu1);
        fmpz_init_fft(&mu2);
        fmpz_init_fft(&mu3);

        fmpz_init_fft(&mu_lambda_aux);

        #pragma omp for
        for (int rep = 0; rep < delta; rep++){
            EVP_DigestInit_ex2(mdctx_com1, MD_COM_AUX, NULL);
            EVP_DigestInit_ex2(mdctx_com2, MD_COM_AUX, NULL);

            fmpz_sampler_zq_poly_FFT(mu1);
            fmpz_sampler_zq_poly_FFT(mu2);

            prover_initial_aux(A, B, r1, e_primes1, mdctx_com1, mdctx_com2, g1[rep], seed1[rep], e_tilde1[rep], mu1, rho1[rep], y1[rep]);
            prover_initial_aux(A, B, r2, e_primes2, mdctx_com1, mdctx_com2, g2[rep], seed2[rep], e_tilde2[rep], mu2, rho2[rep], y2[rep]);

            fmpz_mod_poly_mult_FFT(mu3, lambda1, mu1);
            fmpz_mod_poly_mult_FFT(mu_lambda_aux, lambda2, mu2);
            fmpz_mod_poly_add_FFT(mu3, mu3, mu_lambda_aux);

            prover_initial_aux(A, B, r3, e_primes3, mdctx_com1, mdctx_com2, g3[rep], seed3[rep], e_tilde3[rep], mu3, rho3[rep], y3[rep]);

            RAND_priv_bytes(d1[rep],LAMB);
            EVP_DigestUpdate(mdctx_com1, d1[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com1, md_value_com1[rep], NULL);

            RAND_priv_bytes(d2[rep],LAMB);
            EVP_DigestUpdate(mdctx_com2, d2[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com2, md_value_com2[rep], NULL);


            EVP_MD_CTX_reset(mdctx_com1);
            EVP_MD_CTX_reset(mdctx_com2);

        }

        fmpz_clear_fft(mu1);
        fmpz_clear_fft(mu2);
        fmpz_clear_fft(mu3);

        fmpz_clear_fft(mu_lambda_aux);

        EVP_MD_CTX_free(mdctx_com1);
        EVP_MD_CTX_free(mdctx_com2);
    }

    for(int i = 0; i < KAPPA + 1; i++){
        free(e_primes1[i]);
        free(e_primes2[i]);
        free(e_primes3[i]);
    }
    free(e_primes1);
    free(e_primes2);
    free(e_primes3);

    // Obtain first challenge from initial commitments via a XOF
    fmpz_t * alpha = (fmpz_t *) malloc(delta*sizeof(fmpz_t));
    unsigned char alphaSeed[SEEDOL];
    const unsigned char ** aux_coms[2] = {(const unsigned char **) md_value_com1, (const unsigned char **) md_value_com2};
    const fmpz_mod_poly_t * statement[5*K+2];
    for (int i = 0; i < K; i++){
        statement[i] = A[i];
        statement[i+K] = B[i];
        statement[i+2*K] = C1[i];
        statement[i+3*K] = C2[i];
        statement[i+4*K] = C3[i];
    }
    statement[5*K] = lambda1;
    statement[5*K+1] = lambda2;
    alphas(alpha, alphaSeed, SEEDOL, delta, delta, statement, 5*K+2, aux_coms, 2);

    #pragma omp parallel for
    for (int rep = 0; rep < delta; rep++){
        prover_middle_aux(alpha[rep], e_tilde1[rep], g1[rep]);
        prover_middle_aux(alpha[rep], e_tilde2[rep], g2[rep]);
        prover_middle_aux(alpha[rep], e_tilde3[rep], g3[rep]);
    }

    // Obtain second challenge from previous messages via a XOF
    unsigned char round = 1;
    unsigned char bSeed[SEEDOL];
    unsigned char * secondChallenge = malloc((int) ceil(delta/8.0)*sizeof(unsigned char));
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, &round, sizeof(unsigned char));
    for(int j = 0; j < 5*K+2; j++) {
        hash_poly_FFT(MD_CTX, statement[j]);
    }
    EVP_DigestUpdate(MD_CTX, alphaSeed, SEEDOL);
    for (int rep = 0; rep < delta; rep++){
        for(int j = 0; j < KAPPA + 1; j++){
            hash_vector(MD_CTX, g1[rep][j], 2*N*K);
            hash_vector(MD_CTX, g2[rep][j], 2*N*K);
            hash_vector(MD_CTX, g3[rep][j], 2*N*K);
        }
    }
    EVP_DigestFinalXOF(MD_CTX, bSeed, SEEDOL);
    EVP_MD_CTX_reset(MD_CTX);
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, bSeed, SEEDOL);
    EVP_DigestFinalXOF(MD_CTX, secondChallenge, (int) ceil(delta/8.0));
    EVP_MD_CTX_reset(MD_CTX);

    #pragma omp parallel for
    for (int rep = 0; rep < delta; rep++){
        if (!(secondChallenge[rep/8] & (1 << rep%8))) {
            prover0_aux(r1, alpha[rep], rho1[rep], s1[rep]);
            prover0_aux(r2, alpha[rep], rho2[rep], s2[rep]);
            prover0_aux(r3, alpha[rep], rho3[rep], s3[rep]);
        }
    }

    for (int rep = 0; rep < delta; rep++){
        fmpz_clear(alpha[rep]);
    }
    free(alpha);
    free(secondChallenge);
    fmpz_clear_array_fft(rho1,delta);
    fmpz_clear_array_fft(rho2,delta);
    fmpz_clear_array_fft(rho3,delta);

    return 1;
}

/**********************************************************************
* Name:  verifier_linear
*
* Description:
*    given a commitment key, 3 commitments and a nizkpok of 3 openings satisfying a linear relation verifies the latter
*
* Arguments (read only, EXCEPT elements of g that are overwritten during the verification):
*    const fmpz_mod_poly_t ** A, B: commitment key, a pair of K polynomials
*    const fmpz_mod_poly_t ** Ci: i-th commitment A * m + B * r + e, K polynomials
*    const fmpz_mod_poly_t * lambda1, lambda2: linear coefficients, a pair of polynomials
*    int delta: repetitions to obtain soundness
* - Common elements:
*    unsigned char ** com1: delta auxiliary commirments com1 = Hash( {seedi} || {yi} || d1 )
*    unsigned char ** com2: delta auxiliary commirments com2 = Hash( {fij permuted} || {e'ij permuted} || d2 )
*    mpz_t *** gi: i-th array of delta vectors of coefficients of K polynomials pi_ij(fij+alpha·e'ij)
* - Answers for verifier 0:
*    unsigned char ** seedi: i-th array of delta seeds that determine the permutations pi_j
*    fmpz_mod_poly_t *** yi: i-th array of delta vectors of polynomials A * mui + B * rhoi + I' * sum(2^l * f_ij)
*    fmpz_mod_poly_t ** si: i-th array of delta polynomials rhoi + alpha * ri
*    unsigned char ** d1: array of delta openings for com1
* - Answers for verifier 1:
*    unsigned long *** e_tildei: i-th array of delta e'i, permuted bits of the errors
*    unsigned char ** d2: array of delta openings for com2
* - Auxiliary elements:
*    const int * pos: position of elements of A invertible modulo each of the factors of x^n+1
*
***********************************************************************/
int verifier_linear(const fmpz_mod_poly_t ** A,  const fmpz_mod_poly_t ** B,  const fmpz_mod_poly_t ** C1,  const fmpz_mod_poly_t ** C2,  const fmpz_mod_poly_t ** C3, const fmpz_mod_poly_t * lambda1, const fmpz_mod_poly_t * lambda2, int delta, const unsigned char ** com1, const unsigned char ** com2, mpz_t *** g1, mpz_t *** g2, mpz_t *** g3, const unsigned char ** seed1, const unsigned char ** seed2, const unsigned char ** seed3, const fmpz_mod_poly_t *** y1, const fmpz_mod_poly_t *** y2, const fmpz_mod_poly_t *** y3, const fmpz_mod_poly_t ** s1, const fmpz_mod_poly_t ** s2, const fmpz_mod_poly_t ** s3, const unsigned char ** d1, const unsigned long *** e_tilde1, const unsigned long *** e_tilde2, const unsigned long *** e_tilde3, const unsigned char ** d2, int * pos){

    // Obtain first challenge from initial commitments via a XOF
    fmpz_t * alpha = (fmpz_t *) malloc(delta*sizeof(fmpz_t));
    unsigned char alphaSeed[SEEDOL];
    const unsigned char ** aux_coms[2] = {com1, com2};
    const fmpz_mod_poly_t * statement[5*K+2];
    for (int i = 0; i < K; i++){
        statement[i] = A[i];
        statement[i+K] = B[i];
        statement[i+2*K] = C1[i];
        statement[i+3*K] = C2[i];
        statement[i+4*K] = C3[i];
    }
    statement[5*K] = lambda1;
    statement[5*K+1] = lambda2;
    alphas(alpha, alphaSeed, SEEDOL, delta, delta, statement, 5*K+2, aux_coms, 2);

    // Obtain second challenge from previous messages via a XOF
    unsigned char round = 1;
    unsigned char bSeed[SEEDOL];
    unsigned char * secondChallenge = malloc((int) ceil(delta/8.0)*sizeof(unsigned char));
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, &round, sizeof(unsigned char));
    for(int j = 0; j < 5*K+2; j++) {
        hash_poly_FFT(MD_CTX, statement[j]);
    }
    EVP_DigestUpdate(MD_CTX, alphaSeed, SEEDOL);
    for (int rep = 0; rep < delta; rep++){
        for(int j = 0; j < KAPPA + 1; j++){
            hash_vector(MD_CTX, g1[rep][j], 2*N*K);
            hash_vector(MD_CTX, g2[rep][j], 2*N*K);
            hash_vector(MD_CTX, g3[rep][j], 2*N*K);
        }
    }
    EVP_DigestFinalXOF(MD_CTX, bSeed, SEEDOL);
    EVP_MD_CTX_reset(MD_CTX);
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, bSeed, SEEDOL);
    EVP_DigestFinalXOF(MD_CTX, secondChallenge, (int) ceil(delta/8.0));
    EVP_MD_CTX_reset(MD_CTX);

    int res = 1;

    EVP_MD_CTX *mdctx_com;
    unsigned char md_value_com[md_len_com];

    fmpz_mod_poly_t ** z1, ** z2, ** z3;

    #pragma omp parallel private(mdctx_com,md_value_com,z1,z2,z3) reduction(&:res)
    {
        mdctx_com = EVP_MD_CTX_new();
        fmpz_init_array_fft(&z1,K);
        fmpz_init_array_fft(&z2,K);
        fmpz_init_array_fft(&z3,K);

        #pragma omp for
        for (int rep = 0; rep < delta; rep++){
            EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
            if (!(secondChallenge[rep/8] & (1 << rep%8))) {
                /* verifier 0 */
                res &= verifier0_aux(seed1[rep], alpha[rep], g1[rep], y1[rep], s1[rep], A, B, C1, mdctx_com, z1, pos);
                res &= verifier0_aux(seed2[rep], alpha[rep], g2[rep], y2[rep], s2[rep], A, B, C2, mdctx_com, z2, pos);
                res &= verifier0_aux(seed3[rep], alpha[rep], g3[rep], y3[rep], s3[rep], A, B, C3, mdctx_com, z3, pos);
                EVP_DigestUpdate(mdctx_com, d1[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                res &= verify_aux_com("com1", com1[rep], md_value_com, md_len_com);

                /* testing linearity */

                for (int j=0; j<D; j++) {
                    fmpz_mod_poly_mulmod_preinv(z1[pos[j]][j],lambda1[j],z1[pos[j]][j],modulus[j], modulus_inv[j], ctx);
                    fmpz_mod_poly_mulmod_preinv(z2[pos[j]][j],lambda2[j],z2[pos[j]][j],modulus[j], modulus_inv[j], ctx);
                    fmpz_mod_poly_add(z1[pos[j]][j], z1[pos[j]][j], z2[pos[j]][j],ctx);
                }

                for (int j=0; j<D; j++) {
                 res &= fmpz_mod_poly_equal(z1[pos[j]][j],z3[pos[j]][j],ctx);
                }

            } else {
                /* verifier 1 */
                res &= verifier1_aux(alpha[rep], g1[rep], e_tilde1[rep], mdctx_com);
                res &= verifier1_aux(alpha[rep], g2[rep], e_tilde2[rep], mdctx_com);
                res &= verifier1_aux(alpha[rep], g3[rep], e_tilde3[rep], mdctx_com);
                EVP_DigestUpdate(mdctx_com, d2[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                res &= verify_aux_com("com2", com2[rep], md_value_com, md_len_com);
            }
            EVP_MD_CTX_reset(mdctx_com);
        }

        EVP_MD_CTX_free(mdctx_com);
        fmpz_clear_array_fft(z1,K);
        fmpz_clear_array_fft(z2,K);
        fmpz_clear_array_fft(z3,K);
    }

    for (int rep = 0; rep < delta; rep++){
        fmpz_clear(alpha[rep]);
    }
    free(alpha);
    free(secondChallenge);

    return res;
}

/**********************************************************************
* Name:  prover_multiplicative
*
* Description:
*    given a commitment key and 3 commitment openings satisfying a multiplicative relation computes a nizkpok of such openings
*
* Arguments (read only):
*    const fmpz_mod_poly_t ** A, B: commitment key, a pair of K polynomials
*    const fmpz_mod_poly_t ** Ci: i-th commitment A * m + B * r + e, K polynomials
*    const fmpz_mod_poly_t * mi: i-th commited message, a polynomial
*    const fmpz_mod_poly_t * ri: i-th commitment randomness, a uniformly random polynomial
*    const mpz_t * ei: i-th commitment randomness error, coefficients of a vector of polynomials following a discrete gaussian bounded by B
*    int delta: repetitions to obtain soundness
*
* Populates:
* - Common elements:
*    unsigned char ** md_value_com1: array of delta auxiliary commitments com1 = Hash( {seedi} || {yi} || d1 )
*    unsigned char ** md_value_com2: array of delta auxiliary commitments com2 = Hash( mu3 || mu_x || mu_plus || d2 )
*    unsigned char ** md_value_com3: array of delta auxiliary commitments com3 = Hash( {fij permuted} || {e'ij} permuted || d3 )
*    unsigned char ** md_value_com4: array of delta auxiliary commitments com4 = Hash( mu_x + m_x || mu_plus + m_plus || d4 )
*    mpz_t *** gi: i-th array of delta vectors of coefficients of K polynomials pi_ij(fij+alpha·e'ij)
*    unsigned char ** md_value_com5: array of delta auxiliary commitments com5 = Hash( beta * mu_x + alpha * beta * mu_plus + alpha^2 * mu3 || d5 )
* - Answers for verifier 0:
*    unsigned char ** seedi: i-th array of delta seeds that determine the permutations pi_j
*    fmpz_mod_poly_t *** yi: i-th  array of delta vectors of polynomials A * mui + B * rhoi + I' * sum(2^l * f_ij)
*    fmpz_mod_poly_t ** si: i-th  array of delta polynomials rhoi + alpha * ri
*    fmpz_mod_poly_t ** t_x: array of delta polynomials mu_x + m_x
*    fmpz_mod_poly_t ** t_plus: array of delta polynomials mu_plus + m_plus
*    unsigned char ** d1, d4: array of delta openings for md_value_com1 and md_value_com4
* - Answers for both verifier 0 and 1:
*    unsigned char ** d5: array of delta openings for md_value_com5
* - Answers for verifier 1:
*    unsigned long *** e_tildei: i-th array of delta e'i, permuted bits of the errors
*    fmpz_mod_poly_t ** mu3, mu_x, mu_plus: array of delta polynomials for masking m3, m_x and m_plus
*    unsigned char ** d2, d3: array of delta openings for md_value_com2 and md_value_com3
*
***********************************************************************/
int prover_multiplicative(const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t ** C1, const fmpz_mod_poly_t ** C2,  const fmpz_mod_poly_t ** C3, const fmpz_mod_poly_t * m1, const fmpz_mod_poly_t * m2, const fmpz_mod_poly_t * m3, const fmpz_mod_poly_t * r1, const fmpz_mod_poly_t * r2, const fmpz_mod_poly_t * r3, const mpz_t * e1, const mpz_t * e2, const mpz_t * e3, int delta, unsigned char ** md_value_com1, unsigned char ** md_value_com2, unsigned char ** md_value_com3, unsigned char ** md_value_com4, mpz_t *** g1, mpz_t *** g2, mpz_t *** g3, unsigned char ** md_value_com5, unsigned char ** seed1, unsigned char ** seed2, unsigned char ** seed3, fmpz_mod_poly_t *** y1, fmpz_mod_poly_t *** y2, fmpz_mod_poly_t *** y3, fmpz_mod_poly_t ** s1, fmpz_mod_poly_t ** s2, fmpz_mod_poly_t ** s3, fmpz_mod_poly_t ** t_x, fmpz_mod_poly_t ** t_plus, unsigned char ** d1, unsigned char ** d4, unsigned char ** d5, unsigned long *** e_tilde1, unsigned long *** e_tilde2, unsigned long *** e_tilde3, fmpz_mod_poly_t ** mu3, fmpz_mod_poly_t ** mu_x, fmpz_mod_poly_t ** mu_plus, unsigned char ** d2, unsigned char ** d3) {


    EVP_MD_CTX *mdctx_com1, *mdctx_com2, *mdctx_com3, *mdctx_com4, *mdctx_com5;

    fmpz_mod_poly_t * mu1, * mu2;

    fmpz_mod_poly_t ** rho1, ** rho2, ** rho3;
    fmpz_init_array_fft(&rho1,delta);
    fmpz_init_array_fft(&rho2,delta);
    fmpz_init_array_fft(&rho3,delta);

    fmpz_mod_poly_t * m_x, * m_plus;


    unsigned long ** e_primes1 = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
    unsigned long ** e_primes2 = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
    unsigned long ** e_primes3 = (unsigned long **)  malloc (sizeof(unsigned long *) * (KAPPA + 1));
    for(int i = 0; i < KAPPA + 1; i++){
        e_primes1[i] = (unsigned long *) malloc(N*K/4);
        e_primes2[i] = (unsigned long *) malloc(N*K/4);
        e_primes3[i] = (unsigned long *) malloc(N*K/4);
    }
    expand(e1, e_primes1);
    expand(e2, e_primes2);
    expand(e3, e_primes3);

    #pragma omp parallel private(mdctx_com1,mdctx_com2,mdctx_com3,mdctx_com4,mdctx_com5,mu1,mu2,m_x,m_plus)
    {
        mdctx_com1 = EVP_MD_CTX_new();
        mdctx_com2 = EVP_MD_CTX_new();
        mdctx_com3 = EVP_MD_CTX_new();
        mdctx_com4 = EVP_MD_CTX_new();

        fmpz_init_fft(&mu1);
        fmpz_init_fft(&mu2);

        fmpz_init_fft(&m_x);
        fmpz_init_fft(&m_plus);

        #pragma omp for
        for (int rep = 0; rep < delta; rep++){
            EVP_DigestInit_ex2(mdctx_com1, MD_COM_AUX, NULL);
            EVP_DigestInit_ex2(mdctx_com2, MD_COM_AUX, NULL);
            EVP_DigestInit_ex2(mdctx_com3, MD_COM_AUX, NULL);
            EVP_DigestInit_ex2(mdctx_com4, MD_COM_AUX, NULL);

            fmpz_sampler_zq_poly_FFT(mu1);
            fmpz_sampler_zq_poly_FFT(mu2);
            fmpz_sampler_zq_poly_FFT(mu3[rep]);

            prover_initial_aux(A, B, r1, e_primes1, mdctx_com1, mdctx_com3, g1[rep], seed1[rep], e_tilde1[rep], mu1, rho1[rep], y1[rep]);
            prover_initial_aux(A, B, r2, e_primes2, mdctx_com1, mdctx_com3, g2[rep], seed2[rep], e_tilde2[rep], mu2, rho2[rep], y2[rep]);
            prover_initial_aux(A, B, r3, e_primes3, mdctx_com1, mdctx_com3, g3[rep], seed3[rep], e_tilde3[rep], mu3[rep], rho3[rep], y3[rep]);


            RAND_priv_bytes(d1[rep],LAMB);
            EVP_DigestUpdate(mdctx_com1, d1[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com1, md_value_com1[rep], NULL);
            EVP_MD_CTX_reset(mdctx_com1);


            RAND_priv_bytes(d3[rep],LAMB);
            EVP_DigestUpdate(mdctx_com3, d3[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com3, md_value_com3[rep], NULL);
            EVP_MD_CTX_reset(mdctx_com3);


            fmpz_mod_poly_mult_FFT(m_x,mu1,mu2);

            fmpz_mod_poly_mult_FFT(mu1,mu1,m2);
            fmpz_mod_poly_mult_FFT(mu2,mu2,m1);
            fmpz_mod_poly_add_FFT(m_plus, mu1, mu2);

            fmpz_sampler_zq_poly_FFT(mu_x[rep]);
            fmpz_sampler_zq_poly_FFT(mu_plus[rep]);


            hash_poly_FFT(mdctx_com2,mu3[rep]);
            hash_poly_FFT(mdctx_com2,mu_x[rep]);
            hash_poly_FFT(mdctx_com2,mu_plus[rep]);

            RAND_priv_bytes(d2[rep],LAMB);
            EVP_DigestUpdate(mdctx_com2, d2[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com2, md_value_com2[rep], NULL);
            EVP_MD_CTX_reset(mdctx_com2);


            fmpz_mod_poly_add_FFT(t_x[rep], mu_x[rep], m_x);
            fmpz_mod_poly_add_FFT(t_plus[rep], mu_plus[rep], m_plus);

            hash_poly_FFT(mdctx_com4,t_x[rep]);
            hash_poly_FFT(mdctx_com4,t_plus[rep]);


            RAND_priv_bytes(d4[rep],LAMB);
            EVP_DigestUpdate(mdctx_com4, d4[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com4, md_value_com4[rep], NULL);
            EVP_MD_CTX_reset(mdctx_com4);
        }

        fmpz_clear_fft(mu1);
        fmpz_clear_fft(mu2);
        fmpz_clear_fft(m_x);
        fmpz_clear_fft(m_plus);

        EVP_MD_CTX_free(mdctx_com1);
        EVP_MD_CTX_free(mdctx_com2);
        EVP_MD_CTX_free(mdctx_com3);
        EVP_MD_CTX_free(mdctx_com4);
    }

    for(int i = 0; i < KAPPA + 1; i++){
        free(e_primes1[i]);
        free(e_primes2[i]);
        free(e_primes3[i]);
    }
    free(e_primes1);
    free(e_primes2);
    free(e_primes3);

    // Obtain first challenge from initial commitments via a XOF
    fmpz_t * alphaBeta = (fmpz_t *) malloc(2*delta*sizeof(fmpz_t));
    unsigned char alphaBetaSeed[SEEDM];
    const unsigned char ** aux_coms[4] = {(const unsigned char **) md_value_com1, (const unsigned char **) md_value_com2, (const unsigned char **) md_value_com3, (const unsigned char **) md_value_com4};
    const fmpz_mod_poly_t * statement[5*K];
    for (int i = 0; i < K; i++){
        statement[i] = A[i];
        statement[i+K] = B[i];
        statement[i+2*K] = C1[i];
        statement[i+3*K] = C2[i];
        statement[i+4*K] = C3[i];
    }
    alphas(alphaBeta, alphaBetaSeed, SEEDM, delta, 2*delta, statement, 5*K, aux_coms, 4);
    fmpz_t * alpha = &alphaBeta[0];
    fmpz_t *  beta = &alphaBeta[delta];



    fmpz_mod_poly_t * mu_3_aux, * t_x2, * t_plus2, * com5;


    #pragma omp parallel private(mdctx_com5, mu_3_aux, t_x2, t_plus2, com5)
    {
        mdctx_com5 = EVP_MD_CTX_new();
        fmpz_init_fft(&mu_3_aux);
        fmpz_init_fft(&t_x2);
        fmpz_init_fft(&t_plus2);
        fmpz_init_fft(&com5);
        #pragma omp for
        for (int rep = 0; rep < delta; rep++){
            prover_middle_aux(alpha[rep], e_tilde1[rep], g1[rep]);
            prover_middle_aux(alpha[rep], e_tilde2[rep], g2[rep]);
            prover_middle_aux( beta[rep], e_tilde3[rep], g3[rep]);


            fmpz_mod_poly_scalar_mul_fmpz_FFT(t_x2, mu_x[rep], beta[rep]);
            fmpz_mod_poly_scalar_mul_fmpz_FFT(t_plus2, mu_plus[rep], beta[rep]);
            fmpz_mod_poly_scalar_mul_fmpz_FFT(t_plus2, t_plus2, alpha[rep]);
            fmpz_mod_poly_scalar_mul_fmpz_FFT(mu_3_aux, mu3[rep], alpha[rep]);
            fmpz_mod_poly_scalar_mul_fmpz_FFT(mu_3_aux, mu_3_aux, alpha[rep]);
            fmpz_mod_poly_add_FFT(com5, t_x2, t_plus2);
            fmpz_mod_poly_add_FFT(com5, com5, mu_3_aux);

            EVP_DigestInit_ex2(mdctx_com5, MD_COM_AUX, NULL);
            hash_poly_FFT(mdctx_com5,com5);
            RAND_priv_bytes(d5[rep],LAMB);
            EVP_DigestUpdate(mdctx_com5, d5[rep], LAMB);
            EVP_DigestFinal_ex(mdctx_com5, md_value_com5[rep], NULL);
            EVP_MD_CTX_reset(mdctx_com5);
        }
        fmpz_clear_fft(mu_3_aux);
        fmpz_clear_fft(t_x2);
        fmpz_clear_fft(t_plus2);
        fmpz_clear_fft(com5);

        EVP_MD_CTX_free(mdctx_com5);
    }



    // Obtain second challenge from previous messages via a XOF
    unsigned char round = 1;
    unsigned char bSeed[SEEDM];
    unsigned char * secondChallenge = malloc((int) ceil(delta/8.0)*sizeof(unsigned char));
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, &round, sizeof(unsigned char));
    for(int j = 0; j < 5*K; j++) {
        hash_poly_FFT(MD_CTX, statement[j]);
    }
    EVP_DigestUpdate(MD_CTX, alphaBetaSeed, SEEDM);
    for (int rep = 0; rep < delta; rep++){
        for(int j = 0; j < KAPPA + 1; j++){
            hash_vector(MD_CTX, g1[rep][j], 2*N*K);
            hash_vector(MD_CTX, g2[rep][j], 2*N*K);
            hash_vector(MD_CTX, g3[rep][j], 2*N*K);
        }
    }
    EVP_DigestFinalXOF(MD_CTX, bSeed, SEEDM);
    EVP_MD_CTX_reset(MD_CTX);
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, bSeed, SEEDM);
    EVP_DigestFinalXOF(MD_CTX, secondChallenge, (int) ceil(delta/8.0));
    EVP_MD_CTX_reset(MD_CTX);

    #pragma omp parallel for
    for (int rep = 0; rep < delta; rep++){
        if (!(secondChallenge[rep/8] & (1 << rep%8))) {
            /* prover 0 */
            prover0_aux(r1, alpha[rep], rho1[rep], s1[rep]);
            prover0_aux(r2, alpha[rep], rho2[rep], s2[rep]);
            prover0_aux(r3,  beta[rep], rho3[rep], s3[rep]);
        }
    }

    free(alphaBeta);
    free(secondChallenge);
    fmpz_clear_array_fft(rho1,delta);
    fmpz_clear_array_fft(rho2,delta);
    fmpz_clear_array_fft(rho3,delta);

  return 1;
}

/**********************************************************************
* Name:  verifier_multiplicative
*
* Description:
*    given a commitment key, 3 commitments and a nizkpok of 3 openings satisfying a multiplicative relation verifies the latter
*
* Arguments (read only, EXCEPT elements of g, t_x, t_plus, mu_plus and mu3 that are overwritten during the verification):
*    const fmpz_mod_poly_t ** A, B: commitment key, a pair of K polynomials
*    const fmpz_mod_poly_t ** Ci: i-th commitment A * m + B * r + e, K polynomials
*    int delta: repetitions to obtain soundness
* - Common elements:
*    unsigned char ** com1: array of delta auxiliary commitments com1 = Hash( {seedi} || {yi} || d1 )
*    unsigned char ** com2: array of delta auxiliary commitments com2 = Hash( mu3 || mu_x || mu_plus || d2 )
*    unsigned char ** com3: array of delta auxiliary commitments com3 = Hash( {fij permuted} || {e'ij} permuted || d3 )
*    unsigned char ** com4: array of delta auxiliary commitments com4 = Hash( mu_x + m_x || mu_plus + m_plus || d4 )
*    unsigned char ** com5: array of delta auxiliary commitments com5 = Hash( beta * mu_x + alpha * beta * mu_plus + alpha^2 * mu3 || d5 )
*    mpz_t *** gi: i-th array of delta vectors of coefficients of K polynomials pi_ij(fij+alpha·e'ij)
* - Answers for verifier 0:
*    unsigned char ** seedi: i-th array of delta seeds that determine the permutations pi_j
*    fmpz_mod_poly_t *** yi: i-th array of delta vectors of polynomials A * mui + B * rhoi + I' * sum(2^l * f_ij)
*    fmpz_mod_poly_t ** si: i-th array of delta polynomials rhoi + alpha * ri
*    fmpz_mod_poly_t ** t_x: array of delta polynomials mu_x + m_x
*    fmpz_mod_poly_t ** t_plus: array of delta polynomials mu_plus + m_plus
*    unsigned char ** d1, d4: array of delta openings for com1 and com4
* - Answers for both verifier 0 and 1:
*    unsigned char ** d5: array of delta openings for com5
* - Answers for verifier 1:
*    unsigned long *** e_tildei: i-th array of delta e'i, permuted bits of the errors
*    fmpz_mod_poly_t ** mu3, mu_x, mu_plus: array of delta polynomials for masking m3, m_x and m_plus
*    unsigned char ** d2, d3: array of delta openings for com2 and com3
* - Auxiliary elements:
*    const fmpz_mod_poly_t * a_inv: inverse of the pos-th polynomial in A
*    const int * pos: position of elements of A invertible modulo each of the factors of x^n+1
*
***********************************************************************/
int verifier_multiplicative(const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t ** C1, const fmpz_mod_poly_t ** C2, const fmpz_mod_poly_t ** C3, int delta, const unsigned  char ** com1, const unsigned  char ** com2, const unsigned  char ** com3, const unsigned  char ** com4, const unsigned  char ** com5, mpz_t *** g1, mpz_t *** g2, mpz_t *** g3, const unsigned char ** seed1, const unsigned char ** seed2, const unsigned char ** seed3, const fmpz_mod_poly_t *** y1, const fmpz_mod_poly_t *** y2, const fmpz_mod_poly_t *** y3, const fmpz_mod_poly_t ** s1, const fmpz_mod_poly_t ** s2, const fmpz_mod_poly_t ** s3, fmpz_mod_poly_t ** t_x, fmpz_mod_poly_t ** t_plus, const unsigned char ** d1,  const unsigned char ** d4, const unsigned char ** d5, const unsigned long *** e_tilde1, const unsigned long *** e_tilde2, const unsigned long *** e_tilde3, fmpz_mod_poly_t ** mu3, const fmpz_mod_poly_t ** mu_x, fmpz_mod_poly_t ** mu_plus, const unsigned char ** d2, const unsigned char ** d3, const fmpz_mod_poly_t * a_inv, const int * pos) {

    // Obtain first challenge from initial commitments via a XOF
    fmpz_t * alphaBeta = (fmpz_t *) malloc(2*delta*sizeof(fmpz_t));
    unsigned char alphaBetaSeed[SEEDM];
    const unsigned char ** aux_coms[4] = {com1, com2, com3, com4};
    const fmpz_mod_poly_t * statement[5*K];
    for (int i = 0; i < K; i++){
        statement[i] = A[i];
        statement[i+K] = B[i];
        statement[i+2*K] = C1[i];
        statement[i+3*K] = C2[i];
        statement[i+4*K] = C3[i];
    }
    alphas(alphaBeta, alphaBetaSeed, SEEDM, delta, 2*delta, statement, 5*K, aux_coms, 4);
    fmpz_t * alpha = &alphaBeta[0];
    fmpz_t *  beta = &alphaBeta[delta];

    // Obtain second challenge from previous messages via a XOF
    unsigned char round = 1;
    unsigned char bSeed[SEEDM];
    unsigned char * secondChallenge = malloc((int) ceil(delta/8.0)*sizeof(unsigned char));
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, &round, sizeof(unsigned char));
    for(int j = 0; j < 5*K; j++) {
        hash_poly_FFT(MD_CTX, statement[j]);
    }
    EVP_DigestUpdate(MD_CTX, alphaBetaSeed, SEEDM);
    for (int rep = 0; rep < delta; rep++){
        for(int j = 0; j < KAPPA + 1; j++){
            hash_vector(MD_CTX, g1[rep][j], 2*N*K);
            hash_vector(MD_CTX, g2[rep][j], 2*N*K);
            hash_vector(MD_CTX, g3[rep][j], 2*N*K);
        }
    }
    EVP_DigestFinalXOF(MD_CTX, bSeed, SEEDM);
    EVP_MD_CTX_reset(MD_CTX);
    EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
    EVP_DigestUpdate(MD_CTX, bSeed, SEEDM);
    EVP_DigestFinalXOF(MD_CTX, secondChallenge, (int) ceil(delta/8.0));
    EVP_MD_CTX_reset(MD_CTX);


    fmpz_mod_poly_t ** z1, ** z2, ** z3;
    fmpz_mod_poly_t * t_1, * t_2, * t_3;
    fmpz_mod_poly_t *c5, *c5aux;

    int res = 1;

    EVP_MD_CTX *mdctx_com;
    unsigned char md_value_com[md_len_com];


    #pragma omp parallel private(mdctx_com,md_value_com,z1,z2,z3,t_1,t_2,t_3,c5,c5aux) reduction(&:res)
    {
        mdctx_com = EVP_MD_CTX_new();

        fmpz_init_array_fft(&z1,K);
        fmpz_init_array_fft(&z2,K);
        fmpz_init_array_fft(&z3,K);
        fmpz_init_fft(&t_1);
        fmpz_init_fft(&t_2);
        fmpz_init_fft(&t_3);

        fmpz_init_fft(&c5);
        fmpz_init_fft(&c5aux);

        #pragma omp for
        for (int rep = 0; rep < delta; rep++){
            if (!(secondChallenge[rep/8] & (1 << rep%8))) {
                /* verifier 0 */
                EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
                res &= verifier0_aux(seed1[rep], alpha[rep], g1[rep], y1[rep], s1[rep], A, B, C1, mdctx_com, z1, pos);
                res &= verifier0_aux(seed2[rep], alpha[rep], g2[rep], y2[rep], s2[rep], A, B, C2, mdctx_com, z2, pos);
                res &= verifier0_aux(seed3[rep],  beta[rep], g3[rep], y3[rep], s3[rep], A, B, C3, mdctx_com, z3, pos);
                EVP_DigestUpdate(mdctx_com, d1[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                EVP_MD_CTX_reset(mdctx_com);
                res &= verify_aux_com("com1", com1[rep], md_value_com, md_len_com);

                EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
                hash_poly_FFT(mdctx_com,t_x[rep]);
                hash_poly_FFT(mdctx_com,t_plus[rep]);
                EVP_DigestUpdate(mdctx_com, d4[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                EVP_MD_CTX_reset(mdctx_com);
                res &= verify_aux_com("com4", com4[rep], md_value_com, md_len_com);

                for (int j=0; j<D; j++) {
                    fmpz_mod_poly_mulmod_preinv(t_1[j], a_inv[j], z1[pos[j]][j], modulus[j], modulus_inv[j], ctx);
                    fmpz_mod_poly_mulmod_preinv(t_2[j], a_inv[j], z2[pos[j]][j], modulus[j], modulus_inv[j], ctx);
                    fmpz_mod_poly_mulmod_preinv(t_3[j], a_inv[j], z3[pos[j]][j], modulus[j], modulus_inv[j], ctx);
                }
                fmpz_mod_poly_scalar_mul_fmpz_FFT(c5, t_x[rep], beta[rep]);
                // c5 = beta * t_x
                fmpz_mod_poly_mult_FFT(c5aux,t_1,t_2);
                fmpz_mod_poly_scalar_mul_fmpz_FFT(c5aux, c5aux, beta[rep]);
                fmpz_mod_poly_sub_FFT(c5, c5, c5aux);
                // c5 = beta * t_x - beta * t1 * t2
                fmpz_mul(beta[rep],beta[rep],alpha[rep]);
                fmpz_mod(beta[rep],beta[rep],fQ);
                // beta[rep] contains alpha * beta
                fmpz_mul(alpha[rep],alpha[rep],alpha[rep]);
                fmpz_mod(alpha[rep],alpha[rep],fQ);
                // alpha[rep] contains alpha * alpha
                fmpz_mod_poly_scalar_mul_fmpz_FFT(c5aux, t_3, alpha[rep]);
                fmpz_mod_poly_add_FFT(c5, c5, c5aux);
                // c5 = beta * t_x + alpha^2 * t3 - beta * t1 * t2
                fmpz_mod_poly_scalar_mul_fmpz_FFT(c5aux, t_plus[rep], beta[rep]);
                fmpz_mod_poly_add_FFT(c5, c5, c5aux);
                // c5 = beta * t_x + alpha * beta * t_plus + alpha^2 * t3 - beta * t1 * t2

                EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
                hash_poly_FFT(mdctx_com,c5);
                EVP_DigestUpdate(mdctx_com, d5[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                EVP_MD_CTX_reset(mdctx_com);
                res &= verify_aux_com("com5", com5[rep], md_value_com, md_len_com);

            } else {
                /* verifier 1 */
                EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
                hash_poly_FFT(mdctx_com,mu3[rep]);
                hash_poly_FFT(mdctx_com,mu_x[rep]);
                hash_poly_FFT(mdctx_com,mu_plus[rep]);

                EVP_DigestUpdate(mdctx_com, d2[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                EVP_MD_CTX_reset(mdctx_com);

                res &= verify_aux_com("com2", com2[rep], md_value_com, md_len_com);

                EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
                res &= verifier1_aux(alpha[rep], g1[rep], e_tilde1[rep], mdctx_com);
                res &= verifier1_aux(alpha[rep], g2[rep], e_tilde2[rep], mdctx_com);
                res &= verifier1_aux( beta[rep], g3[rep], e_tilde3[rep], mdctx_com);
                EVP_DigestUpdate(mdctx_com, d3[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                EVP_MD_CTX_reset(mdctx_com);

                res &= verify_aux_com("com3", com3[rep], md_value_com, md_len_com);

                fmpz_mod_poly_scalar_mul_fmpz_FFT(c5, mu_x[rep], beta[rep]);
                // c5 = beta * mu_x
                fmpz_mul(beta[rep],beta[rep],alpha[rep]);
                fmpz_mod(beta[rep],beta[rep],fQ);
                // beta[rep] contains alpha * beta
                fmpz_mul(alpha[rep],alpha[rep],alpha[rep]);
                fmpz_mod(alpha[rep],alpha[rep],fQ);
                // alpha[rep] contains alpha * alpha
                fmpz_mod_poly_scalar_mul_fmpz_FFT(mu_plus[rep], mu_plus[rep], beta[rep]);
                fmpz_mod_poly_add_FFT(c5, c5, mu_plus[rep]);
                // c5 = beta * mu_x + alpha * beta * mu_plus
                fmpz_mod_poly_scalar_mul_fmpz_FFT(mu3[rep], mu3[rep], alpha[rep]);
                fmpz_mod_poly_add_FFT(c5, c5, mu3[rep]);
                // c5 = beta * mu_x + alpha * beta * mu_plus + alpha^2 * mu3

                EVP_DigestInit_ex2(mdctx_com, MD_COM_AUX, NULL);
                hash_poly_FFT(mdctx_com,c5);
                EVP_DigestUpdate(mdctx_com, d5[rep], LAMB);
                EVP_DigestFinal_ex(mdctx_com, md_value_com, NULL);
                EVP_MD_CTX_reset(mdctx_com);

                res &= verify_aux_com("com5", com5[rep], md_value_com, md_len_com);

            }
        }
        fmpz_clear_array_fft(z1,K);
        fmpz_clear_array_fft(z2,K);
        fmpz_clear_array_fft(z3,K);

        fmpz_clear_fft(t_1);
        fmpz_clear_fft(t_2);
        fmpz_clear_fft(t_3);

        fmpz_clear_fft(c5);
        fmpz_clear_fft(c5aux);

        EVP_MD_CTX_free(mdctx_com);
    }

    free(alphaBeta);
    free(secondChallenge);

    return res;
}
