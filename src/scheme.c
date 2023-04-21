#include <stdio.h>
#include "aux.h"
#include "sampler.h"
#include "param.h"
#include "ZigguratO.h"


/**********************************************************************
* Name:        key_gen
*
* Description: Generates the public key (a,b)
*
* Arguments:   fmpz_mod_poly_t * a: vector of polynomials in Rq of size K
*              fmpz_mod_poly_t * b: vector of polynomials in Rq of size K
***********************************************************************/

void key_gen(fmpz_mod_poly_t ** a, fmpz_mod_poly_t ** b, fmpz_mod_poly_t * a_inv, int * pos){

    for(int i = 0; i < K; i++){
        fmpz_sampler_zq_poly_FFT(a[i]);
        fmpz_sampler_zq_poly_FFT(b[i]);
    }

    if (a_inv != NULL) {
        for(int j = 0; j<D;j++) {
            pos[j] = 0;
            int invertible = 0;
            while(!invertible){
                invertible = fmpz_mod_poly_invmod(a_inv[j], a[pos[j]][j], modulus[j], ctx);
                if (!invertible) {
                    pos[j] += 1;
                }
                if (pos[j] == K) {
                    #ifdef DEBUG
                        red();
                        printf("!!! NOT INVERTIBLE A\n");
                        reset();
                    #endif
                    fmpz_clear_fft(a_inv);
                    key_gen(a,b,a_inv,pos);
                    invertible = 1;
                }
            }
        }
    }
}

/**********************************************************************
* Name:        commit
*
* Description: Produces a commitment c = (am + br + e) and an opening d = (m,r,e).
*
* Arguments:
*           fmpz_mod_poly_t * m: a polynomial in Rq corresponding to the message.
*           fmpz_mod_poly_t * a: vector a of the public key.
*          fmpz_mod_poly_t * b: vector b of the public key.
*           fmpz_mod_poly_t * c: the commitment.
*           mpz_t * e: error term.
*          fmpz_mod_poly_t * r: a random polynomial in Rq
***********************************************************************/

void commit(fmpz_mod_poly_t * m, fmpz_mod_poly_t ** a, fmpz_mod_poly_t ** b, fmpz_mod_poly_t ** c, mpz_t * e, fmpz_mod_poly_t * r) {

    fmpz_mod_poly_t ** mult1;
    fmpz_init_array_fft(&mult1,K);
    fmpz_mod_poly_t ** mult2;
    fmpz_init_array_fft(&mult2,K);


    fmpz_sampler_zq_poly_FFT(r);

    for(int i = 0; i < K; i++){
        fmpz_mod_poly_mult_FFT(mult1[i], m, a[i]);
        fmpz_mod_poly_mult_FFT(mult2[i], r, b[i]);

        fmpz_mod_poly_add_FFT(c[i], mult1[i], mult2[i]);
    }

    mpz_t aux1;
    mpz_init_set(aux1,BOUND);

    mpz_t aux2;
    mpz_init(aux2);
    mpz_add_ui(aux2,BOUND,1);

    mpz_t res;
    mpz_init(res);

    mpz_t res_aux;
    mpz_init(res_aux);

    mpf_t t;
    mpf_init_set_ui(t,TAILCUT);

    mpf_mul(t,t,SIGMA);
    mpz_set_f(res_aux,t); //floor(t*sigma)

    do {
        for(int i = 0; i < N*K; i++){
            ZigguratO(res, NUMBEROFRECTANGLES, SIGMA, XBar, YBar, W);
            if(mpz_cmp(res,res_aux) != 0){ //if res = floor(t*sigma) resample
                mpz_init_set(e[i],res);
            } else {
                i--;
            }
            
        }
        mpz_p_norm(e, N*K, aux2);
    } while(mpz_cmp(aux2,aux1) > 0 || mpz_cmp(aux2,aux1) == 0);

    fmpz_mod_poly_t ** e_split;
    fmpz_init_array_fft(&e_split,K);

    for(int i = 0; i < K; i++){
        for(int j = 0; j < N; j++){
            fmpz_mod_poly_set_coeff_mpz(e_split[i][j/(N/D)], j%(N/D), e[i*N+j], ctx);
        }
        fmpz_p_fcr(e_split[i],0,D,ALPHA,OMEGA,D,fQ);
        fmpz_mod_poly_add_FFT(c[i], c[i], e_split[i]);
    }

    fmpz_clear_array_fft(mult1,K);
    fmpz_clear_array_fft(mult2,K);
    fmpz_clear_array_fft(e_split,K);

    mpz_clear(aux1);
    mpz_clear(aux2);
    mpz_clear(res);
}

/**********************************************************************
* Name:        verify
*
* Description: from a commitment c and an opening d = (m,r,e)
*               accepts if c == am + br + e.
*
* Arguments:
*           fmpz_mod_poly_t * c: the commitment.
*           fmpz_mod_poly_t * m: a polynomial in Rq corresponding to the message.
*          fmpz_mod_poly_t * r: a random polynomial in Rq
*           mpz_t * e: error term.
*           fmpz_mod_poly_t * a: vector a of the public key.
*          fmpz_mod_poly_t * b: vector b of the public key.
*
* Returns 1 if c == am + br + e, 0 otherwise.
***********************************************************************/

int verify(fmpz_mod_poly_t ** c, fmpz_mod_poly_t * m, fmpz_mod_poly_t * r, mpz_t * e, fmpz_mod_poly_t ** a, fmpz_mod_poly_t ** b) {

    mpz_t aux;
    mpz_init_set_si(aux,pow(2,KAPPA));

    int result;

    mpz_t norm;
    mpz_init(norm);

    mpz_p_norm(e, N*K, norm);

    if(mpz_cmp(norm,aux) > 0 || mpz_cmp(norm,aux) == 0){
        #ifdef DEBUG
            red();
            printf("||e|| >= 2^%zu\n", KAPPA);
            reset();
        #endif
        result = 0;
    } else {
        #ifdef DEBUG
            green();
            printf("||e|| < 2^%zu\n", KAPPA);
            reset();
        #endif
        result = 1;
    }

    int sum = 1;

    fmpz_mod_poly_t ** e_split;
    fmpz_init_array_fft(&e_split,K);

    for(int i = 0; i < K; i++){
        for(int j = 0; j < N; j++){
            fmpz_mod_poly_set_coeff_mpz(e_split[i][j/(N/D)], j%(N/D), e[i*N+j], ctx);
        }
        fmpz_p_fcr(e_split[i],0,D,ALPHA,OMEGA,D,fQ);
    }


    fmpz_mod_poly_t ** mult1;
    fmpz_init_array_fft(&mult1,K);
    fmpz_mod_poly_t ** mult2;
    fmpz_init_array_fft(&mult2,K);
    fmpz_mod_poly_t ** res;
    fmpz_init_array_fft(&res,K);

    for(int i = 0; i < K; i++) {
        fmpz_mod_poly_mult_FFT(mult1[i],m,a[i]);
        fmpz_mod_poly_mult_FFT(mult2[i],r,b[i]);
        fmpz_mod_poly_add_FFT(res[i],mult1[i],mult2[i]);
        fmpz_mod_poly_add_FFT(res[i],res[i],e_split[i]);
    }

    for(int i = 0; i < K; i++) {
        for(int j=0; j<D; j++) {
            sum &= fmpz_mod_poly_equal(c[i][j],res[i][j],ctx);
        }
    }



    fmpz_clear_array_fft(mult1,K);
    fmpz_clear_array_fft(mult2,K);

    fmpz_clear_array_fft(res,K);
    fmpz_clear_array_fft(e_split,K);

    mpz_clear(aux);
    mpz_clear(norm);

    if(sum != 0){
        #ifdef DEBUG
            green();
            printf("Correct form of the commitment\n");
            reset();
        #endif
    } else {
        #ifdef DEBUG
            red();
            printf("Incorrect form of the commitment\n");
            reset();
        #endif
    }

    result &= sum;

    return result;
}
