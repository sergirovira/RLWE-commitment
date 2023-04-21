#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/evp.h>
#include <ctype.h>
#include <gmp.h>
#include "ds_benchmark.h"
#include "time.h"
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod_poly.h"
#include "flint/fmpz_mod.h"
#include "sampler.h"
#include "permutations.h"
#include "param.h"

/* PRINTING FUNCTIONS */

void red () {
    printf("\033[1;31m");
}

void green() {
    printf("\033[0;32m");
}

void reset () {
    printf("\033[0m");
}

/*************************************************
* Name:        print_hex
*
* Description: prints in hexadecimal an array of unsigned chars
*
* Arguments:   char name[]: name of the array
*              unsigned char * a: array of characters to print
*              int length: number of characters
**************************************************/
void print_hex(char name[], const unsigned char * a, int length) {
    printf("%s: ",name);
    for(int i = 0; i < length; i++) {
        printf("%02x", a[i]);
    }
    printf("\n");
}

/*************************************************
* Name:        test_message
*
* Description: returns 1 if both auxiliary commitments are equal, 0 otherwise
*
* Arguments:   char name[]: name of the verification
*              int result: non-zero if everithing ok, 0 otherwise
**************************************************/
void test_message(char name[], int result) {
    if (result) {
        green();
        printf("[PASS] All %s tests passed\n\n",name);
        reset();
    } else {
        red();
        printf("[FAIL] Some %s tests failed\n\n",name);
        reset();
    }
}

void print_fmpz(char name[], fmpz_t a) {
    printf("%s = ",name);
    fmpz_print(a);
    printf("\n");
}

void mpz_p_print(mpz_t a[], int l, char * p){
    printf("%s = [", p);
    for(int i = 0; i < l-1; i++){
        gmp_printf ("%Zd,", a[i]);
    }
    gmp_printf ("%Zd]",a[l-1]);
    printf("\n");
}

void mpz_p_tofile(mpz_t a[], int l, char * f){
    FILE * fp;
    /* open the file for writing*/
    fp = fopen (f,"a");

    for(int i = 0; i < l-1; i++){
        mpz_out_str(fp, 10, a[i]);
        fprintf(fp,",");
    }
    mpz_out_str(fp, 10, a[l-1]);
    gmp_fprintf (fp,"\n");

    /* close the file*/
    fclose (fp);
}

void mpf_p_tofile(mpf_t a[], int l, char * f){
   FILE * fp;
   /* open the file for writing*/
   fp = fopen (f,"a");

   for(int i = 0; i < l-1; i++){
        mpf_out_str(fp, 10, 0, a[i]);
        fprintf(fp,",");
    }
    mpf_out_str(fp, 10, 0, a[l-1]);
    gmp_fprintf (fp,"\n");

    /* close the file*/
    fclose (fp);
}

void tofile_double_2(double p, double c, char * f){
    FILE * fp;
    /* open the file for writing*/
    fp = fopen (f,"a");
    fprintf(fp,"%g, %g \n",p,c);

    /* close the file*/
    fclose (fp);
}

/* INIT AND CLEAR FUNCTIONS */

void fmpz_init_fft(fmpz_mod_poly_t ** a) {
    *a = NULL;
    *a = (fmpz_mod_poly_t *) malloc(D*sizeof(fmpz_mod_poly_t));
    for (int j=0; j<D; j++) {
        fmpz_mod_poly_init((*a)[j],ctx);
    }
}

void fmpz_clear_fft(fmpz_mod_poly_t * a) {
    for (int j=0; j<D; j++ ) {
        fmpz_mod_poly_clear(a[j],ctx);
    }
    free(a);
}

void fmpz_init_array_fft(fmpz_mod_poly_t *** a, unsigned int k) {
    *a = (fmpz_mod_poly_t **) malloc(k*sizeof(fmpz_mod_poly_t*));
    for (int i=0; i<k; i++) {
        (*a)[i] = (fmpz_mod_poly_t *) malloc(D*sizeof(fmpz_mod_poly_t));
        for (int j=0; j<D; j++ ) {
            fmpz_mod_poly_init((*a)[i][j],ctx);
        }
    }
}

void fmpz_clear_array_fft(fmpz_mod_poly_t ** a, unsigned int k) {
    for (int i=0; i<k; i++) {
        for (int j=0; j<D; j++ ) {
            fmpz_mod_poly_clear(a[i][j],ctx);
        }
        free(a[i]);
    }
    free(a);
}

void mpz_p_init_array(mpz_t a[], int l){
    for (int i = 0; i < l; i++) {
        mpz_init(a[i]);
    }
}

void mpf_p_init_array(mpf_t a[], int l){
    for (int i = 0; i < l; i++) {
        mpf_init(a[i]);
    }
}


void mpz_p_free_array(mpz_t a[], int l){
    for (int i = 0; i < l; i++) {
            mpz_clear(a[i]);
    }
}

void mpf_p_free_array(mpf_t a[], int l){
    for (int i = 0; i < l; i++) {
            mpf_clear(a[i]);
    }
}


void mpz_p_free_matrix(mpz_t ** a, int l, int t){
    for (int i = 0; i < l; i++) {
        mpz_p_free_array(a[i],t);
        free(a[i]);
    }
    free(a);
}


/* BIT FUNCTIONS */

// index / (8 * sizeof(unsigned long) ) is the index of the long that contains the i-th bit
// index % (8 * sizeof(unsigned long) ) if the index of the i-th bit in the binary decomposition of the long that contains it
// << is the left shift operator, so that 0001 << 2 becomes 0100
// ~ is the bitwise not operator, so that ~0001 becomes 1110

int get_bit(const unsigned long * v, int index){
    // returns the i-th bit from the array v
    return (v[index / (8 * sizeof(unsigned long) )] & ((unsigned long) 1 << (index % (8 * sizeof(unsigned long) )))) ? 1 : 0;
}

void set_bit(unsigned long * v, int index, int bit){
    // sets the i-th bit at the array v
    /*
    printf("we want the %d-th bit to be %d\n",index,bit);
    printf("that corresponds to the %d-th bit from the %d-th unsigned long integer\n",index % (8 * sizeof(unsigned long) ),index / (8 * sizeof(unsigned long) ));
    unsigned long aux[1] = {((unsigned long) 1 << (index % (8 * sizeof(unsigned long) )))};
    printf("we have to operate with %lu \n",aux[0]);
    printf("that has binary decomposition: ");
    print_bit_array(aux, 8 * sizeof(unsigned long));
     */
    if (bit) {
        v[index / (8 * sizeof(unsigned long) )] |= ((unsigned long) 1 << (index % (8 * sizeof(unsigned long) )));
    } else {
        v[index / (8 * sizeof(unsigned long) )] &= ~((unsigned long) 1 << (index % (8 * sizeof(unsigned long) )));
    }
}

void expand(const mpz_t * e, unsigned long ** e_primes){
    /* Compute extensions e_prime_j, i.e, compute e_prime_j = (e_bar_j || w) such that e_prime_j
    has the same number of 0's and 1's */

    for(int j = 0; j < KAPPA + 1; j++){
        int s = 0;
        for (int i = 0; i < N*K; i++) {
            if ((mpz_tstbit(e[i],j) & (j < KAPPA)) || ((mpz_sgn(e[i]) >= 0) & (j == KAPPA))) {
                set_bit(e_primes[j], i, 1);
                s++;
            } else {
                set_bit(e_primes[j], i, 0);
            }
        }
        for(int z = 0; z < s; z++){
            set_bit(e_primes[j], N*K + z, 0);
        }
        for(int z = 0; z < N*K - s; z++){
            set_bit(e_primes[j], N*K + s + z, 1);
        }
    }
}

/* ARRAY FUNCTIONS */

void fmpz_mod_poly_mult_FFT(fmpz_mod_poly_t * c, const fmpz_mod_poly_t * a, const fmpz_mod_poly_t * b){
    for (int j = 0; j < D; j++){
        fmpz_mod_poly_mulmod_preinv(c[j], a[j], b[j], modulus[j], modulus_inv[j], ctx);
    }
}

void fmpz_mod_poly_add_FFT(fmpz_mod_poly_t * c, const fmpz_mod_poly_t * a, const fmpz_mod_poly_t * b){
    for (int j = 0; j < D; j++){
        fmpz_mod_poly_add(c[j], a[j], b[j], ctx);
    }
}

void fmpz_mod_poly_sub_FFT(fmpz_mod_poly_t * c, const fmpz_mod_poly_t * a, const fmpz_mod_poly_t * b){
    for (int j = 0; j < D; j++){
        fmpz_mod_poly_sub(c[j], a[j], b[j], ctx);
    }
}

void fmpz_mod_poly_scalar_mul_fmpz_FFT(fmpz_mod_poly_t * b, fmpz_mod_poly_t * a, const fmpz_t alpha){
    for (int j=0; j<D; j++) {
        fmpz_mod_poly_scalar_mul_fmpz(b[j], a[j], alpha, ctx);
    }
}

void mpz_p_mod(mpz_t a[], int l, mpz_t b[]){
    for(int i = 0; i < l; i++){
        mpz_mod(b[i],a[i],Q);
    }
}

int max(int x, int y) {
    return x ^ ((x ^ y) & -(x < y));
}

int min(int x, int y) {
    return y ^ ((x ^ y) & -(x < y));
}

void mpz_p_add(mpz_t a[], mpz_t b[], int l1, int l2, mpz_t c[]){

    int m1 = max(l1,l2);
    int m2 = min(l1,l2);

    for(int i = 0; i < m2; i++){
        mpz_add(c[i],a[i],b[i]);
    }

    if(m1 == l1){
        for(int i = m2; i < m1; i++){
            mpz_set(c[i],a[i]);
        }
    } else {
        for(int i = m2; i < m1; i++){
            mpz_set(c[i],b[i]);
        }
    }
    mpz_p_mod(c,m1,c);
}

/* OTHER AUXILIARY FUNCTIONS */

int vector_check(const unsigned long * vec){
        int sum = 0;

        for(int i = 0; i < 2*N*K; i++){
            sum += get_bit(vec,i);
        }

        if(sum == N*K){
            return 1;
        } else {
            return 0;
        }
}

/*************************************************
* Name:        lattice_membership
*
* Description: checks if a vector of polynomials a belongs to a lattice defined by A
*
* Arguments:   mpz_t ** a: vector of polynomials
*              mpz_t ** A: vector of polynomials defining the lattice
*
* Returns 1 if a belongs to the lattice defined by A or 0 otherwise
**************************************************/
int lattice_membership(fmpz_mod_poly_t ** a, const fmpz_mod_poly_t ** A, int * pos) {
    fmpz_mod_poly_t test1, test2;
    fmpz_mod_poly_init(test1,ctx);
    fmpz_mod_poly_init(test2,ctx);


    int test = 1;

    for(int j = 0; j < D; j++){
        for(int i = 0; i < K; i++){
            if (i != pos[j]) {
                fmpz_mod_poly_mulmod_preinv(test1, A[pos[j]][j], a[i][j], modulus[j], modulus_inv[j], ctx);
                fmpz_mod_poly_mulmod_preinv(test2, A[i][j], a[pos[j]][j], modulus[j], modulus_inv[j], ctx);
                test &= fmpz_mod_poly_equal(test1, test2, ctx);
            }
        }
    }

    fmpz_mod_poly_clear(test1,ctx);
    fmpz_mod_poly_clear(test2,ctx);

    return test;
}

/*************************************************
* Name:        hash_poly_FFT
*
* Description: hash D fmpz_mod_poly
*
* Arguments:   fmpz_mod_poly * f: polynomial
*              EVP_MD_CTX * mdctx_com: digest context
*
*
* Returns 1 if the hash is computed without errors
**************************************************/
int hash_poly_FFT(EVP_MD_CTX * mdctx_com, const fmpz_mod_poly_t * f) {
    char * ibuf = (char *) malloc((mpz_sizeinbase(Q,2)+8-1)/8);
    size_t size;
    mpz_t coeff;
    mpz_init(coeff);
    for(int j = 0; j < N; j++) {
        fmpz_mod_poly_get_coeff_mpz(coeff, f[j/(N/D)], j%(N/D), ctx);
        mpz_export(ibuf,&size,1,1,0,0,coeff);
        if(!EVP_DigestUpdate(mdctx_com, ibuf, size)){
            #ifdef DEBUG
                printf("Error computing hash");
            #endif
            exit(0);
        }
    }
    free(ibuf);
    mpz_clear(coeff);
    return 1;
}

/*************************************************
* Name:        hash_vector
*
* Description: hash a vector of mpz_t
*
* Arguments:   mpz_t * v: a vector of mpz_t
*              int * N: size of vector v
*              EVP_MD_CTX * mdctx_com: digest context
*
*
* Returns 1 if the hash is computed without errors
**************************************************/
int hash_vector(EVP_MD_CTX * mdctx_com, mpz_t * v, int N) {
    size_t size = (mpz_sizeinbase(Q,2)+8-1)/8;
    char * aux_buf = (char *) calloc(N*size,sizeof(char));
    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        mpz_mod(v[i],v[i],Q);
        mpz_export(aux_buf+i*size,NULL,1,1,0,0,v[i]);
    }
    if(!EVP_DigestUpdate(mdctx_com, aux_buf, N*size)){
        #ifdef DEBUG
            printf("Error computing hash");
        #endif
        exit(0);
    }
    free(aux_buf);
    return 1;
}


/*************************************************
* Name:        verify_aux_com
*
* Description: returns 1 if both auxiliary commitments are equal, 0 otherwise
*
* Arguments:   char name[]: name of the commitment
*              unsigned char * com1: first auxiliary commitment
*              unsigned char * com2: second auxiliary commitment
*              int length: length of the commitments
**************************************************/
int verify_aux_com(char name[], const unsigned char * com1, const unsigned char * com2, int length) {
    int res = 1;
    for (int i = 0; i < length; i++) {
        res &= com1[i] == com2[i];
    }
    #ifdef DEBUG
        if(!res){
                red();
                printf("Error in the verification of %s\n",name);
                reset();
        } else {
                green();
                printf("Commitment %s is correct\n",name);
                reset();
        }
    #endif
    return res;
}




void root_of_unity(fmpz_t omega, fmpz_t q, int d){
    // omega = (1)^{1/n}
    int test = 0;
    fmpz_t u;
    fmpz_init(u);
    fmpz_t aux;
    fmpz_init(aux);
    fmpz_t exp;
    fmpz_init(exp);
    fmpz_sub_ui(exp, q, 1);
    fmpz_fdiv_q_2exp(exp,exp,1);
    while (!test) {
        fmpz_sampler_zq(u);
        fmpz_powm(aux,u,exp,q);
        fmpz_add_ui(aux, aux, 1);
        test = fmpz_divisible(aux,q);
    }
    fmpz_sub_ui(exp, q, 1);
    fmpz_fdiv_q_2exp(exp,exp,round(log2(d)));
    fmpz_powm(omega,u,exp,q);
    fmpz_clear(u);
    fmpz_clear(aux);
    fmpz_clear(exp);
}

int bitSize(int a) {
    int bits = 1;
    while (a >>= 1) ++ bits;
    return bits;
}

void mpz_p_norm(mpz_t a[], int l, mpz_t maximum){
    mpz_abs(maximum,a[0]);
    for (int i = 1; i < l; i++) {
        if (mpz_cmpabs(a[i],maximum) > 0) {
            mpz_abs(maximum,a[i]);
        }
    }
}

int fmpz_p_fcr(fmpz_mod_poly_t * a, int i, int j, const fmpz_t alpha, const fmpz_t omega, int r, const fmpz_t m){

    if(r == 1){
        return 1;
    }

    int k = (int)floor((i + j) / 2);

    fmpz_t exp;
    fmpz_init(exp);
    fmpz_powm_ui(exp, alpha, (int)floor(r/2), m);

    for (int z = k; z < j; z++) {
        fmpz_mod_poly_scalar_mul_fmpz(a[z],a[z],exp,ctx);
    }

    for (int z = i; z < k; z++) {
        fmpz_mod_poly_add(a[z],a[z],a[k + z - i],ctx);
    }

    for (int z = k; z < j; z++) {
        fmpz_mod_poly_add(a[z],a[z],a[z],ctx);
        fmpz_mod_poly_neg(a[z],a[z],ctx);
    }

    for (int z = k; z < j; z++) {
        fmpz_mod_poly_add(a[z],a[z],a[i + z - k],ctx);
    }

    fmpz_t exp_omega;
    fmpz_init(exp_omega);
    fmpz_powm_ui(exp_omega, omega, 2, m);

    fmpz_t alpha_omega;
    fmpz_init(alpha_omega);
    fmpz_mod_mul(alpha_omega,alpha,omega,ctx);

    fmpz_p_fcr(a,i,k,alpha,exp_omega,(int)floor(r/2),m);
    fmpz_p_fcr(a,k,j,alpha_omega,exp_omega,(int)floor(r/2),m);

    fmpz_clear(exp);
    fmpz_clear(exp_omega);
    fmpz_clear(alpha_omega);


    return 1;
}



int reverseBits(int i, int d){
    int j = 0;
    for(int z = 0; z < d; z++){
        j *= 2;
        j += i%2;
        i = (int)floor(i/2);
    }
    return j;
}
