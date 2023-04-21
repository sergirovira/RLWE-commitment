#include <stdint.h>
#include <gmp.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod_poly.h"
#include "flint/fmpz_mod.h"
#include <openssl/evp.h>

void red();
void reset();
void green();
void print_hex(char name[], const unsigned char * a, int length);

void print_fmpz(char name[], fmpz_t a);
int verify_aux_com(char name[], const unsigned char * com1, const unsigned char * com2, int length);
void test_message(char name[], int result);
int lattice_membership(fmpz_mod_poly_t ** a, const fmpz_mod_poly_t ** A, const int * pos);
int vector_check(const unsigned long * vec);
int get_bit(const unsigned long * v, int index);
void set_bit(unsigned long v[], int index, int bit);
void expand(const mpz_t * e, unsigned long ** e_primes);
int hash_poly_FFT(EVP_MD_CTX * mdctx_com, const fmpz_mod_poly_t * f);
int hash_vector(EVP_MD_CTX * mdctx_com, mpz_t * v, int N);
void root_of_unity(fmpz_t omega, fmpz_t q, int d);

int bitSize(int a);
int reverseBits(int i, int d);


void fmpz_init_fft(fmpz_mod_poly_t ** a);
void fmpz_clear_fft(fmpz_mod_poly_t * a);
void fmpz_init_array_fft(fmpz_mod_poly_t *** a, unsigned int k);
void fmpz_clear_array_fft(fmpz_mod_poly_t ** a, unsigned int k);
void fmpz_to_fft(fmpz_mod_poly_t a,fmpz_mod_poly_t * aa);
void fft_to_fmpz(fmpz_mod_poly_t * aa,fmpz_mod_poly_t a);
void fmpz_array_to_fft(fmpz_mod_poly_t * a,fmpz_mod_poly_t ** aa, unsigned int k);
void fft_array_to_fmpz(fmpz_mod_poly_t ** aa,fmpz_mod_poly_t * a, unsigned int k);
void fmpz_mod_poly_mult(fmpz_mod_poly_t a, fmpz_mod_poly_t b, int l, int n, fmpz_mod_poly_t c);
void fmpz_mod_poly_mult_FFT(fmpz_mod_poly_t * c, const fmpz_mod_poly_t * a, const fmpz_mod_poly_t * b);
void fmpz_mod_poly_add_FFT(fmpz_mod_poly_t * c, const fmpz_mod_poly_t * a, const fmpz_mod_poly_t * b);
void fmpz_mod_poly_sub_FFT(fmpz_mod_poly_t * c, const fmpz_mod_poly_t * a, const fmpz_mod_poly_t * b);
void fmpz_mod_poly_scalar_mul_fmpz_FFT(fmpz_mod_poly_t * b, const fmpz_mod_poly_t * a, const fmpz_t alpha);

int fmpz_p_fcr(fmpz_mod_poly_t * a, int i, int j, const fmpz_t alpha, const fmpz_t omega, int r, const fmpz_t m);

void mpz_p_init_array(mpz_t a[], int l);
void mpz_p_add(mpz_t a[], mpz_t b[], int l1, int l2, mpz_t c[]);

void mpz_p_tofile(mpz_t a[], int l, char * f);
void mpf_p_tofile(mpf_t a[], int l, char * f);

void mpz_p_free_array(mpz_t a[], int l);
void mpz_p_free_matrix(mpz_t ** a, int l, int t);
void mpf_p_init_array(mpf_t a[], int l);
void mpf_p_free_array(mpf_t a[], int l);

void mpz_p_norm(mpz_t a[], int l, mpz_t maximum);

void tofile_double_2(double p, double c, char * f);
