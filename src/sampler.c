#include <stdint.h>
#include <stdio.h>
#include <openssl/rand.h>
#include <openssl/evp.h>
#include <sys/time.h>
#include <time.h>
#include <gmp.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod_poly.h"
#include "flint/fmpz_mod.h"
#include "sampler.h"
#include "param.h"
#include "aux.h"
#include <openssl/err.h>

/* to cast chars to integers */
#include "string.h"


void sampler_zq(mpz_t random){
    unsigned char buffer[BQ];
    do {
        RAND_priv_bytes(buffer, BQ);
        buffer[0] &= auxSampler;
        mpz_import(random, BQ, 1, sizeof(buffer[0]), 1, 0, buffer); // most significant word first, most significant byte first
    }
    while(mpz_cmp(random,Q)>=0);
}

void sampler_zq_list(mpz_t * random, size_t l){
    double expansion = 1.2;
    size_t bSize = (size_t) l*expansion;
    size_t counter = 0;
    unsigned char * buffer;
    if( BQ != 0 && bSize > (((SIZE_MAX<PTRDIFF_MAX) ? SIZE_MAX : PTRDIFF_MAX)-1) / BQ ) {
        // the size of the buffer would overflow
        buffer = (unsigned char *) malloc((SIZE_MAX<PTRDIFF_MAX) ? SIZE_MAX : PTRDIFF_MAX);
    } else {
        buffer = (unsigned char *) malloc(BQ*bSize);
    }

    RAND_priv_bytes(buffer, BQ*bSize);

    for(int i = 0; i < l; i++){
        do {
            if (counter >= bSize) {
                counter = 0;
                RAND_priv_bytes(buffer, BQ*bSize);
            }
            buffer[BQ*counter] &= auxSampler;
            mpz_import(random[i],BQ,1,sizeof(buffer[0]),1,0,buffer+BQ*counter);
            counter++;
        } while (mpz_cmp(random[i],Q) >= 0);
    }
    free(buffer);
}

/* deterministically sample an fmpz_t ussing pseudorandomnes obtained from seed and counter and stored in buffer at position bufferPosition */
void fmpz_sampler_zq_B(fmpz_t random, unsigned char * buffer,size_t bufferLength, size_t *bufferPosition,unsigned char * seed, int SEEDS, unsigned char *counter){
    mpz_t aux;
    mpz_init(aux);
    do {
        if ((*bufferPosition) + BQ >= bufferLength) {
            refreshBuffer(buffer,bufferLength,seed,SEEDS,counter);
            (*bufferPosition) = 0;
        }
        buffer[*bufferPosition] &= auxSampler;
        mpz_import(aux, BQ, 1, sizeof(buffer[0]), 1, 0, buffer+(*bufferPosition)); // most significant word first, most significant byte first
        (*bufferPosition) += BQ;
    }
    while(mpz_cmp(aux,Q)>=0);
    fmpz_set_mpz(random, aux);

    mpz_clear(aux);
}

void fmpz_sampler_zq_poly_FFT(fmpz_mod_poly_t * random){
    mpz_t aux;
    mpz_init(aux);

    double expansion = 1.2;
    size_t bSize = (size_t) N*expansion;
    size_t counter = 0;
    unsigned char buffer[BQ*bSize];
    RAND_priv_bytes(buffer, BQ*bSize);

    for(int i = 0; i < N; i++){
        do {
            if (counter >= bSize) {
                counter = 0;
                RAND_priv_bytes(buffer, BQ*bSize);
            }
            buffer[BQ*counter] &= auxSampler;
            mpz_import(aux,BQ,1,sizeof(buffer[0]),1,0,buffer+BQ*counter);
            counter++;
        } while (mpz_cmp(aux,Q) >= 0);
        fmpz_mod_poly_set_coeff_mpz(random[i/(N/D)], i%(N/D), aux, ctx);
    }
    mpz_clear(aux);
}

void fmpz_sampler_zq(fmpz_t random){
    mpz_t aux;
    mpz_init(aux);
    sampler_zq(aux);
    fmpz_set_mpz(random, aux);

    mpz_clear(aux);
}

void samplerAuxiliaryParameters(int *bQ, int *BQ, unsigned char *auxSampler){
    mpz_sub_ui(Q,Q,1);
    *BQ = mpz_sizeinbase(Q,256);
    *bQ = mpz_sizeinbase(Q,2);
    mpz_add_ui(Q,Q,1);

    *auxSampler = 1;
    for (int i=0; i<(*bQ%8);i++){
        *auxSampler*=2;
    }
    *auxSampler-=1;
    if (*auxSampler==0) {*auxSampler-=1;}
}


void refreshBuffer(unsigned char * buffer, size_t bufferLength, unsigned char * seed, int SEEDS, unsigned char *counter) {
    // Initialize the XOF function
    #pragma omp critical
    {
        EVP_DigestInit_ex2(MD_CTX, MD_SHAKE, NULL);
        // Use the seed to get the randomness and store it in a buffer of unsigned chars
        EVP_DigestUpdate(MD_CTX, seed, SEEDS);
        EVP_DigestUpdate(MD_CTX, counter, 1);
        EVP_DigestFinalXOF(MD_CTX, buffer, bufferLength);
    }
    (*counter)++;
    if (*counter == 0) {
        red();
        printf("ERROR: counter overflow\n");
        reset();
        exit(0);
    }
}


unsigned int randUnsignedIntB(int a, int b, unsigned char * buffer,size_t bufferLength, size_t *bufferPosition,unsigned char * seed, unsigned char *counter) {
    /* returns an integer uniformly at random from the interval [a,b] */
    unsigned int j;
    int sizeb = bitSize(b-a);
    int size = (sizeb-1)/8+1;
    unsigned int reject = (1<<(8*size)) - ((1<<(8*size))%(b-a+1));
    do
    {
        if ((*bufferPosition) + size >= bufferLength)
        {
            refreshBuffer(buffer,bufferLength,seed,LAMB,counter);
            (*bufferPosition) = 0;
        }
        j = 0;

        memcpy(&j, buffer+(*bufferPosition), size);

        (*bufferPosition) += size;
    } while (j >= reject);
    j %= (b-a+1);
    j += a;
    return j;
}

unsigned int randUnsignedInt(int a, int b) {
    /* returns an integer uniformly at random from the interval [a,b] */
    unsigned int j;
    int sizeb = bitSize(b-a);
    int size = (sizeb-1)/8+1;
    unsigned char buffer[sizeb];

    unsigned int aux = 1;
    for (int e=sizeb; e<(sizeof(int)*8);e++) {
        aux <<= 1;
        aux ^= 1; //aux++;
    }
    aux = ~(aux << sizeb);

    int response;
    unsigned long err;

    do
    {
        response = RAND_priv_bytes(buffer,sizeb);
        err = ERR_get_error();

        if(response!=1){
            printf("RNG failed, error code: %lu", err);
            exit(0);
        }

        j = 0;

        memcpy(&j, buffer, size);

        j &= aux;
        j += a;
    } while (j > b);
    return j;
}

void sampler_mpz(mpz_t random, mpz_t module){
    mpz_sub_ui(module,module,1);
    int bytes = mpz_sizeinbase(module,256);
    int bits = mpz_sizeinbase(module,2);
    mpz_add_ui(module,module,1);

    int response;
    unsigned long err;

    int auxSampler = 1;
    for (int i=0; i<(bits%8);i++){
        auxSampler*=2;
    }
    auxSampler-=1;
    if (auxSampler==0) {auxSampler-=1;}
    unsigned char buffer[bytes];
    do {
        response = RAND_priv_bytes(buffer, bytes);
        err = ERR_get_error();

        if(response!=1){
            printf("RNG failed, error code: %lu", err);
            exit(0);
        }

        buffer[0] &= auxSampler;
        mpz_import(random, bytes, 1, sizeof(buffer[0]), 1, 0, buffer); // most significant word first, most significant byte first
    }
    while(mpz_cmp(random,module)>=0);
}
