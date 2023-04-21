/* libraries used by the permutation algorithm */
#include <openssl/evp.h> // to get a XOF
#include <openssl/rand.h> // to sample the seeds
#include <gmp.h> // to work with large integers

/* to handle bits */
#include "aux.h"

/* to use cryptographically secure samplers */
#include "sampler.h"

/* to get the SHA3 contexts */
#include "param.h"


EVP_MD_CTX *MD_CTX;
const EVP_MD *MD_SHAKE;

double expansionFactor = 4.5;


/*
 * In order to permute an array of size N the Fisher-Yates algorithm samples indexes uniformly at random from [i,N-1] for every 0≤i<N-1.
 * Equivalently we only need to sample from [0,N-i-1] and add i.
 * To get a uniform distribution in [0,N-i-1] we sample an integer from [0,2^⌈log2(N-i)⌉], repeating if it is greater or equal than N-i.
 * The integer in [0,⌈2^log2(N-i)⌉] is casted from ⌈log256(N-i)⌉ unsigned chars pseudorandomly obtained from a seed using SHAKE128 or SHAKE256.
 */

/* permuting and inverting a permutation with a CSPRNG */
void permute(mpz_t v[],int N,unsigned char * seed) {
    /*
     * Shuffles an array with Fisher-Yates using an array of lamb unsigned chars (128 bits if lamb is 16) as source of randomness for a XOF function.
     * @param v Vector of mpz_t to be shuffled.
     * @param N Size of the array v as an int.
     * @param seed Seed for the XOF as an array of lamb unsigned chars.
     */

    size_t bufferLength = bitSize(N)*N/expansionFactor; // this should be optimized
    unsigned char * buffer = (unsigned char *) malloc(bufferLength * sizeof(unsigned char));
    size_t bufferPosition = 0;
    unsigned char counter = 0;
    refreshBuffer(buffer,bufferLength,seed,&counter);

    // Shuffle the array
    for(int i = 0; i < N-1; i++){
        // Sample an integer index j uniformly from the interval [i,N-1]
        unsigned int j = randUnsignedIntB(i, N-1, buffer,bufferLength, &bufferPosition,seed, &counter);
        // Swap the i-th and j-th elements
        mpz_swap(v[i],v[j]);
    }

    free(buffer);
}

void permuteInv(mpz_t v[],int N,unsigned char * seed) {
    /*
     * Shuffles an array inverting the permutation produced with Fisher-Yates using an array of lamb unsigned chars (128 bits if lamb is 16) as source of randomness for a XOF function.
     * @param v Vector of mpz_t to be shuffled.
     * @param N Size of the array v as an int.
     * @param seed Seed for the XOF as an array of lamb unsigned chars.
     */

    size_t bufferLength = bitSize(N)*N/expansionFactor; // this should be optimized
    unsigned char * buffer = (unsigned char *) malloc(bufferLength * sizeof(unsigned char));
    size_t bufferPosition = 0;
    unsigned char counter = 0;
    refreshBuffer(buffer,bufferLength,seed,&counter);

    // Precompute all the swaps of the Fisher-Yates shuffle
    unsigned int js[N];
    for(int i = 0; i < N-1; i++){
        // Sample an integer index j uniformly from the interval [i,N-1]
        unsigned int j = randUnsignedIntB(i, N-1, buffer,bufferLength, &bufferPosition,seed, &counter);
        // Swap the i-th and j-th elements
        js[i] = j;
    }
    // Shuffle the array in reversed order
    for(int i = N-2; i >=0; i--){
        mpz_swap(v[i],v[js[i]]);
    }

    free(buffer);
}
/* END of permuting and inverting a permutation with a CSPRNG from the seed */

/* permuting and inverting a permutation of bits with a CSPRNG */
void permuteBit(unsigned long v[],int N,unsigned char * seed) {
    /*
     * Shuffles an array of N bits, encoded as unsigned longs,
     * with Fisher-Yates using an array of lamb unsigned chars (128 bits if lamb is 16) as source of randomness for a XOF function.
     * @param v Vector of bits encoded as unsigned longs to be shuffled.
     * @param N Size of the array v as an int.
     * @param seed Seed for the XOF as an array of lamb unsigned chars.
     *
     * Requires: ctx Message Digest context from OpenSSL implementation of the SHA3 family.
     * Requires: md_shake Message Digest function from OpenSSL implementation of the SHA3 family, either SHAKE128 or SHAKE256.
     */

    size_t bufferLength = bitSize(N)*N/expansionFactor; // this should be optimized
    unsigned char * buffer = (unsigned char *) malloc(bufferLength * sizeof(unsigned char));
    size_t bufferPosition = 0;
    unsigned char counter = 0;
    refreshBuffer(buffer,bufferLength,seed,&counter);

    int bitAux;

    // Shuffle the array
    for(int i = 0; i < N-1; i++){
        // Sample an integer index j uniformly from the interval [i,N-1]
        unsigned int j = randUnsignedIntB(i, N-1, buffer,bufferLength, &bufferPosition,seed, &counter);
        // Swap the i-th and j-th elements
        bitAux = get_bit(v,i);
        set_bit(v,i,get_bit(v,j));
        set_bit(v,j,bitAux);
    }

    free(buffer);
}

void permuteInvBit(unsigned long v[],int N,unsigned char * seed) {
    /*
     * Shuffles an array of N bits, encoded as unsigned longs,
     * inverting the permutation produced with Fisher-Yates using an array of lamb unsigned chars (128 bits if lamb is 16) as source of randomness for a XOF function.
     * @param v Vector of bits encoded as unsigned longs to be shuffled.
     * @param N Size of the array v as an int.
     * @param seed Seed for the XOF as an array of lamb unsigned chars.
     *
     * Requires: ctx Message Digest context from OpenSSL implementation of the SHA3 family.
     * Requires: md_shake Message Digest function from OpenSSL implementation of the SHA3 family, either SHAKE128 or SHAKE256.
     */

    size_t bufferLength = bitSize(N)*N/expansionFactor; // this should be optimized
    unsigned char * buffer = (unsigned char *) malloc(bufferLength * sizeof(unsigned char));
    size_t bufferPosition = 0;
    unsigned char counter = 0;
    refreshBuffer(buffer,bufferLength,seed,&counter);

    // Precompute all the swaps of the Fisher-Yates shuffle
    unsigned int js[N];
    for(int i = 0; i < N-1; i++){
        // Sample an integer index j uniformly from the interval [i,N-1]
        unsigned int j = randUnsignedIntB(i, N-1, buffer,bufferLength, &bufferPosition,seed, &counter);
        // Swap the i-th and j-th elements
        js[i] = j;
    }
    // Shuffle the array in reversed order
    int bitAux;
    for(int i = N-2; i >=0; i--){
        bitAux = get_bit(v,i);
        set_bit(v,i,get_bit(v,js[i]));
        set_bit(v,js[i],bitAux);
    }

    free(buffer);
}
/* END of permuting and inverting a permutation of bits with a CSPRNG */
