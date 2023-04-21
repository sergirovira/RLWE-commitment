#ifndef PARAM_H
#define PARAM_H

#include <gmp.h>
#include <openssl/evp.h>
#include <openssl/rand.h>

//Q,N,K and KAPPA correspond to the values with the same names in the paper

// SHA3 family contexts
extern EVP_MD_CTX * MD_CTX;
extern const EVP_MD * MD_SHAKE;

extern const EVP_MD *MD_COM_AUX;
extern int md_len_com;

// the security bytes
extern int LAMB;
extern int SEEDOL;
extern int SEEDM;

// Parameter filename
#define FILEPARAMETERS "generatedParams.csv"


// Params for the sampler
extern int BQ;
extern int bQ;
extern unsigned char auxSampler;

extern mpz_t Q;
extern mpz_t BOUND;

extern fmpz_t fQ;
extern fmpz_t ALPHA;
extern fmpz_t OMEGA;
extern fmpz_t IALPHA;
extern fmpz_t IOMEGA;
extern fmpz_t IDOS;

extern fmpz_mod_ctx_t ctx;
extern fmpz_mod_poly_t * modulus;
extern fmpz_mod_poly_t * modulus_inv;

extern size_t N;
extern size_t K;
extern size_t KAPPA;
extern int D;


#define MESSAGEDIGEST "SHA256"
#define XOF "SHAKE128"

extern mpf_t SIGMA;

// NUMBEROFRECTANGLES, TAILCUT and W are DZ parameters, they should not be changed.

#define NUMBEROFRECTANGLES 10
#define TAILCUT 10
#define W 10

extern mpz_t XBar[NUMBEROFRECTANGLES + 1];
extern mpf_t YBar[NUMBEROFRECTANGLES + 1];


// ITERATIONS_SCHEME determines the number of iterations that will be run for the protocols.
// ITERATIONS_SAMPLING determines the number of samples that will be drawned from both uniform and gaussian distributions.
// ITERATIONS_OPERATIONS determines the number of polynomial multiplications to be performed.

#define ITERATIONS_SCHEME 1
#define ITERATIONS_SAMPLING 1000000
#define ITERATIONS_OPERATIONS 100

//Vervose determines if the output of the execution should print all the details of the different operations performed (0 no, 1 yes)

#undef VERBOSE
#define DEBUG
#undef DEBUG

//ALL THE PARAMETERS BELOW THIS LINE SHOULD NOT BE CHANGED

#define FILECOMMITMENT "commitment"
#define FILEVERIFICATION "verifier"
#define FILEKEYGEN "keygen"

#define FILEPROVEROPENING "proveropening"
#define FILEVERIFIEROPENING "verifieropening"
#define FILEPROVERLINEAR "proverlinear"
#define FILEVERIFIERLINEAR "verifierlinear"
#define FILEPROVERMULTIPLICATIVE "provermultiplicative"
#define FILEVERIFIERMULTIPLICATIVE "verifiermultiplicative"

// #define FILEZIGGURAT "benchmarks/ziggurat.txt"

#define PI 3.14159265358979323846

#endif
