int prover_initial_aux(const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t * r, unsigned long ** e_primes, EVP_MD_CTX * mdctx_com1, EVP_MD_CTX * mdctx_com2, mpz_t ** f_permuted, unsigned char * seed, unsigned long ** e_tilde, fmpz_mod_poly_t * mu, fmpz_mod_poly_t * rho, fmpz_mod_poly_t ** y);
int prover_middle_aux(const fmpz_t alpha, unsigned long ** e_tilde, mpz_t ** f_permuted);
int prover0_aux(const fmpz_mod_poly_t * r, fmpz_t alpha, fmpz_mod_poly_t * rho, fmpz_mod_poly_t * rho_r_scaled);
int verifier0_aux(const unsigned char seed[], const fmpz_t alpha, mpz_t * gs[KAPPA + 1], const fmpz_mod_poly_t ** y, const fmpz_mod_poly_t * s, const fmpz_mod_poly_t ** A, const fmpz_mod_poly_t ** B, const fmpz_mod_poly_t ** C, EVP_MD_CTX *mdctx_com, fmpz_mod_poly_t ** gs_inv, const int * pos);
int verifier1_aux(const fmpz_t alpha, mpz_t * gs[KAPPA + 1], const unsigned long * e_tilde[KAPPA + 1], EVP_MD_CTX * mdctx_com);
