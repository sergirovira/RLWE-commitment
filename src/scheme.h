void key_gen(fmpz_mod_poly_t ** a, fmpz_mod_poly_t ** b, fmpz_mod_poly_t * a_inv, int * pos);
void commit(fmpz_mod_poly_t * m, fmpz_mod_poly_t ** a, fmpz_mod_poly_t ** b, fmpz_mod_poly_t ** c, mpz_t * e, fmpz_mod_poly_t * r);
int verify(fmpz_mod_poly_t ** c, fmpz_mod_poly_t * m, fmpz_mod_poly_t * r, mpz_t * e, fmpz_mod_poly_t ** a, fmpz_mod_poly_t ** b);
