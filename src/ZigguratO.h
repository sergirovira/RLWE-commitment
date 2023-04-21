extern int XBAR[];
extern double YBAR[];

double get_sample(int m, double t, double sigma, int w);
int ZigguratO(mpz_t res, int m, mpf_t sigma, mpz_t * XBar, mpf_t * YBar, int w);
int compute_rectangles(int m, mpf_t t, mpf_t sigma, mpf_t * XBar, mpz_t * X, mpf_t * Y);
