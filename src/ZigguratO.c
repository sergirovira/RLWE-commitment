#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "aux.h"
#include <mpfr.h>
#include "sampler.h"
#include "param.h"
#include <openssl/err.h>
#include <openssl/rand.h>
#include <inttypes.h>

/*A FEW AUXILIARY FUNCTIONS*/

void p_inverse(mpf_t res, mpf_t sigma, mpf_t y){

    mpf_pow_ui(res,sigma,2); //res = sigma^2
    mpf_mul_ui(res,res,2); //res = 2*sigma^2

    mpf_t aux2;
    mpf_init(aux2);
    mpf_set_d(aux2,log(mpf_get_d(y))); //aux2 = log(y)

    mpf_mul(res,res,aux2); //res = 2*sigma^2 * log(y)

    mpf_t aux3;
    mpf_init(aux3);
    mpf_set_d(aux3,-1); //aux3 = -1

    mpf_mul(res,res,aux3);

    mpf_sqrt(res,res);

    mpf_clear(aux2);
    mpf_clear(aux3);
}

void p(mpf_t res, mpf_t sigma, mpz_t x){

    mpfr_t sigma_aux;
    mpfr_t x_aux;
    mpfr_t aux;
    mpfr_t res_aux;

    mpfr_init(sigma_aux);
    mpfr_init(x_aux);
    mpfr_init(aux);
    mpfr_init(res_aux);

    mpfr_set_f(sigma_aux,sigma,MPFR_RNDU);
    mpfr_set_z(x_aux,x,MPFR_RNDU);

    mpfr_pow_ui(sigma_aux,sigma_aux,2,MPFR_RNDU);
    mpfr_pow_ui(x_aux,x_aux,2,MPFR_RNDU);
    mpfr_set_d(aux,-0.5,MPFR_RNDU);

    mpfr_div(res_aux,x_aux,sigma_aux,MPFR_RNDU);
    mpfr_mul(res_aux,res_aux,aux,MPFR_RNDU);

    mpfr_exp(res_aux,res_aux,MPFR_RNDU);

    mpfr_get_f(res,res_aux,MPFR_RNDU);

    mpfr_clear(sigma_aux);
    mpfr_clear(x_aux);
    mpfr_clear(aux);
    mpfr_clear(res_aux);

    mpfr_free_cache();

}

void S_compute(mpf_t res, int m, int c, mpf_t sigma){

    mpf_t aux;
    mpf_init(aux);
    mpf_set_d(aux,m*sqrt(PI/2));

    mpf_mul_ui(res,sigma,c);
    mpf_div(res,res,aux);

    mpf_clear(aux);
}

void print_array(double * X, int m){
    printf("[");
    for(int i = 0; i <= m; i++){
        printf("%0.20f,", X[i]);
    }
    printf("]");
}

void floorarray(mpz_t * XBar, mpf_t * X, int m){
    for(int i = 0; i <= m; i++){
        mpz_set_f(XBar[i],X[i]);
    }
}

void checksize(double * X, double * Y, int m){
    double mult;
    for(int i = 1; i <= m; i++){
        mult = (1 + X[i])*(Y[i-1] - Y[i]);
        printf("%f\n", mult);
    }
}

//Converts 4 bytes into an integer value
unsigned int encode(unsigned char buf[], int slice){
    return (buf[slice] << 24 | buf[slice + 1] << 16 | buf[slice + 2] << 8 | buf[slice + 3]);
}

/*END OF AUXILIARY FUNCTIONS*/

/*COMPUTING THE COORDINTES (X,Y) OF THE RECTANGLES FOR ZIGGURATO*/

int compute_rectangles(int m, mpf_t t, mpf_t sigma, mpf_t * XBar, mpz_t * X, mpf_t * Y)
{
    mpf_t X_aux[m+1];
    mpf_p_init_array(X_aux,m+1);

    mpf_t t_aux;
    mpf_init(t_aux);


    mpf_add_ui(t_aux,t,1); //t_aux = t + 1

    int c;

    mpf_t S;
    mpf_init(S);


    mpf_set_ui(Y[0],0);

    int correct_partition = 0;
    int correct_semipartition = 0;
    int it = 0;

    mpf_mul(XBar[m],t,sigma); //X[m] = t*sigma
    mpf_mul(X_aux[m],t,sigma);

    mpf_t upperbound;
    mpf_init(upperbound);

    mpf_mul(upperbound,t_aux,sigma); //upperbound = (t+1)*sigma


    while(!correct_partition && mpf_cmp(upperbound,XBar[m])){
        it++;

        c = 0;
        while(-1*mpf_cmp_ui(Y[0],1) > 0 || !correct_semipartition){

            c++;

            correct_semipartition = 1;

            S_compute(S,m,c,sigma);

            mpf_set_ui(Y[m],0);
            mpf_set_ui(XBar[0],0);

            mpf_floor(X_aux[m],XBar[m]);
            mpf_add_ui(X_aux[m],X_aux[m],1);

            mpf_div(Y[m-1],S,X_aux[m]);

            if(mpf_cmp_ui(Y[m-1],1) > 0){
                correct_semipartition = 0;
                break;
            } else {
                p_inverse(XBar[m-1],sigma,Y[m-1]);
            }

            for(int i = m-2; i >= 1; i--){

                mpf_floor(X_aux[i+1],XBar[i+1]);
                mpf_add_ui(X_aux[i+1],X_aux[i+1],1);
                mpf_div(Y[i],S,X_aux[i+1]);
                mpf_add(Y[i],Y[i],Y[i+1]);
                if(mpf_cmp_ui(Y[i],1) > 0){
                    correct_semipartition = 0;
                } else {
                    p_inverse(XBar[i],sigma,Y[i]);
                }
            }

            mpf_floor(X_aux[1],XBar[1]);
            mpf_add_ui(X_aux[1],X_aux[1],1);
            mpf_div(Y[0],S,X_aux[1]);
            mpf_add(Y[0],Y[0],Y[1]);
        }

        if(!isnan(mpf_get_d(Y[0]))){
            correct_partition = 1;
        } else {
            mpf_add_ui(XBar[m],XBar[m],1);
            mpf_set_ui(Y[0],0);
            c = 0;
        }
    }

    mpf_mul(t_aux,t_aux,sigma);

    if(mpf_cmp(t_aux,XBar[m]) == 0){
        printf("\nDZ: NO VALID PARTITION FOUND \n");
        return 0;
    } else {
        floorarray(X,XBar,m);
    }


    mpf_clear(t_aux);
    mpf_p_free_array(X_aux,m+1);
    mpf_clear(upperbound);
    mpf_clear(S);

    return 1;

}
/*OPTIMIZATION*/

void sLine(mpf_t res, int i, mpz_t x1, mpz_t x2, mpf_t y1, mpf_t y2, mpz_t x){

    mpf_t div;
    mpf_t yhat1;
    mpf_t yhat2;

    mpf_t x1_aux;
    mpf_t x2_aux;
    mpf_t x_aux;

    mpf_init(div);
    mpf_init(yhat1);
    mpf_init(yhat2);

    mpf_init(x1_aux);
    mpf_init(x2_aux);
    mpf_init(x_aux);

    mpf_set_z(x1_aux,x1);
    mpf_set_z(x2_aux,x2);
    mpf_set_z(x_aux,x);

    if(mpz_cmp(x1,x2) == 0){
        mpf_set_si(res,-1);
    } else {
        if(i > 1){
            mpf_set(yhat1,y1);
        } else {
            mpf_set_ui(yhat1,1);
        }
        mpf_set(yhat2,y2);

        mpf_sub(x1_aux,x2_aux,x1_aux);
        mpf_sub(res,yhat2,yhat1);

        mpf_div(res,res,x1_aux);

        mpf_sub(x2_aux,x_aux,x2_aux);

        mpf_mul(res,res,x2_aux);
    }

    mpf_clear(x1_aux);
    mpf_clear(x2_aux);
    mpf_clear(x_aux);
    mpf_clear(div);
    mpf_clear(yhat1);
    mpf_clear(yhat2);

}

/*DISCRETE GAUSSIAN SAMPLING*/
int ZigguratO(mpz_t res, int m, mpf_t sigma, mpz_t * XBar, mpf_t * YBar, int w){
    int i;
    int s;
    mpz_t x;
    mpz_t x_aux;
    int y;
    mpf_t ybar;
    mpz_t Xbar_aux;
    int S[2] = {-1,1};

    unsigned char buf[1]; //buffer to store 1 random byte
    int response;
    unsigned long err;
    unsigned int rand_1, rand_2, rand_3, rand_4; //the 4 random integer values extracted from buf

    mpz_init(x);
    mpz_init(x_aux);
    mpz_init(Xbar_aux);
    mpf_init(ybar);


    while(1) {

        response = RAND_priv_bytes(buf, sizeof(buf)); //generate 16 random bytes and store them into buff
        err = ERR_get_error();

        if(response!=1){
            printf("RNG failed, error code: %lu", err);
            exit(0);
        }

        rand_1 = randUnsignedInt(0, m-1);
        rand_2 = ((buf[0] & (1 << 0)) != 0); //get the last bit of buf
        rand_3 = ((buf[0] & (1 << 1)) != 0);
        rand_4 = randUnsignedInt(0, (int)(pow(2,w))-1);

        i = rand_1 + 1;  // 2: Choose rectangle, sign and value
        s = S[rand_2];

        mpz_set(Xbar_aux,XBar[i]);
        mpz_add_ui(Xbar_aux,Xbar_aux,1);


        sampler_mpz(x, Xbar_aux);

        mpz_set(x_aux,x);
        mpz_sub_ui(x_aux,x_aux,1);

        if(mpz_cmp_ui(x,0) > 0 && mpz_cmp(x,XBar[i-1]) >= 0){ //3: If x is in a left rectangle return sx
            mpz_mul_si(res,x,s);
            mpz_clear(x);
            mpz_clear(x_aux);
            mpz_clear(Xbar_aux);
            mpf_clear(ybar);
            return 1;
        }


        if(mpz_cmp_ui(x,0) == 0){ // 5: If x = 0 return sx with probability 1/2
            int b = rand_3;
            if(b == 0){
                mpz_mul_si(res,x,s);
                mpz_clear(x);
                mpz_clear(x_aux);
                mpz_clear(Xbar_aux);
                mpf_clear(ybar);
                return 1;
            }
        }

        y = rand_4; // 10: Rejection area of R_i

        mpf_sub(ybar,YBar[i-1],YBar[i]);
        mpf_mul_ui(ybar,ybar,y);

        mpf_t down;
        mpf_t up;
        mpf_t value;
        mpf_init(down);
        mpf_init(up);
        mpf_init(value);

        sLine(down,i,XBar[i-1],XBar[i],YBar[i-1],YBar[i],x);
        sLine(up,i,XBar[i-1],XBar[i],YBar[i-1],YBar[i],x_aux);

        mpf_mul_ui(down,down,pow(2,w));
        mpf_mul_ui(up,up,pow(2,w));

        p(value,sigma,x);

        mpf_sub(value,value,YBar[i]);
        mpf_mul_ui(value,value,pow(2,w));


        mpf_t xbar_aux;
        mpf_init(xbar_aux);

        mpf_set_z(xbar_aux, XBar[i]);
        mpf_add_ui(xbar_aux,xbar_aux,1);


        if(mpf_cmp(xbar_aux,sigma) <= 0){  //11: Concave-down case
            if(mpf_cmp(ybar,down) <= 0 || mpf_cmp(ybar,value) <= 0){
                mpz_mul_si(res,x,s);

                mpf_clear(xbar_aux);
                mpf_clear(down);
                mpf_clear(up);
                mpf_clear(value);
                mpz_clear(x);
                mpz_clear(x_aux);
                mpz_clear(Xbar_aux);
                mpf_clear(ybar);

                return 1;
            }
        }

        mpf_set_z(xbar_aux,XBar[i-1]);

        if(mpf_cmp(sigma,xbar_aux) <= 0){ // 14: Concave-up case
            if(mpf_cmp(ybar,up) < 0 && mpf_cmp(ybar,value) <= 0){
                mpz_mul_si(res,x,s);

                mpf_clear(xbar_aux);
                mpf_clear(down);
                mpf_clear(up);
                mpf_clear(value);
                mpz_clear(x);
                mpz_clear(x_aux);
                mpz_clear(Xbar_aux);
                mpf_clear(ybar);


                return 1;
            }
        }

        if(mpf_cmp(ybar,value) <= 0){
            mpz_mul_si(res,x,s);

            mpf_clear(xbar_aux);
            mpf_clear(down);
            mpf_clear(up);
            mpf_clear(value);
            mpz_clear(x);
            mpz_clear(x_aux);
            mpz_clear(Xbar_aux);
            mpf_clear(ybar);

            return 1;
        }

        mpf_clear(xbar_aux);
        mpf_clear(down);
        mpf_clear(up);
        mpf_clear(value);
    }
}
