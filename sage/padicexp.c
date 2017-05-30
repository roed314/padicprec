#include <math.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include "flint/flint.h"
#include "flint/padic.h"

unsigned long primes[10] = { 2, 3, 5, 7, 11, 13, 17, 23, 29, 31 };
unsigned long precs[10] = { 20, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 1000000 };

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


/* p-adic exponential */
void padicexp(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo) {
    double tme;

    unsigned long i, N, saveN, Np, tmp, trunc, step;
    mpz_t f, arg, trunc_mod, h, hpow, mpz_tmp, d, inv;
    mpz_t denominator;
    mpz_t *num, *denom;

    mpz_init(mpz_tmp);
    mpz_init(arg);
    mpz_set_ui(ans, 1);
    mpz_init(denominator);
    mpz_set_ui(denominator, 1);

    if (p == 2) {
        mpz_fdiv_r_ui(mpz_tmp, a, 4);
    } else {
        mpz_fdiv_r_ui(mpz_tmp, a, p);
    }
    if (mpz_cmp_ui(mpz_tmp, 0) != 0) {
       puts("Argh...");
       return;
    }
    mpz_set(arg,a);

    /* Where do we need to truncate the Taylor expansion */
    if (p == 2) {
        N = prec;
    } else {
        N = (prec*(p-1)) / (p-2);
    }
    saveN = N;

    /* We allocate memory and initialize variables */
    mpz_init(f);
    mpz_init(h); mpz_init(hpow);
    mpz_init(d); mpz_init(inv);
    num = (mpz_t*)malloc((N+1)*sizeof(mpz_t));
    denom = (mpz_t*)malloc((N+1)*sizeof(mpz_t));
    for (i = 0; i <= N; i++) {
        mpz_init(num[i]);
        mpz_init(denom[i]);
    }

    if (p == 2) {
        trunc = 4;
        mpz_init_set_ui(trunc_mod, p);
        mpz_mul_ui(trunc_mod, trunc_mod, p);
        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
    } else {
        trunc = 2;
        mpz_init_set_ui(trunc_mod, p);
        mpz_mul_ui(trunc_mod, trunc_mod, p);
    }
    while(1) {
        mpz_fdiv_r(f, arg, trunc_mod);
        mpz_sub(arg, arg, f);

        if (mpz_cmp_ui(f, 0) != 0) {

            /* We compute the Taylor expansion of exp(f)
               For now, computations are carried out over the rationals */
            mpz_set_ui(num[0], 1);
            mpz_set_ui(denom[0], 1);
            for (i = 1; i <= N; i++) {
                mpz_set_ui(num[i], 1);
                mpz_set_ui(denom[i], i);
            }
            step = 1;
            mpz_set(h, f);
            mpz_set(hpow, h);

            while(1) {
                for (i = 0; i <= N - step; i += step << 1) {
                    mpz_mul(mpz_tmp, hpow, num[i+step]);
                    mpz_mul(num[i], num[i], denom[i+step]);
                    mpz_add(num[i], num[i], mpz_tmp);
                    mpz_mul(denom[i], denom[i], denom[i+step]);
                }
                step <<= 1;
                if (step > N) break;
                mpz_mul(hpow, hpow, hpow);
            }

            /* We simplify the fraction */
            Np = N; tmp = 0;
            while(Np > 0) { Np /= p; tmp += Np; }
            mpz_ui_pow_ui(d, p, tmp);
            mpz_divexact(num[0], num[0], d);
            mpz_divexact(denom[0], denom[0], d);

            /* We add this contribution to exp(f) */
            mpz_mul(ans, ans, num[0]);
            mpz_fdiv_r(ans, ans, modulo);
            mpz_mul(denominator, denominator, denom[0]);
            mpz_fdiv_r(denominator, denominator, modulo);

        }

        if (trunc > prec) break;

        /* We update the variables for the next step */
        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        trunc <<= 1;
        N >>= 1;
    }

    /* We coerce the result from Q to Zp */
    mpz_gcdext(d, inv, NULL, denominator, modulo);
    mpz_mul(ans, ans, inv);
    mpz_fdiv_r(ans, ans, modulo);

    /* We clear memory */
    mpz_clear(arg);
    mpz_clear(denominator);
    mpz_clear(f);
    mpz_clear(trunc_mod);
    mpz_clear(h);
    mpz_clear(hpow);
    mpz_clear(mpz_tmp);
    mpz_clear(d);
    mpz_clear(inv);
    for (i = 0; i <= saveN; i++) {
        mpz_clear(num[i]);
        mpz_clear(denom[i]);
    }
    free(num);
    free(denom);
}


void padiclog_inverse(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, unsigned long precinit, mpz_t modulo) {
    unsigned long i, N, saveN, Np, tmp, trunc, step;
    mpz_t f, arg, li, bi, trunc_mod, h, hpow, mpz_tmp;
    mpz_t d, inv;
    mpz_t *num, *denom;

    N = 1 + prec;
    while(1) {
      tmp = 1 + prec + (unsigned long)(log(N)/log(p));
      if (tmp == N) break;
      N = tmp;
    }
    saveN = N;

    mpz_init(mpz_tmp);
    mpz_init(arg);
    mpz_set(arg, ans);
    mpz_set_ui(ans, 1);
    mpz_init(li);

    // We first compute l_1 = log(ans) but
    // stop the computation at some finite level
    // We also update ans so that ans = exp(l_1)

    trunc = 2;
    mpz_init_set_ui(trunc_mod, p);
    mpz_mul_ui(trunc_mod, trunc_mod, p);
    mpz_init(f);
    mpz_init(h); mpz_init(hpow);
    mpz_init(d); mpz_init(inv);
    num = (mpz_t*)malloc(N*sizeof(mpz_t));
    denom = (mpz_t*)malloc(N*sizeof(mpz_t));
    for (i = 0; i < N; i++) {
        mpz_init(num[i]);
        mpz_init(denom[i]);
    }

    while(1) {
        mpz_fdiv_r(f, arg, trunc_mod);

        if (mpz_cmp_ui(f, 1) != 0) {

            mpz_mul(ans, ans, f);
            mpz_fdiv_r(ans, ans, modulo);
            mpz_ui_sub(f, 2, f);
            mpz_mul(arg, arg, f);

            for (i = 0; i < N; i++) {
                mpz_set_ui(num[i], 1);
                mpz_set_ui(denom[i], i+1);
            }
            step = 1;
            mpz_ui_sub(h, 1, f);
            mpz_set(hpow, h);
            while(1) {
                for (i = 0; i < N - step; i += step << 1) {
                    mpz_mul(mpz_tmp, hpow, num[i+step]);
                    mpz_mul(mpz_tmp, mpz_tmp, denom[i]);
                    mpz_mul(num[i], num[i], denom[i+step]);
                    mpz_add(num[i], num[i], mpz_tmp);
                    mpz_mul(denom[i], denom[i], denom[i+step]);
                }
                step <<= 1;
                if (step >= N) break;
                mpz_mul(hpow, hpow, hpow);
            }

            Np = N; tmp = 0;
            while(Np > 0) { Np /= p; tmp += Np; }
            mpz_ui_pow_ui(d, p, tmp);
            mpz_divexact(mpz_tmp, num[0], d);
            mpz_mul(mpz_tmp, h, mpz_tmp);

            mpz_divexact(denom[0], denom[0], d);
            mpz_gcdext(d, inv, NULL, denom[0], modulo);
            mpz_mul(mpz_tmp, mpz_tmp, inv);

            mpz_add(li, li, mpz_tmp);

        }

        if (trunc > precinit) break;

        mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        trunc <<= 1;
        N >>= 1;
    }

    mpz_gcdext(d, inv, NULL, ans, modulo);
    mpz_mul(ans, ans, inv);

    // Here comes the Newton iteration
    //   a_(i+1) = a_i * (1 + b_i)
    //   l_(i+1) = l_i + log(1 + b_i)
    //   b_(i+1) = (a - l_(i+1)) correctly truncated
    // NB: the value of a_i is stored in the variable ans

    // Initialization of b_1
    N = 1 + prec/precinit;
    while(1) {
        tmp = 1 + prec/precinit + (unsigned long)(log(N)/log(p));
        if (tmp == N) break;
        N = tmp;
    }
    if (p == 2) trunc = 2*precinit - 1;
    else trunc = precinit << 1;

    mpz_ui_pow_ui(trunc_mod, p, trunc);
    mpz_init(bi);
    mpz_sub(bi, a, li);
    mpz_fdiv_r(bi, bi, trunc_mod);

    while(1) {
        if (mpz_cmp_ui(bi, 0) != 0) {

            // We set a_(i+1) = a_i * (1 + b_i)
            mpz_add_ui(mpz_tmp, bi, 1);
            mpz_mul(ans, ans, mpz_tmp);
            mpz_fdiv_r(ans, ans, modulo);

            // We compute l_(i+1) = l_i + log(1 + b_i)
            for (i = 0; i < N; i++) {
                mpz_set_ui(num[i], 1);
                mpz_set_ui(denom[i], i+1);
            }
            step = 1;
            mpz_neg(h, bi);
            mpz_set(hpow, h);
            while(1) {
                for (i = 0; i < N - step; i += step << 1) {
                    mpz_mul(mpz_tmp, hpow, num[i+step]);
                    mpz_mul(mpz_tmp, mpz_tmp, denom[i]);
                    mpz_mul(num[i], num[i], denom[i+step]);
                    mpz_add(num[i], num[i], mpz_tmp);
                    mpz_mul(denom[i], denom[i], denom[i+step]);
                }
                step <<= 1;
                if (step >= N) break;
                mpz_mul(hpow, hpow, hpow);
            }

            Np = N; tmp = 0;
            while(Np > 0) { Np /= p; tmp += Np; }
            mpz_ui_pow_ui(d, p, tmp);
            mpz_divexact(mpz_tmp, num[0], d);
            mpz_mul(mpz_tmp, h, mpz_tmp);

            mpz_divexact(denom[0], denom[0], d);
            mpz_gcdext(d, inv, NULL, denom[0], modulo);
            mpz_mul(mpz_tmp, mpz_tmp, inv);

            mpz_sub(li, li, mpz_tmp);
        }

        if (trunc > prec) break;

        // We update N, trunc and trunc_mod
        if (p == 2) {
            N = 1 + prec/trunc;
            while(1) {
                tmp = 1 + prec/trunc + (unsigned long)(log(N)/log(p));
                if (tmp == N) break;
                N = tmp;
            }
            trunc = 2*trunc - 1;
            mpz_mul(trunc_mod, trunc_mod, trunc_mod);
            mpz_divexact_ui(trunc_mod, trunc_mod, p);
        } else {
            trunc <<= 1;
            N >>= 1;
            mpz_mul(trunc_mod, trunc_mod, trunc_mod);
        }

        // We set b_(i+1) = a - l_(i+1)  mod  trunc_mod
        mpz_sub(bi, a, li);
        mpz_fdiv_r(bi, bi, trunc_mod);

    }

    mpz_fdiv_r(ans, ans, modulo);

    /* We clear memory */
    mpz_clear(arg);
    mpz_clear(f);
    mpz_clear(trunc_mod);
    mpz_clear(h);
    mpz_clear(hpow);
    mpz_clear(mpz_tmp);
    mpz_clear(d);
    mpz_clear(inv);
    mpz_clear(li);
    mpz_clear(bi);
    for (i = 0; i < saveN; i++) {
        mpz_clear(num[i]);
        mpz_clear(denom[i]);
    }
    free(num);
    free(denom);
}



int main(){
  unsigned long p;
  unsigned long prec;
  int i, j, count;

  double tme_tmp, tme, tme_inverse, tme_flint;
  mpz_t modulo, a, ans, ans_inverse;
  gmp_randstate_t state;

  fmpz_t p_flint;
  padic_t a_flint, ans_flint;
  padic_ctx_t ctx;

  mpz_init(modulo);
  mpz_init(a);
  fmpz_init(p_flint);

  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));
  mpz_init(ans);

  for (i = 0; i < 10; i++) { p = primes[i];
  for (j = 0; j < 10; j++) { prec = precs[j];

  printf("   p = %u\n", p);
  printf("prec = %u\n", prec);
  mpz_ui_pow_ui(modulo, p, prec);
  fmpz_set_ui(p_flint, p);
  padic_ctx_init(ctx, p_flint, 0, 0, PADIC_TERSE);

  tme = tme_flint = 0;
  count = 0;
  while (tme < 60 & count < 100000) {
    mpz_urandomm(a, state, modulo);
    if (p == 2) {
      mpz_mul_ui(a, a, p*p);
    } else {
      mpz_mul_ui(a, a, p);
    }

    tme_tmp = get_wall_time();
    mpz_init(ans);
    padicexp(ans, a, p, prec, modulo);
    tme += get_wall_time() - tme_tmp;

    tme_tmp = get_wall_time();
    mpz_init(ans_inverse);
    mpz_set_ui(ans_inverse, 1);
    if (p == 2) {
      padiclog_inverse(ans_inverse, a, p, prec, 2, modulo);
    } else {
      padiclog_inverse(ans_inverse, a, p, prec, 1, modulo);
    }
    tme_inverse += get_wall_time() - tme_tmp;

    if (mpz_cmp(ans, ans_inverse) != 0) {
      puts("Oups...");
    }

    padic_init2(a_flint, prec);
    padic_set_mpz(a_flint, a, ctx);
    tme_tmp = get_wall_time();
    padic_init2(ans_flint, prec);
    padic_exp(ans_flint, a_flint, ctx);
    tme_flint += get_wall_time() - tme_tmp;

    count++;
  }

  printf("Count: %u\n", count);
  printf("MyLog: %fs\n", tme);
  printf("MyInv: %fs\n", tme_inverse);
  printf("Flint: %fs\n---\n", tme_flint);

  padic_ctx_clear(ctx);

  }
  }

/*
  if (prec < 1000) {
    printf("     a = ");
    mpz_out_str(stdout, 10, a);
    printf("\n");
    printf("exp(a) = ");
    mpz_out_str(stdout, p, ans);
    printf("\n");
  }
*/
  return 0;
}
