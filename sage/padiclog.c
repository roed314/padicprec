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

void padiclog(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo) {
  unsigned long i, N, saveN, Np, tmp, trunc, step;
  mpz_t f, arg, trunc_mod, h, hpow, tmp2;
  mpz_t d, inv;
  mpz_t *num;
  mpz_t *denom;

  mpz_init(tmp2);
/*
  mpz_fdiv_r_ui(tmp2, a, p);
  if (mpz_cmp_ui(tmp2, 1) != 0) {
    puts("Pas congru");
    return;
  }
*/

  mpz_init_set(arg, a);

  N = prec;
  while(1) {
    tmp = prec + (unsigned long)(log(N)/log(p));
    if (tmp == N) break;
    N = tmp;
  }
  saveN = N;

  mpz_set_ui(ans, 0);

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

    mpz_ui_sub(f, 2, f);
    mpz_mul(arg, arg, f);

    for (i = 0; i < N; i++) {
      mpz_set_ui(num[i], 1);
      mpz_set_ui(denom[i], i+1);
    }
    step = 1;
    mpz_ui_sub(h, 1, f);
    mpz_set(hpow, h);
    while(step < N) {
      for (i = 0; i < N - step; i += step << 1) {
        mpz_mul(tmp2, hpow, num[i+step]);
        mpz_mul(tmp2, tmp2, denom[i]);
        mpz_mul(num[i], num[i], denom[i+step]);
        mpz_add(num[i], num[i], tmp2);
        mpz_mul(denom[i], denom[i], denom[i+step]);
      }
      step <<= 1;
      mpz_mul(hpow, hpow, hpow);
    }

    Np = N; tmp = 0;
    while(Np > 0) { Np /= p; tmp += Np; }
    mpz_ui_pow_ui(d, p, tmp);
    mpz_divexact(tmp2, num[0], d);
    mpz_mul(tmp2, h, tmp2);

    mpz_divexact(denom[0], denom[0], d);
    mpz_gcdext(d, inv, NULL, denom[0], modulo);
    mpz_mul(tmp2, tmp2, inv);

    mpz_add(ans, ans, tmp2);

    }

    if (trunc > prec) break;

    mpz_mul(trunc_mod, trunc_mod, trunc_mod);
    trunc <<= 1;
    N >>= 1;
  }

  mpz_fdiv_r(ans, ans, modulo);

    /* We clear memory */
    mpz_clear(arg);
    mpz_clear(f);
    mpz_clear(trunc_mod);
    mpz_clear(h);
    mpz_clear(hpow);
    mpz_clear(tmp2);
    mpz_clear(d);
    mpz_clear(inv);
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

  double tme_tmp, tme, tme_flint;
  mpz_t modulo, a, ans;
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
    mpz_mul_ui(a, a, p);
    mpz_add_ui(a, a, 1);

    tme_tmp = get_wall_time();
    padiclog(ans, a, p, prec, modulo);
    tme += get_wall_time() - tme_tmp;

    padic_init2(a_flint, prec);
    padic_set_mpz(a_flint, a, ctx);
    tme_tmp = get_wall_time();
    padic_init2(ans_flint, prec);
    padic_log(ans_flint, a_flint, ctx);
    padic_clear(a_flint);
    padic_clear(ans_flint);
    tme_flint += get_wall_time() - tme_tmp;

    count++;
  }

  printf("Count: %u\n", count);
  printf("MyLog: %fs\n", tme);
  printf("Flint: %fs\n---\n", tme_flint);

  padic_ctx_clear(ctx);

  // printf("Gain: %f%%\n", 100*(1 - tme/tme_flint));

  }
  }

/*
  if (prec < 100) {
    printf("log = ");
    mpz_out_str(stdout, p, ans);
    printf("\n");
  }
*/
  return 0;
}
