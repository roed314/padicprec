#include <math.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void padiclog(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec) {
  unsigned long i, N, Np, tmp, trunc, step;
  mpz_t modulo;
  mpz_t r, f, arg, trunc_mod, h, hpow, tmp2;
  mpz_t d, inv;
  mpz_t *num;
  mpz_t *denom;

  mpz_init(tmp2);
  mpz_fdiv_r_ui(tmp2, a, p);
  if (mpz_cmp_ui(tmp2, 1) != 0) {
    puts("Pas congru");
    return;
  }

  mpz_init(modulo);
  mpz_ui_pow_ui(modulo, p, prec);

  mpz_init_set(arg, a);

  N = prec;
  while(1) {
    tmp = prec + (unsigned long)(log(N)/log(p));
    if (tmp == N) break;
    N = tmp;
  }

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

    if (trunc > prec) break;

    mpz_mul(trunc_mod, trunc_mod, trunc_mod);
    trunc <<= 1;
    N >>= 1;
  }

  mpz_fdiv_r(ans, ans, modulo);
}

int main(){
  unsigned long p = 2;
  unsigned long prec = 1600000;

  double tme;
  mpz_t modulo, a, ans;
  gmp_randstate_t state;


  mpz_init(modulo);
  mpz_ui_pow_ui(modulo, p, prec-1);
  mpz_init(a);
  gmp_randinit_default(state);
  mpz_urandomm(a, state, modulo);
  mpz_mul_ui(a, a, p);
  mpz_add_ui(a, a, 1);

  tme = get_wall_time();
  mpz_init(ans);
  padiclog(ans, a, p, prec);
  tme = get_wall_time() - tme;

  printf("Wall time: %fs\n", tme);
  if (prec < 100) {
    printf("log = ");
    mpz_out_str(stdout, p, ans);
    printf("\n");
  }
  return 0;
}
