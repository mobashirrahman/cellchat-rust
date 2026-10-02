/* Reference oracle for the F80 port: performs the same operations with the real x87
 * 80-bit `long double` that R uses for `mean`, and prints the bit-exact double result.
 *
 * Build:  gcc -O2 -o f80ref gen_f80_ref.c && ./f80ref > tests/fixtures/f80_ref.txt
 *
 * The output is the ground truth the Rust port is differentially tested against. Hand
 * reasoning about the 128-bit accumulator layout found nothing; comparing against the
 * real thing found the `sub_same_sign` exponent-underflow bug immediately.
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/* R's own generator: Mersenne-Twister via LCG is not needed here; a small explicit
 * LCG keeps the reference independent of the port's RNG. */
static uint64_t s_ = 0x853c49e6748fea9bULL;
static uint64_t nx(void) {
    s_ ^= s_ << 13; s_ ^= s_ >> 7; s_ ^= s_ << 17; return s_;
}
static double u01(void) { return (double)(nx() >> 11) * (1.0 / 9007199254740992.0); }

int main(void) {
    printf("#x87 long double reference; %zu-bit mantissa\n", (size_t)(LDBL_MANT_DIG));
    setvbuf(stdout, NULL, _IOFBF, 1 << 20);

    /* 1. Pairwise add/sub over a spread of magnitudes, exponents and signs. */
    for (int i = 0; i < 3000; i++) {
        double a, b;
        switch (i % 6) {
            case 0: a = u01();                b = u01();                break;
            case 1: a = u01() * 1e-8;         b = u01() * 1e8;          break;
            case 2: a = (double)(int64_t)(nx() % 1000000); b = u01() * 1000.0; break;
            case 3: a = 1.0 / (1.0 + u01());  b = 1.0 / (1.0 + u01());  break;  /* near 1 */
            case 4: a = ldexpl(u01(), 40);     b = ldexpl(u01(), -40);   break;  /* wide gap */
            default: a = u01() * 1e16;        b = 1.0;                  break;  /* cancellation */
        }
        if (i & 1) { a = -a; b = -b; }
        long double la = a, lb = b;
        printf("add\t%.17g\t%.17g\t%.17g\n", a, b, (double)(la + lb));
        printf("sub\t%.17g\t%.17g\t%.17g\n", a, b, (double)(la - lb));
    }

    static const int lens[] = {1, 2, 3, 4, 5, 8, 16, 25, 40, 64, 100, 257, 1000, 5000, 50000};

    /* 1b. Products, including chains: R's `prod` is an LDOUBLE accumulation, and
     *    computeExpr_coreceptor/_agonist/_antagonist end in apply(..., 2, prod). */
    for (int i = 0; i < 3000; i++) {
        double a, b;
        switch (i % 5) {
            case 0: a = 1.0 + u01();           b = 1.0 + u01();           break;
            case 1: a = ldexpl(u01(), -30);     b = ldexpl(u01(), 30);     break;
            case 2: a = 1.0 / (1.0 + u01());    b = 1.0 / (1.0 + u01());    break;
            case 3: a = 1e300;                  b = 1e-300;                break;
            default: a = 2.0;                   b = ldexpl(u01(), -20);     break;
        }
        if (i & 1) b = -b;
        long double la = a, lb = b;
        printf("mul\t%.17g\t%.17g\t%.17g\n", a, b, (double)(la * lb));
    }
    /* 1c. Chained products, mirroring prod(). */
    for (unsigned li = 0; li < sizeof(lens)/sizeof(lens[0]); li++) {
        int n = lens[li];
        for (int rep = 0; rep < 8; rep++) {
            double *v = malloc((size_t)n * sizeof(double));
            long double acc = 1.0L;
            for (int k = 0; k < n; k++) {
                double x = (rep % 2) ? (1.0 + u01()) : (1.0 / (1.0 + u01()));
                v[k] = x; acc *= x;
            }
            printf("prod\t");
            for (int k = 0; k < n; k++) printf("%s%.17g", k ? "," : "", v[k]);
            printf("\t%.17g\n", (double)acc);
            free(v);
        }
    }

    /* 2. Chained sums of realistic lengths: this is R's `real_mean` pass 1, and the
     *    case that actually matters (the second pass is covered by the stats corpus). */
    for (unsigned li = 0; li < sizeof(lens)/sizeof(lens[0]); li++) {
        int n = lens[li];
        for (int rep = 0; rep < 12; rep++) {
            double *v = malloc((size_t)n * sizeof(double));
            long double acc = 0.0L;
            for (int k = 0; k < n; k++) {
                double x;
                switch (rep % 5) {
                    case 0:  x = u01(); break;                              /* [0,1)   */
                    case 1:  x = (nx() % 2) ? 0.0 : u01(); break;           /* sparse  */
                    case 2:  x = ldexpl(u01(), (int)(nx() % 60) - 30); break;/* wide    */
                    case 3:  x = (double)(int64_t)(nx() % 100) - 50.0; break;/* signed  */
                    default: x = (nx() % 10 == 0) ? 1e16 : 1.0 / (1.0 + u01()); break;
                }
                v[k] = x;
                acc += x;
            }
            printf("chain\t");
            for (int k = 0; k < n; k++) printf("%s%.17g", k ? "," : "", v[k]);
            printf("\t%.17g\n", (double)acc);
            free(v);
        }
    }
    /* 3. Adversarial chains: interleaved catastrophic cancellation.
     *
     * The `chain` rows above use `u01()`-scale values, and at 80 bits almost none of them change
     * answer when the sum is reassociated -- measured at 2 of 156 on the row set above, because
     * 11 extra mantissa bits absorb most reordering. That is a real and reassuring property, but
     * it is not evidence that the order is right, because the order is barely observable there.
     *
     * These rows are built to make it observable. The pattern `+B, +s, -B, +s` loses every `s` to
     * the `B` that precedes it when accumulated left to right, and a tree reduction recovers some
     * of them. Both answers are computed by the same x87 type, so the recorded value is R's
     * answer for R's order, and the port has to reproduce it.
     */
    {
        static const double bases[] = {1e16, 1e17, 1e300, 1e-300, 9.007199254740992e15};
        static const double smalls[] = {1.0, 0.5, 3.0, 1e-8, 7.0};
        for (int bi = 0; bi < (int)(sizeof(bases)/sizeof(bases[0])); bi++) {
            for (int si = 0; si < (int)(sizeof(smalls)/sizeof(smalls[0])); si++) {
                for (int reps = 2; reps <= 12; reps++) {
                    int n = reps * 4;
                    double *v = malloc((size_t)n * sizeof(double));
                    long double acc = 0.0L;
                    for (int k = 0; k < n; k++) {
                        switch (k % 4) {
                            case 0: v[k] =  bases[bi]; break;
                            case 1: v[k] =  smalls[si]; break;
                            case 2: v[k] = -bases[bi]; break;
                            default: v[k] = smalls[si]; break;
                        }
                        acc += v[k];
                    }
                    printf("chain_cancel\t");
                    for (int k = 0; k < n; k++) printf("%s%.17g", k ? "," : "", v[k]);
                    printf("\t%.17g\n", (double)acc);
                    free(v);
                }
            }
        }
    }
    return 0;
}
