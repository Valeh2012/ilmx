#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct{
  int64_t coeffs[PQMX_N];
} poly;


#define poly_tobytes PQMX_NAMESPACE(poly_tobytes)
void poly_tobytes(uint8_t r[PQMX_POLYBYTES], const poly *a);
#define poly_frombytes PQMX_NAMESPACE(poly_frombytes)
void poly_frombytes(poly *r, const uint8_t a[PQMX_POLYBYTES]);

#define poly_frommsg PQMX_NAMESPACE(poly_frommsg)
void poly_frommsg(poly *r, const uint8_t msg[PQMX_INDCPA_MSGBYTES]);
#define poly_tomsg PQMX_NAMESPACE(poly_tomsg)
void poly_tomsg(uint8_t msg[PQMX_INDCPA_MSGBYTES], const poly *r);


#define poly_uniform PQMX_NAMESPACE(poly_uniform)
void poly_uniform(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce);
#define poly_uniform_alpha PQMX_NAMESPACE(poly_uniform_alpha)
void poly_uniform_alpha(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce);
#define constant_poly_uniform_ntt PQMX_NAMESPACE(constant_poly_uniform_ntt)
void constant_poly_uniform_ntt(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce);
#define poly_nonuniform PQMX_NAMESPACE(poly_nonuniform)
void poly_nonuniform(poly *r, uint8_t mode, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce);

#define poly_uniform_delta PQMX_NAMESPACE(poly_uniform_delta)
void poly_uniform_delta(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint16_t nonce);

#define poly_ntt PQMX_NAMESPACE(poly_ntt)
void poly_ntt(poly *r);
#define poly_invntt_tomont PQMX_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *r);
#define poly_basemul_montgomery PQMX_NAMESPACE(poly_basemul_montgomery)
void poly_basemul_montgomery(poly *r, const poly *a, const poly *b);

#define poly_basemul_acc PQMX_NAMESPACE(poly_basemul_acc)
void poly_basemul_acc(int64_t r[4], const poly *a, const poly *b);

#define poly_tomont PQMX_NAMESPACE(poly_tomont)
void poly_tomont(poly *r);

#define poly_inner_prod_mont PQMX_NAMESPACE(poly_inner_prod_mont)
void poly_inner_prod_mont(poly *r, const poly *a, const poly *b);


#define poly_reduce PQMX_NAMESPACE(poly_reduce)
void poly_reduce(poly *r);
#define poly_reduce_mont PQMX_NAMESPACE(poly_reduce_mont)
void poly_reduce_mont(poly *r);

#define poly_add PQMX_NAMESPACE(poly_add)
void poly_add(poly *r, const poly *a, const poly *b);
#define poly_sub PQMX_NAMESPACE(poly_sub)
void poly_sub(poly *r, const poly *a, const poly *b);

#define poly_shift PQMX_NAMESPACE(poly_shift)
void poly_shift(poly *r);
#endif
