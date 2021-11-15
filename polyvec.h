#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"


#define polyvec_tobytes PQMX_NAMESPACE(polyvec_tobytes)
void polyvec_tobytes(uint8_t *r, const poly *a, int vlen);
#define polyvec_frombytes PQMX_NAMESPACE(polyvec_frombytes)
void polyvec_frombytes(poly *r, const uint8_t *a, int vlen);

#define polyvec_ntt PQMX_NAMESPACE(polyvec_ntt)
void polyvec_ntt(poly *r, int vlen);
#define polyvec_invntt_tomont PQMX_NAMESPACE(polyvec_invntt_tomont)
void polyvec_invntt_tomont(poly *r, int vlen);

#define polyvec_basemul_acc_montgomery PQMX_NAMESPACE(polyvec_basemul_acc_montgomery)
void polyvec_basemul_acc_montgomery(poly *r, const poly *a, const poly *b, int vlen);

#define polyvec_reduce PQMX_NAMESPACE(polyvec_reduce)
void polyvec_reduce(poly *r, int vlen);

#define polyvec_reduce_mont PQMX_NAMESPACE(polyvec_reduce_mont)
void polyvec_reduce_mont(poly *r, int vlen);

#define polyvec_add PQMX_NAMESPACE(polyvec_add)
void polyvec_add(poly *r, const poly *a, const poly *b, int vlen);

#endif
