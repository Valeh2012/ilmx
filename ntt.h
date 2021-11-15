#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define zetas PQMX_NAMESPACE(zetas)
extern const int64_t zetas[1024];

#define ntt PQMX_NAMESPACE(ntt)
void ntt(int64_t poly[4096]);

#define invntt PQMX_NAMESPACE(invntt)
void invntt(int64_t poly[4096]);

#define basemul PQMX_NAMESPACE(basemul)
void basemul(int64_t r[4], const int64_t a[4], const int64_t b[4], int64_t zeta);

#define scalar_field_mul PQMX_NAMESPACE(scalar_field_mul)
void scalar_field_mul(int64_t r[4], const int64_t a, const int64_t b[4]);

int64_t fqmul(int64_t a, int64_t b);

#endif
