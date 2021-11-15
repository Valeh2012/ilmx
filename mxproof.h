#ifndef MXPROOF_H
#define MXPROOF_H

#include <stdint.h>
#include "params.h"
#include "rlwe.h"
#include "comm.h"

typedef struct {
  uint8_t hash[PQMX_SYMBYTES];
  commrnd z;
  comm *t;
  void *sp;
} proof;

void mx_proof(proof *p, const uint8_t rho[2*PQMX_SYMBYTES], const poly rlwepk[2], const rlwernd rnd[PQMX_NV], const rlwecp out[PQMX_NV], const rlwecp reenc[PQMX_NV], const int64_t pi[PQMX_NV]);
int mx_proof_verify(const proof *p, const uint8_t rho[2*PQMX_SYMBYTES], poly rlwepk[2], const rlwecp in[PQMX_NV], const rlwecp out[PQMX_NV]);

#endif