#ifndef SHORTNESSPROOF_H
#define SHORTNESSPROOF_H

#include <stdlib.h>
#include "rlwe.h"
#include "comm.h"

typedef struct {
  uint8_t hash[PQMX_SYMBYTES];
  commrnd z;
  poly h;
  comm *t;
} shortness_proof;

void shortness_prove(shortness_proof *p, const uint8_t rho[PQMX_SYMBYTES], poly *msg, const poly pk[2], poly *A, poly *u);

int shortness_verify(shortness_proof *p, const uint8_t rho[PQMX_SYMBYTES], const poly pk[2], poly *A, poly *u);

#endif