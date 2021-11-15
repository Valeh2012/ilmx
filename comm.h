#ifndef COMM_H
#define COMM_H

#include "params.h"
#include "polyvec.h"

typedef struct {
  poly *t0;     // length mu
  poly *tm;     // length m
} comm;

typedef struct {
  poly *s;      // length lambda
  poly *e;      // length mu
  poly *em;     // length m
} commrnd;

typedef struct {
  poly * b0;   // length mu x lambda
  poly * bt;   // length mu x m
  poly * bm;   // length m x lambda
} commkey;

#ifndef SHORT
void expand_commkey(commkey *ck, const uint8_t seed[PQMX_SYMBYTES]);
void commit(comm *t, commrnd *r, const poly *msg, const commkey *ck);
#else
void shortness_expand_commkey(commkey *ck, const uint8_t seed[PQMX_SYMBYTES]);
void shortness_commit(comm *t, commrnd *r, const poly *msg, const commkey *ck);
#endif
#endif