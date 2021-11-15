#ifndef RLWE_H
#define RLWE_H

#include <stdlib.h>
#include "params.h"
#include "poly.h"

typedef struct{
  poly r;
  poly eu;
  poly ev;
}rlwernd;

typedef struct{
    poly u;
    poly v;
}rlwecp;

#define rlwe_genkey PQMX_NAMESPACE(rlwe_keygen)
void rlwe_genkey(poly pk[2], poly * sk);

#define rlwe_enc PQMX_NAMESPACE(rlwe_enc)
void rlwe_enc(rlwecp* cp, const poly *msg, const poly pk[2], rlwernd *rnd);

#define rlwe_dec PQMX_NAMESPACE(rlwe_dec)
void rlwe_dec(poly *msg, const poly * sk, const rlwecp* cp);

#endif