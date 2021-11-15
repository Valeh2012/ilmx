#include "stdlib.h"
#include "string.h"
#include "comm.h"
#include "symmetric.h"
#include "randombytes.h"
#include "util.h"

/*************************************************
* Name:        expand_commkey or shortness_expand_commkey
*
* Description: Expand a commitment key from a seed
*
* Arguments:   - commkey *ck: pointer to output commitment key
*              - const uint8_t seed[PQMX_SYMBYTES]: seed buffer
**************************************************/
#ifndef SHORT
void expand_commkey(commkey *ck, const uint8_t seed[PQMX_SYMBYTES])
{
#else
void shortness_expand_commkey(commkey *ck, const uint8_t seed[PQMX_SYMBYTES])
{
#endif
  unsigned int i,j;
  uint8_t buf[PQMX_SYMBYTES];
  memcpy(buf, seed, PQMX_SYMBYTES);

  ck->b0 = (poly *) aligned_alloc(32, PQMX_MU*PQMX_LAMBDA*sizeof(poly));
  ck->bt = (poly *) aligned_alloc(32, PQMX_MU*PQMX_M*sizeof(poly));
  ck->bm = (poly *) aligned_alloc(32, PQMX_M*PQMX_LAMBDA*sizeof(poly));
  
  for(i=0;i<PQMX_MU;i++) {
    for(j=0;j<PQMX_LAMBDA;j++) {
      poly_uniform(&ck->b0[i] + j, buf, (i<<8) + j);
      poly_ntt(&ck->b0[i] + j);
    }
  }
  prf(buf, PQMX_SYMBYTES, buf, (uint8_t)(PQMX_MU << 8) + PQMX_LAMBDA);

  for(i=0;i<PQMX_MU;i++) {
    for(j=0;j<PQMX_M;j++) {
      poly_uniform(&ck->bt[i] + j, buf, (i<<8) + j);
      poly_ntt(&ck->bt[i] + j);
    }
  }
  prf(buf, PQMX_SYMBYTES, buf, (uint8_t)(PQMX_MU << 8) + PQMX_M);
  
  for(i=0;i<PQMX_M;i++) {
    for(j=0;j<PQMX_LAMBDA;j++) {
      poly_uniform(&ck->bm[i] + j, buf, (i<<8) + j);
      poly_ntt(&ck->bm[i] + j);
    }
  }
}

/*************************************************
* Name:        commit or shortness_commit
*
* Description: Commit to a message using commitment key
*
* Arguments:   - comm *t: pointer to output commitment
*              - const poly *msg: pointer to message polynomial to be committed
*              - const commkey *ck: pointer to commitment key
**************************************************/
#ifndef  SHORT
void commit(comm *t, commrnd *r, const poly *msg, const commkey *ck) {
#else
void shortness_commit(comm *t, commrnd *r, const poly *msg, const commkey *ck) {
#endif
  int i;
  uint32_t nonce = 0;
  uint8_t seed[PQMX_SYMBYTES];
  poly* tag = (poly *) aligned_alloc(32, PQMX_MU*sizeof(poly));

  r->s = (poly *) aligned_alloc(32, PQMX_LAMBDA*sizeof(poly));
  r->e = (poly *) aligned_alloc(32, PQMX_MU*sizeof(poly));
  r->em = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));

  t->t0 = (poly *) aligned_alloc(32, PQMX_MU*sizeof(poly));
  t->tm = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));

  memset(seed,0,PQMX_SYMBYTES);
  randombytes(seed,PQMX_SYMBYTES);
  for(i=0;i<PQMX_LAMBDA;i++) {
    poly_nonuniform(&r->s[i], 0xA0, seed, nonce++);
  }
  for(i=0;i<PQMX_MU;i++) {
    poly_nonuniform(&r->e[i], 0xA0, seed, nonce++);
  }
  for(i=0;i<PQMX_M;i++) {
    poly_nonuniform(&r->em[i], 0xA0, seed, nonce++);
  }

  polyvec_ntt(r->s,  PQMX_LAMBDA);
  polyvec_ntt(r->e,  PQMX_MU);
  polyvec_ntt(r->em, PQMX_M);


  for(i=0;i<PQMX_MU;i++){
    polyvec_basemul_acc_montgomery(&t->t0[i],&ck->b0[i*PQMX_LAMBDA],r->s, PQMX_LAMBDA);
    poly_tomont(&t->t0[i]);
  }
  for(i=0;i<PQMX_M;i++){
    polyvec_basemul_acc_montgomery(&t->tm[i],&ck->bm[i*PQMX_LAMBDA],r->s, PQMX_LAMBDA);
    poly_tomont(&t->tm[i]);
  }
  for(i=0;i<PQMX_MU;i++){
    polyvec_basemul_acc_montgomery(&tag[i],&ck->bt[i*PQMX_M],r->em, PQMX_M);
    poly_tomont(&tag[i]);
  }

  polyvec_add(t->t0,t->t0,tag, PQMX_MU);
  polyvec_reduce(t->t0, PQMX_MU);
  polyvec_add(t->t0,t->t0,r->e, PQMX_MU);
  polyvec_reduce(t->t0, PQMX_MU);
  polyvec_add(t->tm,t->tm,r->em, PQMX_M);
  polyvec_reduce(t->tm, PQMX_M);
  polyvec_add(t->tm,t->tm,msg, PQMX_M);
  polyvec_reduce(t->tm, PQMX_M);
  
  free(tag);
}
