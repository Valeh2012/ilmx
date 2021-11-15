#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "randombytes.h"
#include "poly.h"
#include "comm.h"
#include "util.h"
#include "rlwe.h"
#include "shortness_proof.h"

#define NTESTS 100
int main(void)
{
  unsigned int i,j,k;
  uint8_t seed[PQMX_SYMBYTES];

  poly rlwepk[2], rlwesk;
  rlwe_genkey(rlwepk, &rlwesk);
  rlwecp *bb = (rlwecp *)  aligned_alloc(32, PQMX_NV*sizeof(rlwecp)); // ballot box
  rlwernd *rlwernds = (rlwernd *)  aligned_alloc(32, PQMX_NV*sizeof(rlwernd)); // rlwe randomnesses

  uint8_t msgbuf[PQMX_INDCPA_MSGBYTES];
  randombytes(msgbuf, PQMX_INDCPA_MSGBYTES);
  poly msg;
  for(j=0;j<PQMX_NV;j++){
    memset(&msg, 0, sizeof(poly));  
    rlwe_enc(&bb[j], &msg, rlwepk, &rlwernds[j]); //store rlwe randomness used in 0-encryptions
  }

  // generate commitment key and commit ciphertexts
  randombytes(seed, PQMX_SYMBYTES);
  commkey ck;
  commrnd r;
  expand_commkey(&ck, seed);

  comm t1;
  poly *m1 = (poly*) aligned_alloc(32, PQMX_M*sizeof(poly));
  memset(m1, 0, PQMX_M*sizeof(poly));
  
  for(j=0;j<PQMX_NV;j++){
    m1[j] = bb[j].u;
    m1[PQMX_NV+j] = bb[j].v;
  }
  commit(&t1,&r, m1, &ck);

  polyvec_invntt_tomont(r.s,  PQMX_LAMBDA);
  polyvec_invntt_tomont(r.e,  PQMX_MU);
  polyvec_invntt_tomont(r.em, PQMX_M);

  polyvec_reduce_mont(r.s,  PQMX_LAMBDA);
  polyvec_reduce_mont(r.e,  PQMX_MU);
  polyvec_reduce_mont(r.em, PQMX_M);


  polyvec_invntt_tomont(m1, 2*PQMX_NV);
  polyvec_reduce_mont(m1, 2*PQMX_NV);

  // prepare shortness proof secret input
  poly *m  = (poly *) aligned_alloc(32, (4*(PQMX_LAMBDA+5*PQMX_NV)+3)*sizeof(poly));
  memset(m, 0, (4*(PQMX_LAMBDA+5*PQMX_NV)+3 )*sizeof(poly));
  
  for(k=0;k<2*PQMX_NV;k++){
    for(i=0;i<4;i++){
      for(j=0;j<PQMX_N/4;j++){
          m[4*k+i].coeffs[4*j] = r.em[k].coeffs[i*PQMX_N/4 +j];
      }
    }
  }
  for(k=0;k<PQMX_LAMBDA;k++){
    for(i=0;i<4;i++){
      for(j=0;j<PQMX_N/4;j++){
          m[4*(2*PQMX_NV+k)+i].coeffs[4*j] = r.s[k].coeffs[i*PQMX_N/4 +j];
      }
    }
  }
  for(k=0;k<PQMX_NV;k++){
    for(i=0;i<4;i++){
      for(j=0;j<PQMX_N/4;j++){
          m[4*(2*PQMX_NV+PQMX_LAMBDA+k)+i].coeffs[4*j] = rlwernds[k].r.coeffs[i*PQMX_N/4 +j];
      }
    }
    for(i=0;i<4;i++){
      for(j=0;j<PQMX_N/4;j++){
          m[4*(2*PQMX_NV+PQMX_LAMBDA+PQMX_NV+k)+i].coeffs[4*j] = rlwernds[k].eu.coeffs[i*PQMX_N/4 +j];
      }
    }
    for(i=0;i<4;i++){
      for(j=0;j<PQMX_N/4;j++){
          m[4*(2*PQMX_NV+PQMX_LAMBDA+2*PQMX_NV+k)+i].coeffs[4*j] = rlwernds[k].ev.coeffs[i*PQMX_N/4 +j];
      }
    }
  }
  
  uint8_t seed2[PQMX_SYMBYTES];
  randombytes(seed2, PQMX_SYMBYTES);
  shortness_proof p;

  for(j=0;j<2*PQMX_NV;j++){
    polyvec_invntt_tomont(&ck.bm[j], PQMX_LAMBDA);
    polyvec_reduce_mont(&ck.bm[j], PQMX_LAMBDA);

  }

  polyvec_invntt_tomont(t1.tm,2*PQMX_NV);
  polyvec_reduce_mont(t1.tm,2*PQMX_NV);

  polyvec_invntt_tomont(rlwepk,2);
  polyvec_reduce_mont(rlwepk,2);

  // generate proof of shortness of rlwe randomnesses
  clock_t tic = clock();
  shortness_prove(&p, seed2, m, rlwepk, ck.bm, t1.tm);
  clock_t toc = clock();
  printf("shortness proof generated in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);


  //verify the proof
  tic = clock();
  int res = shortness_verify(&p, seed2, rlwepk, ck.bm, t1.tm);
  toc = clock();
  printf("shortness proof %s in %f seconds\n", (res ? "rejected" : "verified"), (double)(toc - tic) / CLOCKS_PER_SEC);

  return 0;
}