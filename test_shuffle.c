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
#include "mxproof.h"


int main(void)
{
  unsigned int j;
  uint8_t seed[2*PQMX_SYMBYTES];

  randombytes(seed, 2*PQMX_SYMBYTES);

  poly rlwepk[2], rlwesk;
  rlwe_genkey(rlwepk, &rlwesk);

  const char *choices[] = {"Yes", "No" };  

  rlwecp *bb = (rlwecp *)  aligned_alloc(32, PQMX_NV*sizeof(rlwecp));   //ballot box
  rlwernd *rlwernds = (rlwernd *)  aligned_alloc(32, PQMX_NV*sizeof(rlwernd)); //rlwe randomnesses

  // Randomly generate ballots, encrpyt them and submit to the ballot box
  poly m;
  uint8_t paddedballot[PQMX_INDCPA_MSGBYTES];
  int chosen[4];
  for(j=0;j<PQMX_NV;j++){
    int choice = secure_random(2);
    chosen[choice]++; // store actual choices to compare with election outcome; 
    padding(paddedballot, (uint8_t*) choices[choice], strlen(choices[choice]));
    poly_frommsg(&m, paddedballot);
    rlwe_enc(&bb[j], &m, rlwepk, &rlwernds[j]);
  }
  
  // Generate secret permutation vector
  int64_t pi[PQMX_NV];
  for(j=0;j<PQMX_NV;j++)
    pi[j] = (int64_t) j+1;
    
  shuffle_array(pi, PQMX_NV);

  // Shuffle ballots using permutation vector and reencrypt
  rlwecp *shuffledbb = (rlwecp *)  aligned_alloc(32, PQMX_NV*sizeof(rlwecp)); // output of mix-net 
  rlwecp *reenc = (rlwecp *)  aligned_alloc(32, PQMX_NV*sizeof(rlwecp));  // rlwe 0-encryptions 
  memset(&m, 0, sizeof(poly));
  for(j=0;j<PQMX_NV;j++){
      shuffledbb[j] = bb[pi[j]-1];

      rlwe_enc(&reenc[j], &m, rlwepk, &rlwernds[j]);
      poly_add(&shuffledbb[j].u, &shuffledbb[j].u, &reenc[j].u);
      poly_reduce(&shuffledbb[j].u);
      poly_add(&shuffledbb[j].v, &shuffledbb[j].v, &reenc[j].v);
      poly_reduce(&shuffledbb[j].v);
  }
  
  polyvec_invntt_tomont(rlwepk,2);
  polyvec_reduce_mont(rlwepk,2);

  // generate Zero-Knowledge proof of mix-net
  proof p;
  clock_t tic = clock();
  mx_proof(&p, seed, rlwepk, rlwernds, shuffledbb, reenc, pi);
  clock_t toc = clock();
  printf("proof generated in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  // verify the proof
  tic = clock();
  int res = mx_proof_verify(&p, seed, rlwepk, bb, shuffledbb);
  toc = clock();
  printf("proof %s in %f seconds\n", (res ? "rejected" : "verified"), (double)(toc - tic) / CLOCKS_PER_SEC);

  // Start tallying
  int tally[2] = {0,0};
  uint8_t ballot[32];
  int i;
  for(j=0;j<PQMX_NV;j++){
    rlwe_dec(&m, &rlwesk, &shuffledbb[j]);
    poly_tomsg(paddedballot, &m);
    remove_padding(ballot, paddedballot);
    for(i=0;i<2;i++){
      if(!strcmp((const char*)ballot, choices[i]))
        tally[i]++;
    }
  }

  printf("Candidate | Actual | Tally \n");
  for(i=0;i<2;i++){
      printf("%9s | %6d | %d \n", choices[i], chosen[i], tally[i]);
  }

  free(bb);
  free(shuffledbb);
  free(rlwernds);
  free(reenc);
  free(p.t->t0);
  free(p.t->tm);
  free(p.z.e);
  free(p.z.em);
  free(p.z.s);
  
  return 0;
}
