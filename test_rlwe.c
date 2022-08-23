#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "randombytes.h"
#include "poly.h"
#include "rlwe.h"

#define NTESTS 1000

int main(void)
{
  int i;
  uint8_t seed[PQMX_SYMBYTES];

  randombytes(seed, PQMX_SYMBYTES);

  poly rlwepk[2], rlwesk;
  rlwe_genkey(rlwepk, &rlwesk);
  rlwernd rnd;
  rlwecp cp;
  
  uint8_t msg[PQMX_INDCPA_MSGBYTES] = "Hello world!";
  printf("Plaintext: %s\n", (const char *)msg);
  poly m;
  poly_frommsg(&m, msg);

  clock_t tic = clock();
  for(i=0;i<NTESTS;i++)
    rlwe_enc(&cp, &m, rlwepk, &rnd);
  clock_t toc = clock();
  printf("rlwe encryption time in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC / NTESTS);
  
  memset(msg, 0, PQMX_INDCPA_MSGBYTES);
  // verify
  tic = clock();
  for(i=0;i<NTESTS;i++)
    rlwe_dec(&m, &rlwesk, &cp);
  toc = clock();
  printf("rlwe decryption time in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC / NTESTS);

  poly_tomsg(msg, &m);
  printf("Decrypted ciphertext: %s\n",(const char *)msg);

  uint8_t msg2[PQMX_INDCPA_MSGBYTES];
  memset(msg2, 0, PQMX_INDCPA_MSGBYTES);

  poly_frommsg(&m, msg);
  poly m2;

  int failures = 0;
  for(i=0;i<NTESTS;i++){
    rlwe_enc(&cp, &m, rlwepk, &rnd);
    rlwe_dec(&m2, &rlwesk, &cp);
    poly_tomsg(msg2, &m2);
    failures += memcmp(msg, msg2, PQMX_INDCPA_MSGBYTES) != 0;
  }
  printf("number of failures: %d\n", failures);
  return 0;
}
