#include <stdint.h>
#include <string.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "symmetric.h"

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for PQMX_POLYBYTES bytes)
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tobytes(uint8_t r[PQMX_POLYBYTES], const poly *a)
{
 memcpy(r, a->coeffs, PQMX_N*sizeof(int64_t));
}

/*************************************************
* Name:        poly_frombytes
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of PQMX_POLYBYTES bytes)
**************************************************/
void poly_frombytes(poly *r, const uint8_t a[PQMX_POLYBYTES])
{
  memcpy(r->coeffs, a, PQMX_N*sizeof(int64_t));
}

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *msg: pointer to input message
**************************************************/
void poly_frommsg(poly *r, const uint8_t msg[PQMX_INDCPA_MSGBYTES])
{
  unsigned int i,j;
  int64_t mask;

#if (PQMX_INDCPA_MSGBYTES != PQMX_N/8)
#error "PQMX_INDCPA_MSGBYTES must be equal to PQMX_N/8 bytes!"
#endif

  for(i=0;i<PQMX_N/8;i++) {
    for(j=0;j<8;j++) {
      mask = -(int64_t)((msg[i] >> j)&1);
      r->coeffs[8*i+j] = mask & ((PQMX_Q+1)/2);
    }
  }
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[PQMX_INDCPA_MSGBYTES], const poly *a)
{
  unsigned int i,j;
  uint64_t t;
  for(i=0;i<PQMX_N/8;i++) {
    msg[i] = 0; 
    for(j=0;j<8;j++) {
      t  = a->coeffs[8*i+j];
      msg[i] |= t << j;
    }
  }
}


/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int64_t *r: pointer to output buffer
*              - unsigned int len: requested number of 64-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 64-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int64_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos, j;

  uint64_t t;
  ctr = pos = 0;
  while(ctr < len && pos + 8 <= buflen) {
    t = buf[pos++];
    for(j=8;j<64;j+=8)
      t |= (uint64_t)buf[pos++] << j;
    t &= (1L << 62)-1;

    if(t < PQMX_Q)
      r[ctr++] = t;
  }
  return ctr;
}



/*************************************************
* Name:        poly_uniform
*
* Description: Generate uniform random polynomial 
*
* Arguments:   - const poly *r: pointer to output polynomial
*              - const uint8_t seed[]: pointer to input buffer 
*              (assumed to be uniformly random bytes) of length PQMX_SYMBYTES
*              - uint32_t nonce: 32-bit nonce 
**************************************************/
#define POLY_UNIFORM_NBLOCKS (PQMX_POLYBYTES + XOF_BLOCKBYTES - 1)/(XOF_BLOCKBYTES)
void poly_uniform(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce)
{
  unsigned int i, ctr, off, buflen;
  xof_state state;
  uint8_t rnd[POLY_UNIFORM_NBLOCKS*XOF_BLOCKBYTES+2];
  memset(rnd,0,sizeof(rnd));
  memset(r->coeffs, 0, PQMX_POLYBYTES);

  buflen = POLY_UNIFORM_NBLOCKS*XOF_BLOCKBYTES;

  xof_absorb(&state, seed, nonce);
  xof_squeezeblocks(rnd, POLY_UNIFORM_NBLOCKS, &state);
  
  ctr = rej_uniform(r->coeffs, PQMX_N, rnd, buflen);
  
  while(ctr<PQMX_N) {
    off = buflen % 3;
    for(i = 0; i < off; ++i)
      rnd[i] = rnd[buflen - off + i];
    xof_squeezeblocks(rnd+off, 1, &state);
    buflen = XOF_BLOCKBYTES + off;
    ctr+= rej_uniform(r->coeffs + ctr, PQMX_N - ctr, rnd, buflen);
  }

  poly_reduce(r);
}

/*************************************************
* Name:        poly_uniform_alpha
*
* Description: Generate uniform random polynomial in S_q
*
* Arguments:   - const poly *r: pointer to output polynomial
*              - const uint8_t seed[]: pointer to input buffer 
*              (assumed to be uniformly random bytes) of length PQMX_SYMBYTES
*              - uint32_t nonce: 32-bit nonce 
**************************************************/
#define POLY_UNIFORM_ALPHA_NBLOCKS (PQMX_POLYBYTES/PQMX_L + XOF_BLOCKBYTES - 1)/(XOF_BLOCKBYTES)
void poly_uniform_alpha(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce)
{
  unsigned int i,ctr,off, buflen;
  xof_state state;
  uint8_t rnd[POLY_UNIFORM_ALPHA_NBLOCKS*XOF_BLOCKBYTES+2];
  memset(rnd,0,sizeof(rnd));
  memset(r->coeffs, 0, PQMX_POLYBYTES);

  buflen = POLY_UNIFORM_ALPHA_NBLOCKS*XOF_BLOCKBYTES;

  xof_absorb(&state, seed, nonce);
  xof_squeezeblocks(rnd, POLY_UNIFORM_ALPHA_NBLOCKS, &state);
  
  ctr = rej_uniform(r->coeffs, PQMX_N/PQMX_L, rnd, buflen);
  
  while(ctr<PQMX_N/PQMX_L) {
    off = buflen % 3;
    for(i = 0; i < off; ++i)
      rnd[i] = rnd[buflen - off + i];
    xof_squeezeblocks(rnd+off, 1, &state);
    buflen = XOF_BLOCKBYTES + off;
    ctr+= rej_uniform(r->coeffs + ctr, PQMX_N - ctr, rnd, buflen);
  }

  poly_reduce(r);
}

/*************************************************
* Name:        constant_poly_uniform_ntt
*
* Description: Generate uniform random constant polynomial (integer) in R_q in NTT domain 
*
* Arguments:   - const poly *r: pointer to output polynomial
*              - const uint8_t seed[]: pointer to input buffer 
*              (assumed to be uniformly random bytes) of length PQMX_SYMBYTES
*              - uint32_t nonce: 32-bit nonce 
**************************************************/
void constant_poly_uniform_ntt(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce)
{
  unsigned int i,ctr=0,buflen;
  xof_state state;
  uint8_t rnd[XOF_BLOCKBYTES];
  memset(rnd,0,sizeof(rnd));
  memset(r->coeffs, 0, PQMX_POLYBYTES);

  buflen = XOF_BLOCKBYTES;
  int64_t coeff=0;
  xof_absorb(&state, seed, nonce);
  
  while(!ctr) {
    xof_squeezeblocks(rnd, 1, &state);
    ctr += rej_uniform(&coeff, 1, rnd, buflen);
  }
  for(i=0;i<PQMX_N;i+=PQMX_N/PQMX_L)
    r->coeffs[i] = coeff;
  
  poly_reduce(r);
}


/*************************************************
* Name:        poly_nonuniform
*
* Description: Generate ternary polynomial with zero coefficient probability mode
*
* Arguments:   - const poly *r: pointer to output polynomial
*              - uint8_t mode: probability mode
*              - const uint8_t seed[]: pointer to input buffer 
*              (assumed to be uniformly random bytes) of length PQMX_SYMBYTES
*              - uint32_t nonce: 32-bit nonce 
**************************************************/
void poly_nonuniform(poly *r, uint8_t mode, const uint8_t seed[PQMX_SYMBYTES], uint32_t nonce)
{
  unsigned int i;
  uint8_t rnd[PQMX_N];
  memset(rnd,0,sizeof(rnd));
  prf(rnd, PQMX_N, seed, nonce);
  int64_t t = 0;
  for(i=0;i<PQMX_N;i++) {
    t = rnd[i] & 1;
    r->coeffs[i] = rnd[i]< mode ? (1 - 2*t) : 0;
  }
}


/*************************************************
* Name:        poly_uniform_delta
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-(delta - 1), delta] by unpacking output stream
*              of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - const poly *r: pointer to output polynomial
*              - const uint8_t seed[]: pointer to input buffer 
*              (assumed to be uniformly random bytes) of length PQMX_SYMBYTES
*              - uint16_t nonce: 16-bit nonce 
**************************************************/
#define POLY_UNIFORM_DELTA_NBLOCKS ((PQMX_POLYCOMPRESSEDBYTES + XOF_BLOCKBYTES - 1)/XOF_BLOCKBYTES)
void poly_uniform_delta(poly *r, const uint8_t seed[PQMX_SYMBYTES], uint16_t nonce)
{
  uint8_t buf[POLY_UNIFORM_DELTA_NBLOCKS*XOF_BLOCKBYTES];
  xof_state state;
  unsigned int i;


  xof_absorb(&state, seed, nonce);
  xof_squeezeblocks(buf, POLY_UNIFORM_DELTA_NBLOCKS, &state);

  

#if PQMX_DELTA == (1L<<47)  

  for(i=0;i<PQMX_N; i++){
    r->coeffs[i] = buf[6*i+0];
    r->coeffs[i] |= (uint64_t) buf[6*i+1] << 8;
    r->coeffs[i] |= (uint64_t) buf[6*i+2] << 16;
    r->coeffs[i] |= (uint64_t) buf[6*i+3] << 24;
    r->coeffs[i] |= (uint64_t) buf[6*i+4] << 32;
    r->coeffs[i] |= (uint64_t) buf[6*i+5] << 40;
    

    r->coeffs[i] = PQMX_DELTA -  r->coeffs[i];
  }


#elif PQMX_DELTA == (1L<<43)  

  for(i=0;i<PQMX_N/2; i++){
    r->coeffs[2*i+0] = buf[11*i+0];
    r->coeffs[2*i+0] |= (uint64_t) buf[11*i+1] << 8;
    r->coeffs[2*i+0] |= (uint64_t) buf[11*i+2] << 16;
    r->coeffs[2*i+0] |= (uint64_t) buf[11*i+3] << 24;
    r->coeffs[2*i+0] |= (uint64_t) buf[11*i+4] << 32;
    r->coeffs[2*i+0] |= (uint64_t) buf[11*i+5] << 40;
    r->coeffs[2*i+0] &= 0xFFFFFFFFFFF;

    r->coeffs[2*i+1] = buf[11*i+5] >> 4;
    r->coeffs[2*i+1] |= (uint64_t) buf[11*i+6] << 4;
    r->coeffs[2*i+1] |= (uint64_t) buf[11*i+7] << 12;
    r->coeffs[2*i+1] |= (uint64_t) buf[11*i+8] << 20;
    r->coeffs[2*i+1] |= (uint64_t) buf[11*i+9] << 28;
    r->coeffs[2*i+1] |= (uint64_t) buf[11*i+10] << 36;
    r->coeffs[2*i+1] &= 0xFFFFFFFFFFF;


    r->coeffs[2*i+0] = PQMX_DELTA -  r->coeffs[2*i+0];
    r->coeffs[2*i+1] = PQMX_DELTA -  r->coeffs[2*i+1];
  }

#endif

}

/*************************************************
* Name:        poly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint64_t *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r)
{
  ntt(r->coeffs);
  poly_reduce(r);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT)
*              of a polynomial in place;
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint64_t *a: pointer to in/output polynomial
**************************************************/
void poly_invntt_tomont(poly *r)
{
  invntt(r->coeffs);
}

/*************************************************
* Name:        poly_basemul_montgomery
*
* Description: Multiplication of two polynomials in NTT domain
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_montgomery(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<PQMX_N/8;i++) {
    basemul(&r->coeffs[8*i], &a->coeffs[8*i], &b->coeffs[8*i], zetas[512+i]);
    basemul(&r->coeffs[8*i+4], &a->coeffs[8*i+4], &b->coeffs[8*i+4], -zetas[512+i]);
  }
}


/*************************************************
* Name:        poly_basemul_acc
*
* Description: Inner product in F_q^4. 
*
* Arguments:   - int64_t *r: pointer to output array
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial array
**************************************************/
void poly_basemul_acc(int64_t r[4], const poly *a, const poly b[4])
{
  unsigned int i,j;
  int64_t tmp[4];
  memset(r, 0, 4*sizeof(int64_t));
  for(i=0;i<4;i++){
    for(j=0; j< PQMX_N/4; j++){
      scalar_field_mul(tmp, a->coeffs[i*PQMX_N/4+j], &b[i].coeffs[4*j]);
      r[0] = barrett_reduce(r[0] + tmp[0]);
      r[1] = barrett_reduce(r[1] + tmp[1]);
      r[2] = barrett_reduce(r[2] + tmp[2]);
      r[3] = barrett_reduce(r[3] + tmp[3]);
    }
  }
}

/*************************************************
* Name:        poly_tomont
*
* Description: Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_tomont(poly *r)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = montgomery_reduce((__int128)r->coeffs[i]*MONT2);
}

/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *r)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = barrett_reduce(r->coeffs[i]);
}

/*************************************************
* Name:        poly_reduce_mont
*
* Description: Applies Montgomery reduction to all coefficients of a polynomial
*              for details of the Montgomery reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_reduce_mont(poly *r)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = montgomery_reduce( (__int128) r->coeffs[i] );
}


/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials; no modular reduction is performed
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials; no modular reduction is performed
*
* Arguments: - poly *r:       pointer to output polynomial+
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<PQMX_N;i++)
    r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

/*************************************************
* Name:        poly_shift
*
* Description: Inplace Shift polynomial coefficients right by one and negate
*              the new leading coefficient. 
*
* Arguments: - poly *r:       pointer to input polynomial
**************************************************/
void poly_shift(poly *r){
  int64_t tmp = r->coeffs[PQMX_N-1];
  unsigned int i;
  for(i=PQMX_N-1;i>0;i--){
    r->coeffs[i] = r->coeffs[i-1];
  }
  r->coeffs[0] = -1*tmp;
}