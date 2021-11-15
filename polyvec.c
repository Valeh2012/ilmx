#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"


/*************************************************
* Name:        polyvec_tobytes
*
* Description: Serialize vector of polynomials
*
* Arguments:   - uint8_t *r: pointer to output byte array
*              - const poly *a: pointer to input vector of polynomials
*              - int vlen: vector dimension 
**************************************************/
void polyvec_tobytes(uint8_t *r, const poly *a, int vlen)
{
  int i;
  for(i=0;i<vlen;i++)
    poly_tobytes(r+i*PQMX_POLYBYTES, &a[i]);
}

/*************************************************
* Name:        polyvec_frombytes
*
* Description: De-serialize vector of polynomials;
*              inverse of polyvec_tobytes
*
* Arguments:   - uint8_t *r:       pointer to output byte array
*              - const poly *a: pointer to input vector of polynomials
*              - int vlen: vector dimension 
**************************************************/
void polyvec_frombytes(poly *r, const uint8_t *a, int vlen)
{
  int i;
  for(i=0;i<vlen;i++)
    poly_frombytes(&r[i], a+i*PQMX_POLYBYTES);
}

/*************************************************
* Name:        polyvec_ntt
*
* Description: Apply forward NTT to all elements of a vector of polynomials
*
* Arguments:   - poly *r: pointer to in/output vector of polynomials
*              - int vlen: vector dimension 
**************************************************/
void polyvec_ntt(poly *r, int vlen)
{
  int i;
  for(i=0;i<vlen;i++)
    poly_ntt(&r[i]);
}

/*************************************************
* Name:        polyvec_invntt_tomont
*
* Description: Apply inverse NTT to all elements of a vector of polynomials
*              and multiply by Montgomery factor 2^16
*
* Arguments:   - poly *r: pointer to in/output vector of polynomials
*              - int vlen: vector dimension 
**************************************************/
void polyvec_invntt_tomont(poly *r, int vlen)
{
  int i;
  for(i=0;i<vlen;i++)
    poly_invntt_tomont(&r[i]);
}

/*************************************************
* Name:        polyvec_basemul_acc_montgomery
*
* Description: Multiply elements of a and b in NTT domain, accumulate into r,
*              and multiply by 2^-16.
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input vector of polynomials
*            - const poly *b: pointer to second input vector of polynomials
*            - int vlen: vector dimension 
**************************************************/
void polyvec_basemul_acc_montgomery(poly *r, const poly *a, const poly *b, int vlen)
{
  int i;
  poly t;

  poly_basemul_montgomery(r, &a[0], &b[0]);
  for(i=1;i<vlen;i++) {
    poly_basemul_montgomery(&t, &a[i], &b[i]);
    poly_add(r, r, &t);
    poly_reduce(r);
  }
  
}

/*************************************************
* Name:        polyvec_reduce
*
* Description: Applies Barrett reduction to each coefficient
*              of each element of a vector of polynomials;
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
*              - int vlen: vector dimension 
**************************************************/
void polyvec_reduce(poly *r, int vlen)
{
  int i;
  for(i=0;i<vlen;i++)
    poly_reduce(&r[i]);
}

/*************************************************
* Name:        polyvec_reduce_mont
*
* Description: Applies Montgomery reduction to each coefficient
*              of each element of a vector of polynomials;
*              for details of the Montgomery reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
*              - int vlen: vector dimension 
**************************************************/
void polyvec_reduce_mont(poly *r, int vlen)
{
  int i;
  for(i=0;i<vlen;i++)
    poly_reduce_mont(&r[i]);
}



/*************************************************
* Name:        polyvec_add
*
* Description: Add vectors of polynomials
*
* Arguments: - polyvec *r: pointer to output vector of polynomials
*            - const poly *a: pointer to first input vector of polynomials
*            - const poly *b: pointer to second input vector of polynomials
*            - int vlen: vector dimension 
**************************************************/
void polyvec_add(poly *r, const poly *a, const poly *b, int vlen)
{
  int i;
  for(i=0;i<vlen;i++)
    poly_add(&r[i], &a[i], &b[i]);
}
