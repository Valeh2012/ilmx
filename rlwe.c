#include "rlwe.h"
#include "randombytes.h"
#include "symmetric.h"
#include "math.h"
#include "reduce.h"

/*************************************************
* Name:        rlwe_genkey
*
* Description: Generate a RLWE public and private keypair
*
* Arguments:   - poly pk[2]: public key polynomials
*              - poly sk: secret key polynomial
**************************************************/
void rlwe_genkey(poly pk[2], poly * sk)
{
    uint8_t seed[PQMX_SYMBYTES];
    unsigned int nonce=0;

    randombytes(seed, PQMX_SYMBYTES);
    hash_g(seed, seed, PQMX_SYMBYTES);  // Because kyber does so

    poly_uniform(&pk[0], seed, nonce); 
    poly_ntt(&pk[0]);

    poly e;
    poly_nonuniform(&e, 0x80, seed, nonce++);
    poly_ntt(&e);

    poly_nonuniform(sk, 0x80, seed, nonce++);
    poly_ntt(sk);

    poly_basemul_montgomery(&pk[1], &pk[0], sk);
    poly_tomont(&pk[1]);
    poly_add(&pk[1], &pk[1], &e);
    poly_reduce(&pk[1]);
    
}

/*************************************************
* Name:        rlwe_enc
*
* Description: RLWE Encryption of a polynomial
*
* Arguments:   - rlwecp *cp: pointer to output ciphertext
*              - const poly *msg: pointer to message polynomial
*              - const poly pk[2]: public key polynomials
*              - rlwernd *rnd: pointer to output rlwe randomnesses
**************************************************/
void rlwe_enc(rlwecp* cp, const poly *msg, const poly pk[2], rlwernd *rnd)
{
    uint8_t seed[PQMX_SYMBYTES];
    unsigned int nonce=0;

    randombytes(seed, PQMX_SYMBYTES);

    poly r,e1,e2;
    poly_nonuniform(&r, 0x80, seed, nonce++);
    poly_nonuniform(&e1, 0x80, seed, nonce++);
    poly_nonuniform(&e2, 0x80, seed, nonce++);

    rnd->r = r;
    rnd->eu = e1;
    rnd->ev = e2;

    poly_ntt(&r);
    poly_ntt(&e1);
    poly_ntt(&e2);

    poly_basemul_montgomery(&cp->u, &pk[0], &r);
    poly_tomont(&cp->u);
    poly_add(&cp->u, &cp->u, &e1);
    poly_reduce(&cp->u);
    
    poly m = *msg;
    poly_ntt(&m);

    poly_basemul_montgomery(&cp->v, &pk[1], &r);
    poly_tomont(&cp->v);
    poly_add(&cp->v, &cp->v, &e2);
    poly_reduce(&cp->v);
    poly_add(&cp->v, &cp->v, &m);
    poly_reduce(&cp->v);
}

/*************************************************
* Name:        rlwe_dec
*
* Description: RLWE Decryption of a polynomial
*
* Arguments:   - poly *msg: pointer to output message polynomial
*              - const *poly sk: pointer to private key polynomial
*              - rlwecp *cp: pointer to input ciphertext
**************************************************/
void rlwe_dec(poly *msg, const poly * sk, const rlwecp* cp)
{
    unsigned int j;
    poly res; 

    poly_basemul_montgomery(&res, sk, &cp->u);
    poly_tomont(&res);
    poly_sub(&res, &cp->v, &res);
    poly_reduce(&res);
    poly_invntt_tomont(&res);
    
    int64_t t;
    for(j=0;j<PQMX_N;j++){
        t = labs(montgomery_reduce((__int128 )res.coeffs[j])); // think better way
        msg->coeffs[j] = ( t < PQMX_Q/4 ) ? 0 : 1;
    }
    
}