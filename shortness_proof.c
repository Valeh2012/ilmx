#include "string.h"
#include "shortness_proof.h"
#include "symmetric.h"
#include "ntt.h"
#include "randombytes.h"
#include "reduce.h"
#include "util.h"

/*************************************************
* Name:        shortness_prove
*
* Description: Generate Zero-Knowledge proof of shortness of rlwe encryption parameters
*
* Arguments:   - shortness_proof *p: pointer to proof
*              - const uint8_t rho[PQMX_SYMBYTES]: pointer to seed buffer
*              - poly *msg: secret message s in NTT domain which is solution to As=u equation
*              - const poly pk[2]: RLWE encryption public key
*              - poly *A: polynomial matrix A in As=u
*              - poly *u: polynomial vector u in As=u
**************************************************/
void shortness_prove(shortness_proof *p, const uint8_t rho[PQMX_SYMBYTES], poly *msg, const poly pk[2], poly *A, poly *u)
{
    unsigned int i,j, nonce=0;
    commkey ck;
    commrnd r;
    xof_state state;

    uint8_t buf[2*PQMX_SYMBYTES];
    uint8_t *thash = buf;
    uint8_t *chash = buf+PQMX_SYMBYTES;

    uint8_t seed[PQMX_SYMBYTES];
    randombytes(seed, PQMX_SYMBYTES);

    comm *t = (comm *) aligned_alloc(32, (PQMX_MU+PQMX_M)*sizeof(poly));
    memset(t,0,(PQMX_MU+PQMX_M)*sizeof(poly));

    shortness_expand_commkey(&ck, rho);
    shortness_commit(t, &r, msg, &ck); 
    poly g;
    poly_uniform(&g, seed, 0);

    for(j=0;j<4;j++)
        g.coeffs[j] = 0;

    poly_ntt(&g);

    poly_add(&t->tm[PQMX_M-3], &t->tm[PQMX_M-3], &g);
    poly_reduce(&t->tm[PQMX_M-3]);

    shake128_init(&state);
    shake128_absorb(&state, rho, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)t->t0, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)t->tm, (PQMX_M-2)*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);
    
    poly *alpha = (poly *) aligned_alloc(32, (PQMX_M-3)*sizeof(poly));
    poly *gamma = (poly *) aligned_alloc(32, (4*2*PQMX_NV)*sizeof(poly));
    poly *atgamma = (poly *) aligned_alloc(32, (PQMX_M-3)*sizeof(poly));
    commrnd y;
    y.s = (poly *) aligned_alloc(32, PQMX_LAMBDA*sizeof(poly));
    y.e = (poly *) aligned_alloc(32, PQMX_MU*sizeof(poly));
    y.em = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));

    poly *by = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));
    poly tmp;
    poly w[PQMX_MU];
    poly three;
    memset(&three, 0, sizeof(poly));
    three.coeffs[0] = 3;
    poly_ntt(&three);
    poly one;
    memset(&one, 0, sizeof(poly));
    one.coeffs[0] = 1;
    poly_ntt(&one);
    poly h,v,f,c;
    do{
        for(i=0;i<PQMX_LAMBDA;i++) {
            poly_uniform_delta(&y.s[i], seed, nonce++);
        }
        for(i=0;i<PQMX_MU;i++) {
            poly_uniform_delta(&y.e[i], seed, nonce++);
        }
        for(i=0;i<PQMX_M;i++) {
            poly_uniform_delta(&y.em[i], seed, nonce++);
        }

        polyvec_ntt(y.s, PQMX_LAMBDA);
        polyvec_ntt(y.e, PQMX_MU);
        polyvec_ntt(y.em, PQMX_M);

        for(j=0;j<PQMX_MU;j++){
            polyvec_basemul_acc_montgomery(&w[j],&ck.b0[j*PQMX_LAMBDA],y.s, PQMX_LAMBDA);
            polyvec_basemul_acc_montgomery(&tmp,&ck.bt[j*PQMX_M], y.em, PQMX_M);
            poly_add(&w[j],&w[j],&tmp);
            poly_tomont(&w[j]);
        }
        polyvec_add(w,w,y.e, PQMX_MU);
        polyvec_reduce(w, PQMX_MU);

        for(j=0;j<PQMX_M;j++){
            polyvec_basemul_acc_montgomery(&by[j], &ck.bm[j*PQMX_LAMBDA], y.s, PQMX_LAMBDA);
            poly_tomont(&by[j]);
        }
        polyvec_add(by, by, y.em, PQMX_M);
        polyvec_reduce(by, PQMX_M);

        // hash thash,w
        shake128_init(&state);
        shake128_absorb(&state, thash, PQMX_SYMBYTES);
        shake128_absorb(&state, (uint8_t*)w, PQMX_MU*sizeof(poly));
        shake128_finalize(&state);
        shake128_squeeze(thash, PQMX_SYMBYTES, &state);
             
        for(j=0;j<PQMX_M-3;j++){
            poly_uniform(&alpha[j], thash, nonce++);
        }
        for(j=0;j<4*2*PQMX_NV;j++){
            poly_uniform(&gamma[j], thash, nonce++);
        }

        poly_add(&t->tm[PQMX_M-2], &t->tm[PQMX_M-2], &by[PQMX_M-1]);
        poly_reduce(&t->tm[PQMX_M-2]);
        for(j=0;j<PQMX_M-3;j++){
            poly_basemul_montgomery(&tmp, &by[j], &by[j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &msg[j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &three);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &alpha[j]);
            poly_tomont(&tmp);
            poly_sub(&t->tm[PQMX_M-2], &t->tm[PQMX_M-2], &tmp);
            poly_reduce(&t->tm[PQMX_M-2]);
        }
        for(j=0;j<PQMX_M-3;j++){
            poly_basemul_montgomery(&tmp, &msg[j], &msg[j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &three);
            poly_tomont(&tmp);
            poly_sub(&tmp, &tmp, &one);
            poly_basemul_montgomery(&tmp, &tmp, &by[j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &alpha[j]);
            poly_tomont(&tmp);
            poly_add(&t->tm[PQMX_M-1], &t->tm[PQMX_M-1], &tmp);
            poly_reduce(&t->tm[PQMX_M-1]);
        }

        v = by[PQMX_M-2];
        for(j=0;j<PQMX_M-3;j++){
            poly_basemul_montgomery(&tmp, &by[j], &by[j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &by[j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &alpha[j]);
            poly_tomont(&tmp);
            poly_add(&v, &v, &tmp);
            poly_reduce(&v);
        }

        unsigned int offset = 4*2*PQMX_NV+4*(PQMX_LAMBDA+PQMX_NV);
        
        for(i=0;i<4*2*PQMX_NV;i++){
            atgamma[i] = gamma[i];
            atgamma[offset+i] = gamma[i];
        }
        offset = 4*2*PQMX_NV;
        int64_t tr[4], tr2[4];
        unsigned int k,o;
        for(i=0;i<4;i++){
            for(j=0;j<PQMX_N/4;j++){
                memset(tr, 0, sizeof(tr));
                for(o=0;o<2*PQMX_NV;o++){
                    poly_basemul_acc(tr2, &A[o], &gamma[4*o]);
                    for(k=0;k<4;k++){
                        tr[k] = barrett_reduce(tr[k]+tr2[k]); 
                    }
                    poly_shift(&A[o]);
                }
                memcpy(&atgamma[offset+i].coeffs[4*j], tr, 4*sizeof(int64_t));
            }
        }

        for(o=0;o<2*PQMX_NV;o++){
            for(i=0;i<PQMX_N;i++){
               A[o].coeffs[i] *= -1;
            }
        }
        offset+=4;

        poly tmp2;
        for(o=0;o<PQMX_NV;o++){
            tmp = pk[0];
            tmp2 = pk[1];
            for(i=0;i<4;i++){
                for(j=0;j<PQMX_N/4;j++){
                    poly_basemul_acc(tr, &tmp, &gamma[8*o]);
                    poly_basemul_acc(tr2, &tmp2, &gamma[8*o+4]);
                    for(k=0;k<4;k++){
                        tr[k] = barrett_reduce(tr[k]+tr2[k]); 
                    }
                    memcpy(&atgamma[offset+i].coeffs[4*j], tr, 4*sizeof(int64_t));
                    poly_shift(&tmp);
                    poly_shift(&tmp2);
                }
            }
            offset+=4;  
        }

        memset(&f, 0, sizeof(poly));
        for(j=0;j<PQMX_M-3;j++){
            poly_basemul_montgomery(&tmp, &msg[j], &atgamma[j]);
            poly_tomont(&tmp);
            poly_add(&f, &f, &tmp);
            poly_reduce(&f);
        }
        
        memset(tr, 0, sizeof(tr));
        for(o=0;o<2*PQMX_NV;o++){
            poly_basemul_acc(tr2, &u[o], &gamma[4*o]);
            for(i=0;i<4;i++){
                tr[i] = barrett_reduce(tr[i]+tr2[i]);
            }
        }
        for(i=0;i<4;i++){
            tr[i] = fqmul(tr[i], 18014398509481984L);
            for(j=0;j<PQMX_N/4;j++)
                f.coeffs[4*j + i] -= tr[i];
        }

        poly_add(&h, &f, &g);
        poly_reduce(&h);

        poly v2;
        memset(&v2, 0, sizeof(poly));
        for(j=0;j<PQMX_M-3;j++){
            poly_basemul_montgomery(&tmp, &by[j], &atgamma[j]);
            poly_tomont(&tmp);
            poly_add(&v2, &v2, &tmp);
            poly_reduce(&v2);
        }
        poly_add(&v2, &v2, &by[PQMX_M-3]);
        poly_reduce(&v2);

        shake128_init(&state);
        shake128_absorb(&state, thash, PQMX_SYMBYTES);
        shake128_absorb(&state, (uint8_t*)&v, sizeof(poly));
        shake128_absorb(&state, (uint8_t*)&v2, sizeof(poly));
        shake128_absorb(&state, (uint8_t*)(t->tm+PQMX_M-2), 2*sizeof(poly));
        shake128_finalize(&state);
        shake128_squeeze(chash, PQMX_SYMBYTES, &state);

        // sample c
        poly_nonuniform(&c, 0x80, chash, 0);
        poly_ntt(&c);
        // z = y + c*r
        for(i=0;i<PQMX_MU;i++){
            poly_basemul_montgomery(&tmp, &c, &r.e[i]);
            poly_tomont(&tmp);
            poly_add(&y.e[i], &tmp, &y.e[i]);
            poly_reduce(&y.e[i]);
        }
        for(i=0;i<PQMX_LAMBDA;i++){
            poly_basemul_montgomery(&tmp, &c, &r.s[i]);
            poly_tomont(&tmp);
            poly_add(&y.s[i], &tmp, &y.s[i]);
            poly_reduce(&y.s[i]);
        }
        for(i=0;i<PQMX_M;i++){
            poly_basemul_montgomery(&tmp, &c, &r.em[i]);
            poly_tomont(&tmp);
            poly_add(&y.em[i], &tmp, &y.em[i]);
            poly_reduce(&y.em[i]);
        }
        
        // do reje ction sampling
        polyvec_invntt_tomont(y.s, PQMX_LAMBDA);
        polyvec_reduce_mont(y.s, PQMX_LAMBDA);
        polyvec_invntt_tomont(y.e, PQMX_MU);
        polyvec_reduce_mont(y.e, PQMX_MU);
        polyvec_invntt_tomont(y.em, PQMX_M);
        polyvec_reduce_mont(y.em, PQMX_M);

        if(polyvec_norm(y.s, 0, PQMX_LAMBDA) < PQMX_DELTA - PQMX_BETA){
            break;
        }
        if(polyvec_norm(y.e, 0, PQMX_MU) < PQMX_DELTA - PQMX_BETA){
            break;
        }
        if(polyvec_norm(y.em, 0, PQMX_M) < PQMX_DELTA - PQMX_BETA){
            break;
        }
    }while(1);

    p->z.e = y.e;
    p->z.em = y.em;
    p->z.s = y.s;
    p->h = h;
    memcpy(p->hash, chash, PQMX_SYMBYTES);
    p->t = t;

    free(by);
    free(alpha);
    free(r.e);
    free(r.em);
    free(r.s);
    free(ck.b0);
    free(ck.bm);
    free(ck.bt);
}

/*************************************************
* Name:        shortness_verify
*
* Description: Verify Zero-Knowledge proof of shortness of rlwe encryption parameters
*
* Arguments:   - shortness_proof *p: pointer to proof
*              - const uint8_t rho[PQMX_SYMBYTES]: pointer to seed buffer
*              - const poly pk[2]: RLWE encryption public key
*              - poly *A: polynomial matrix A in As=u
*              - poly *u: polynomial vector u in As=u
*
* Returns:     - 0: if verified successfully
*              - 1: otherwise
**************************************************/
int shortness_verify(shortness_proof *p, const uint8_t rho[PQMX_SYMBYTES], const poly pk[2], poly *A, poly *u)
{
    if(polyvec_norm(p->z.s, 0, PQMX_LAMBDA) >= PQMX_DELTA - PQMX_BETA){
        return 1;
    }
    if(polyvec_norm(p->z.e, 0, PQMX_MU) >= PQMX_DELTA - PQMX_BETA){
        return 1;
    }
    if(polyvec_norm(p->z.em, 0, PQMX_M) >= PQMX_DELTA - PQMX_BETA){
        return 1;
    }
    
    polyvec_ntt(p->z.s, PQMX_LAMBDA);
    polyvec_ntt(p->z.e, PQMX_MU);
    polyvec_ntt(p->z.em, PQMX_M);

    unsigned int i,j, nonce=0;
    commkey ck;
    xof_state state;

    uint8_t buf[2*PQMX_SYMBYTES];
    uint8_t *thash = buf;
    uint8_t *chash = buf+PQMX_SYMBYTES;

    shortness_expand_commkey(&ck, rho);
    comm *t = p->t;

    shake128_init(&state);
    shake128_absorb(&state, rho, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)t->t0, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)t->tm, (PQMX_M-2)*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);
    
    poly c;
    poly_nonuniform(&c, 0x80, p->hash, 0);
    poly_ntt(&c);

    poly tmp;
    poly w[PQMX_MU];
    for(j=0;j<PQMX_MU;j++){
      polyvec_basemul_acc_montgomery(&w[j],&ck.b0[j*PQMX_LAMBDA],p->z.s, PQMX_LAMBDA);
      polyvec_basemul_acc_montgomery(&tmp,&ck.bt[j*PQMX_M], p->z.em, PQMX_M);
      poly_add(&w[j],&w[j],&tmp);
      poly_tomont(&w[j]);
    }
    polyvec_add(w,w,p->z.e, PQMX_MU);
    polyvec_reduce(w, PQMX_MU);

    for(i=0;i<PQMX_MU;i++){
        poly_basemul_montgomery(&tmp, &c, &t->t0[i]);
        poly_tomont(&tmp);
        poly_sub(&w[i], &w[i], &tmp);
        poly_reduce(&w[i]);
    }

    poly *f = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));
    poly *bz = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));
    for(j=0;j<PQMX_M;j++){
        polyvec_basemul_acc_montgomery(&bz[j], &ck.bm[j*PQMX_LAMBDA], p->z.s, PQMX_LAMBDA);
        poly_tomont(&bz[j]);
    }
    polyvec_add(bz,bz,p->z.em, PQMX_M);
    polyvec_reduce(bz, PQMX_M);

    for(j=0;j<PQMX_M;j++){
        poly_basemul_montgomery(&tmp, &c, &t->tm[j]);
        poly_tomont(&tmp);
        poly_sub(&f[j], &bz[j], &tmp);
        poly_reduce(&f[j]);
    }

    poly *alpha = (poly *) aligned_alloc(32, (PQMX_M-3)*sizeof(poly));
    poly *gamma = (poly *) aligned_alloc(32, (4*2*PQMX_NV)*sizeof(poly));
    poly *atgamma = (poly *) aligned_alloc(32, (PQMX_M-3)*sizeof(poly));
    // hash thash,w
    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)w, PQMX_MU*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);
    
    
    nonce += PQMX_M + PQMX_MU + PQMX_LAMBDA;

    for(j=0;j<PQMX_M-3;j++){
        poly_uniform(&alpha[j], thash, nonce++);
    }
    for(j=0;j<4*2*PQMX_NV;j++){
        poly_uniform(&gamma[j], thash, nonce++);
    }

    poly v, v2, tmp2;
    memset(&v, 0, sizeof(poly));
    memset(&v2, 0, sizeof(poly));
    for(j=0;j<PQMX_M-3;j++){
        poly_add(&tmp, &f[j], &c);
        poly_reduce(&tmp);
        poly_basemul_montgomery(&tmp2,&tmp,&f[j]);
        poly_tomont(&tmp2);
        poly_sub(&tmp, &f[j], &c);
        poly_reduce(&tmp);
        poly_basemul_montgomery(&tmp2, &tmp2, &tmp);
        poly_tomont(&tmp2);
        poly_basemul_montgomery(&tmp2, &tmp2, &alpha[j]);
        poly_tomont(&tmp2);
        poly_add(&v, &v, &tmp2);
        poly_reduce(&v);
    }
    poly_basemul_montgomery(&tmp, &f[PQMX_M-1], &c);
    poly_tomont(&tmp);
    poly_add(&v, &v, &tmp);
    poly_reduce(&v);
    poly_add(&v, &v, &f[PQMX_M-2]);
    poly_reduce(&v);

    unsigned int offset = 4*2*PQMX_NV+4*(PQMX_LAMBDA+PQMX_NV);
        
    for(i=0;i<4*2*PQMX_NV;i++){
        atgamma[i] = gamma[i];
        atgamma[offset+i] = gamma[i];
    }
    offset = 4*2*PQMX_NV;
    int64_t tr[4], tr2[4];
    unsigned int k,o;
    for(i=0;i<4;i++){
        for(j=0;j<PQMX_N/4;j++){
            memset(tr, 0, sizeof(tr));
            for(o=0;o<2*PQMX_NV;o++){
                poly_basemul_acc(tr2, &A[o], &gamma[4*o]);
                for(k=0;k<4;k++){
                    tr[k] = barrett_reduce(tr[k]+tr2[k]); 
                }
                poly_shift(&A[o]);
            }
            memcpy(&atgamma[offset+i].coeffs[4*j], tr, 4*sizeof(int64_t));
        }
    }
    
    for(o=0;o<2*PQMX_NV;o++){
        for(i=0;i<PQMX_N;i++){
            A[o].coeffs[i] *= -1;
        }
    }

    offset+=4;

    for(o=0;o<PQMX_NV;o++){
        tmp = pk[0];
        tmp2 = pk[1];
        for(i=0;i<4;i++){
            for(j=0;j<PQMX_N/4;j++){
                poly_basemul_acc(tr, &tmp, &gamma[8*o]);
                poly_basemul_acc(tr2, &tmp2, &gamma[8*o+4]);
                for(k=0;k<4;k++){
                    tr[k] = barrett_reduce(tr[k]+tr2[k]); 
                }
                memcpy(&atgamma[offset+i].coeffs[4*j], tr, 4*sizeof(int64_t));
                poly_shift(&tmp);
                poly_shift(&tmp2);
            }
        }
        offset+=4;
    }
    
    poly tf;
    memset(&tf, 0, sizeof(poly));
    
    for(j=0;j<PQMX_M-3;j++){
        poly_basemul_montgomery(&tmp, &t->tm[j], &atgamma[j]);
        poly_tomont(&tmp);
        poly_add(&tf, &tf, &tmp);
        poly_reduce(&tf);
    }
    memset(tr, 0, sizeof(tr));
    for(o=0;o<2*PQMX_NV;o++){
        poly_basemul_acc(tr2, &u[o], &gamma[4*o]);
        for(i=0;i<4;i++){
            tr[i] = barrett_reduce(tr[i]+tr2[i]);
        }
    }
    for(i=0;i<4;i++){
        tr[i] = fqmul(tr[i], 18014398509481984L);
        for(j=0;j<PQMX_N/4;j++)
            tf.coeffs[4*j + i] -= tr[i];
    }

    memset(&v2, 0, sizeof(poly));
    for(j=0;j<PQMX_M-3;j++){
        poly_basemul_montgomery(&tmp, &bz[j], &atgamma[j]);
        poly_tomont(&tmp);
        poly_add(&v2, &v2, &tmp);
        poly_reduce(&v2);
    }
    poly_add(&v2, &v2, &bz[PQMX_M-3]);
    poly_reduce(&v2);

    poly_add(&tmp, &tf, &t->tm[PQMX_M-3]);
    poly_reduce(&tmp);
    poly_sub(&tmp, &tmp, &p->h);
    poly_reduce(&tmp);
    poly_basemul_montgomery(&tmp, &tmp, &c);
    poly_tomont(&tmp);

    poly_sub(&v2, &v2, &tmp);
    poly_reduce(&v2);


    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)&v, sizeof(poly));
    shake128_absorb(&state, (uint8_t*)&v2, sizeof(poly));
    shake128_absorb(&state, (uint8_t*)(t->tm+PQMX_M-2), 2*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(chash, PQMX_SYMBYTES, &state);

    free(f);
    free(bz);
    free(alpha);
    free(gamma);
    free(atgamma);
    free(ck.b0);
    free(ck.bm);
    free(ck.bt);

    polyvec_invntt_tomont(p->z.s, PQMX_LAMBDA);
    polyvec_invntt_tomont(p->z.e, PQMX_MU);
    polyvec_invntt_tomont(p->z.em, PQMX_M);

    polyvec_reduce_mont(p->z.s, PQMX_LAMBDA);
    polyvec_reduce_mont(p->z.e, PQMX_MU);
    polyvec_reduce_mont(p->z.em, PQMX_M);

    uint8_t r = 0;
    for(j=0;j<PQMX_SYMBYTES;j++){
        r |= chash[j] ^ p->hash[j];
    }

    return (-(int)r)>>31;
}