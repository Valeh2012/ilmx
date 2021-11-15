#include <string.h>
#include "mxproof.h"
#include "symmetric.h"
#include "randombytes.h"
#include "shortness_proof.h"
#include "util.h"

/*************************************************
* Name:        mx_proof
*
* Description: Generate Zero-Knowledge proof of a shuffle
*
* Arguments:   - proof *p: pointer to output proof
*              - const uint8_t rho[2*PQMX_SYMBYTES]: pointer to seed buffer
*              - const poly rlwepk[2]: RLWE encryption public key
*              - const rlwernd rnd[PQMX_NV]: RLWE encryption randomnesses
*              - const rlwecp out[PQMX_NV]: mix-net output ciphertexts
*              - const rlwecp reenc[PQMX_NV]: RLWE 0-encryption ciphertexts
*              - const int64_t pi[PQMX_NV]: permutation vector used in mix-net
**************************************************/
void mx_proof(proof *p, const uint8_t rho[2*PQMX_SYMBYTES], const poly rlwepk[2], const rlwernd rnd[PQMX_NV], const rlwecp out[PQMX_NV], const rlwecp reenc[PQMX_NV], const int64_t pi[PQMX_NV])
{
    unsigned int i,j,k, nonce=0;
    uint8_t symbuf[3*PQMX_SYMBYTES];
    uint8_t *seed = symbuf;
    uint8_t *thash = symbuf+PQMX_SYMBYTES;
    uint8_t *chash = symbuf+2*PQMX_SYMBYTES;

    xof_state state;
    randombytes(seed, PQMX_SYMBYTES);

    poly *m = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));
    memset(m, 0, PQMX_M*sizeof(poly));
    for(j=0;j<PQMX_NV;j++){
        m[j] = reenc[j].u;
        m[PQMX_NV + j] = reenc[j].v;
        m[2 * PQMX_NV + j].coeffs[0] = pi[j];
        poly_ntt(&m[2 * PQMX_NV + j]);
    }

    commkey ck;
    commrnd r;
    comm *t = (comm *) aligned_alloc(32, (PQMX_MU+PQMX_M)*sizeof(poly));
    memset(t, 0, (PQMX_MU+PQMX_M)*sizeof(poly));
    
    expand_commkey(&ck, rho);
    commit(t, &r, m, &ck);

    // shortness proof
    polyvec_invntt_tomont(r.s,  PQMX_LAMBDA);
    polyvec_invntt_tomont(r.e,  PQMX_MU);
    polyvec_invntt_tomont(r.em, PQMX_M);

    polyvec_reduce_mont(r.s,  PQMX_LAMBDA);
    polyvec_reduce_mont(r.e,  PQMX_MU);
    polyvec_reduce_mont(r.em, PQMX_M);

    poly *m2  = (poly *) aligned_alloc(32, (4*(PQMX_LAMBDA+5*PQMX_NV)+3)*sizeof(poly));
    memset(m2, 0, (4*(PQMX_LAMBDA+5*PQMX_NV)+3 )*sizeof(poly));

    for(k=0;k<2*PQMX_NV;k++){
        for(i=0;i<4;i++){
            for(j=0;j<PQMX_N/4;j++){
                m2[4*k+i].coeffs[4*j] = r.em[k].coeffs[i*PQMX_N/4 +j];
            }
        }
        }
        for(k=0;k<PQMX_LAMBDA;k++){
        for(i=0;i<4;i++){
            for(j=0;j<PQMX_N/4;j++){
                m2[4*(2*PQMX_NV+k)+i].coeffs[4*j] = r.s[k].coeffs[i*PQMX_N/4 +j];
            }
        }
    }

    for(k=0;k<PQMX_NV;k++){
        for(i=0;i<4;i++){
            for(j=0;j<PQMX_N/4;j++){
                m2[4*(2*PQMX_NV+PQMX_LAMBDA+k)+i].coeffs[4*j] = rnd[k].r.coeffs[i*PQMX_N/4 +j];
            }
        }
        for(i=0;i<4;i++){
            for(j=0;j<PQMX_N/4;j++){
                m2[4*(2*PQMX_NV+PQMX_LAMBDA+PQMX_NV+k)+i].coeffs[4*j] = rnd[k].eu.coeffs[i*PQMX_N/4 +j];
            }
        }
        for(i=0;i<4;i++){
            for(j=0;j<PQMX_N/4;j++){
                m2[4*(2*PQMX_NV+PQMX_LAMBDA+2*PQMX_NV+k)+i].coeffs[4*j] = rnd[k].ev.coeffs[i*PQMX_N/4 +j];
            }
        }
    }

    shortness_proof shortness_p;
    for(j=0;j<2*PQMX_NV;j++){
        polyvec_invntt_tomont(&ck.bm[j], PQMX_LAMBDA);
        polyvec_reduce_mont(&ck.bm[j], PQMX_LAMBDA);
    }

    polyvec_invntt_tomont(t->tm,2*PQMX_NV);
    polyvec_reduce_mont(t->tm,2*PQMX_NV);

    shortness_prove(&shortness_p, rho+PQMX_SYMBYTES, m2, rlwepk, ck.bm, t->tm);
    
    free(m2);
    for(j=0;j<2*PQMX_NV;j++){
        polyvec_ntt(&ck.bm[j], PQMX_LAMBDA);
    }

    polyvec_ntt(r.s,  PQMX_LAMBDA);
    polyvec_ntt(r.e,  PQMX_MU);
    polyvec_ntt(r.em, PQMX_M);

    polyvec_ntt(t->tm,2*PQMX_NV);

    //randombytes(seed, PQMX_SYMBYTES);

    shake128_init(&state);
    shake128_absorb(&state, rho, 2*PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)t->t0, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)t->tm, 3*PQMX_NV*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)shortness_p.t->t0, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)shortness_p.t->tm, (4*PQMX_LAMBDA+20*PQMX_NV+3)*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)&shortness_p.h, sizeof(poly));
    shake128_absorb(&state, shortness_p.hash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)shortness_p.z.e, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)shortness_p.z.em, (4*PQMX_LAMBDA+20*PQMX_NV+3)*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)shortness_p.z.s, PQMX_LAMBDA*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);
    

    poly alpha;
    poly_uniform_alpha(&alpha, thash, nonce++);
    poly_ntt(&alpha);

    poly *alphapowpis = (poly *) aligned_alloc(32, PQMX_NV*sizeof(poly));
    memset(alphapowpis,0,PQMX_NV*sizeof(poly));
    alphapowpis[0] = alpha;
    for(j=1;j<PQMX_NV;j++){
        poly_basemul_montgomery(&alphapowpis[j], &alphapowpis[j-1], &alpha);
        poly_tomont(&alphapowpis[j]);
    }

    poly tmp;
    for(j=0;j<PQMX_NV;j++){
        m[3*PQMX_NV+j] = alphapowpis[pi[j]-1];
    }
    polyvec_add(&t->tm[3*PQMX_NV],&t->tm[3*PQMX_NV], &m[3*PQMX_NV], PQMX_NV);
    poly_reduce(&t->tm[3*PQMX_NV+j]);

    for(j=0; j<PQMX_NV;j++){
        poly_basemul_montgomery(&tmp, &m[3*PQMX_NV+j], &reenc[j].u);
        poly_tomont(&tmp);
        m[4*PQMX_NV+j] = tmp;
        poly_add(&t->tm[4*PQMX_NV+j],&t->tm[4*PQMX_NV+j],&tmp);
        poly_reduce(&t->tm[4*PQMX_NV+j]);

        poly_basemul_montgomery(&tmp, &m[3*PQMX_NV+j], &reenc[j].v);
        poly_tomont(&tmp);
        m[5*PQMX_NV+j] = tmp;
        poly_add(&t->tm[5*PQMX_NV+j],&t->tm[5*PQMX_NV+j],&tmp);
        poly_reduce(&t->tm[5*PQMX_NV+j]);
    }

    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)(t->tm+3*PQMX_NV), 3*PQMX_NV*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    poly beta, gamma;
    poly_uniform_alpha(&beta, thash, nonce++);
    poly_ntt(&beta);
    poly_uniform_alpha(&gamma, thash, nonce++);
    poly_ntt(&gamma);

    poly *P = (poly *) aligned_alloc(32, (PQMX_NV+1)*sizeof(poly));    
    memset(P, 0, (PQMX_NV+1)*sizeof(poly));
    P[0].coeffs[0]=1;
    poly_ntt(&P[0]);
    
    for(j=0;j<PQMX_NV;j++){
        poly_basemul_montgomery(&tmp, &m[2*PQMX_NV + j], &beta);
        poly_tomont(&tmp);
        poly_add(&tmp, &tmp, &m[3*PQMX_NV+j]);
        poly_reduce(&tmp);
        poly_sub(&tmp, &tmp, &gamma);
        poly_reduce(&tmp);
        m[6*PQMX_NV+j] = tmp;
        poly_add(&t->tm[6*PQMX_NV+j], &t->tm[6*PQMX_NV+j], &tmp);
        poly_reduce(&t->tm[6*PQMX_NV+j]);

        m[7*PQMX_NV+j] = P[j];
        poly_add(&t->tm[7*PQMX_NV+j], &t->tm[7*PQMX_NV+j], &P[j]);
        poly_reduce(&t->tm[7*PQMX_NV+j]);

        poly_basemul_montgomery(&P[j+1], &P[j], &tmp);
        poly_tomont(&P[j+1]);
        m[8*PQMX_NV+j] = P[j+1];
        poly_add(&t->tm[8*PQMX_NV+j], &t->tm[8*PQMX_NV+j], &P[j+1]);
        poly_reduce(&t->tm[8*PQMX_NV+j]);
    }

    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)(t->tm+6*PQMX_NV), 3*PQMX_NV*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    // generate uniform y
    commrnd y;
    y.s = (poly *) aligned_alloc(32, PQMX_LAMBDA*sizeof(poly));
    y.e = (poly *) aligned_alloc(32, PQMX_MU*sizeof(poly));
    y.em = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));

    poly *by = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));
    poly *epsilon = (poly *) aligned_alloc(32, (4*PQMX_NV+4)*sizeof(poly));
    poly v[4];
    poly w[PQMX_MU];
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
        // w = By

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

        
        for(j=0;j<4*PQMX_NV+4;j++){
            poly_uniform(&epsilon[j], thash, nonce++);
            poly_ntt(&epsilon[j]);
        }

        memset(v, 0, 4*sizeof(poly));

        for(j=0;j<PQMX_NV;j++){
            poly_basemul_montgomery(&tmp, &by[2*PQMX_NV+j], &beta);
            poly_tomont(&tmp);
            poly_add(&tmp, &tmp, &by[3*PQMX_NV+j]);
            poly_reduce(&tmp);
            poly_sub(&tmp, &tmp, &by[6*PQMX_NV+j]);
            poly_reduce(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &epsilon[j]);
            poly_tomont(&tmp);
            poly_add(&v[0], &v[0], &tmp);
            poly_reduce(&v[0]);
        }
        poly vtmp;
        v[1] = by[PQMX_M-1];
        for(j=0;j<PQMX_NV;j++){
            poly_basemul_montgomery(&tmp, &by[6*PQMX_NV+j], &by[7*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &epsilon[PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_add(&v[1], &v[1], &tmp);
            poly_reduce(&v[1]);

            poly_basemul_montgomery(&tmp, &by[3*PQMX_NV+j], &by[j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &epsilon[2*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_add(&v[1], &v[1], &tmp);
            poly_reduce(&v[1]);

            poly_basemul_montgomery(&tmp, &by[3*PQMX_NV+j], &by[PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_basemul_montgomery(&tmp, &tmp, &epsilon[3*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_add(&v[1], &v[1], &tmp);
            poly_reduce(&v[1]);
        }

        for(j=0;j<PQMX_NV;j++){
            poly_basemul_montgomery(&tmp, &by[6*PQMX_NV+j], &m[7*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_sub(&vtmp, &by[8*PQMX_NV+j], &tmp);
            poly_basemul_montgomery(&tmp, &by[7*PQMX_NV+j], &m[6*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_sub(&vtmp, &vtmp, &tmp);
            poly_basemul_montgomery(&vtmp, &vtmp, &epsilon[PQMX_NV+j]);
            poly_tomont(&vtmp);
            poly_add(&t->tm[PQMX_M-1], &t->tm[PQMX_M-1], &vtmp);
            poly_reduce(&t->tm[PQMX_M-1]);

            poly_basemul_montgomery(&tmp, &by[j], &m[3*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_sub(&vtmp, &by[4*PQMX_NV+j], &tmp);
            poly_basemul_montgomery(&tmp, &by[3*PQMX_NV+j], &m[j]);
            poly_tomont(&tmp);
            poly_sub(&vtmp, &vtmp, &tmp);
            poly_basemul_montgomery(&vtmp, &vtmp, &epsilon[2*PQMX_NV+j]);
            poly_tomont(&vtmp);
            poly_add(&t->tm[PQMX_M-1], &t->tm[PQMX_M-1], &vtmp);
            poly_reduce(&t->tm[PQMX_M-1]);

            poly_basemul_montgomery(&tmp, &by[PQMX_NV+j], &m[3*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_sub(&vtmp, &by[5*PQMX_NV+j], &tmp);
            poly_basemul_montgomery(&tmp, &by[3*PQMX_NV+j], &m[PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_sub(&vtmp, &vtmp, &tmp);
            poly_basemul_montgomery(&vtmp, &vtmp, &epsilon[3*PQMX_NV+j]);
            poly_tomont(&vtmp);
            poly_add(&t->tm[PQMX_M-1], &t->tm[PQMX_M-1], &vtmp);
            poly_reduce(&t->tm[PQMX_M-1]);
        }

        memset(&vtmp, 0 , sizeof(poly));
        for(j=0;j<PQMX_NV;j++){
            poly_basemul_montgomery(&tmp, &out[j].u, &by[3*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_add(&vtmp, &vtmp, &tmp);
            poly_reduce(&vtmp);
            poly_sub(&vtmp, &vtmp, &by[4*PQMX_NV+j]);
            poly_reduce(&vtmp);
        }
        poly_basemul_montgomery(&vtmp, &vtmp, &epsilon[4*PQMX_NV]);
        poly_tomont(&vtmp);
        poly_add(&v[2], &v[2], &vtmp);

        memset(&vtmp, 0 , sizeof(poly));
        for(j=0;j<PQMX_NV;j++){
            poly_basemul_montgomery(&tmp, &out[j].v, &by[3*PQMX_NV+j]);
            poly_tomont(&tmp);
            poly_add(&vtmp, &vtmp, &tmp);
            poly_reduce(&vtmp);
            poly_sub(&vtmp, &vtmp, &by[5*PQMX_NV+j]);
            poly_reduce(&vtmp);
        }
        poly_basemul_montgomery(&vtmp, &vtmp, &epsilon[4*PQMX_NV+1]);
        poly_tomont(&vtmp);
        poly_add(&v[2], &v[2], &vtmp);
        poly_reduce(&v[2]);

        poly_basemul_montgomery(&tmp, &epsilon[4*PQMX_NV+2], &by[9*PQMX_NV-1]);
        poly_tomont(&tmp);
        poly_add(&v[3], &v[3], &tmp);
        poly_basemul_montgomery(&tmp, &epsilon[4*PQMX_NV+3], &by[7*PQMX_NV]);
        poly_tomont(&tmp);
        poly_add(&v[3], &v[3], &tmp);
        poly_reduce(&v[3]);

        shake128_init(&state);
        shake128_absorb(&state, thash, PQMX_SYMBYTES);
        shake128_absorb(&state, (uint8_t*)v, 4*sizeof(poly));
        shake128_absorb(&state, (uint8_t*)(t->tm+9*PQMX_NV), sizeof(poly));
        shake128_finalize(&state);
        shake128_squeeze(chash, PQMX_SYMBYTES, &state);
        
        // sample c
        poly c;
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
    
    //rejection sampling -- means check norm
    p->z.e = y.e;
    p->z.em = y.em;
    p->z.s = y.s;
    p->sp = (void *) &shortness_p;
    p->t = t;
    memcpy(p->hash, chash, PQMX_SYMBYTES);

    free(m);
    free(P);
    free(alphapowpis);
    free(by);
    free(epsilon);
    free(r.e);
    free(r.em);
    free(r.s);
    free(ck.b0);
    free(ck.bm);
    free(ck.bt);
}

/*************************************************
* Name:        mx_proof
*
* Description: Verify Zero-Knowledge proof of a shuffle
*
* Arguments:   - proof *p: pointer to proof
*              - const uint8_t rho[2*PQMX_SYMBYTES]: pointer to seed buffer
*              - const poly rlwepk[2]: RLWE encryption public key
*              - const rlwecp in[PQMX_NV]: mix-net input ciphertexts
*              - const rlwecp out[PQMX_NV]: mix-net output ciphertexts
*
* Returns:     - 0: if verified successfully
*              - 1: otherwise
**************************************************/
int mx_proof_verify(const proof *p, const uint8_t rho[2*PQMX_SYMBYTES], poly rlwepk[2], const rlwecp in[PQMX_NV], const rlwecp out[PQMX_NV])
{
    unsigned int i,j, nonce=0;
    uint8_t symbuf[2*PQMX_SYMBYTES];
    uint8_t *thash = symbuf;
    uint8_t *chash = symbuf+PQMX_SYMBYTES;
    xof_state state;

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

    commkey ck;
    comm *t = p->t;
    expand_commkey(&ck, rho);

    shortness_proof *sp = (shortness_proof *) p->sp;
    for(j=0;j<2*PQMX_NV;j++){
        polyvec_invntt_tomont(&ck.bm[j], PQMX_LAMBDA);
        polyvec_reduce_mont(&ck.bm[j], PQMX_LAMBDA);
    }

    polyvec_invntt_tomont(t->tm,2*PQMX_NV);
    polyvec_reduce_mont(t->tm,2*PQMX_NV);

    int res=0; 
    res = shortness_verify(sp, rho + PQMX_SYMBYTES, rlwepk, ck.bm, t->tm);

    for(j=0;j<2*PQMX_NV;j++){
        polyvec_ntt(&ck.bm[j], PQMX_LAMBDA);
    }
    polyvec_ntt(t->tm, 2*PQMX_NV);

    if(res){
        return 1;
    }

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
    for(j=0;j<PQMX_M;j++){
        polyvec_basemul_acc_montgomery(&f[j], &ck.bm[j*PQMX_LAMBDA], p->z.s, PQMX_LAMBDA);
        poly_tomont(&f[j]);
        poly_basemul_montgomery(&tmp, &c, &t->tm[j]);
        poly_tomont(&tmp);
        poly_sub(&f[j], &f[j], &tmp);
        poly_reduce(&f[j]);
    }
    polyvec_add(f,f,p->z.em, PQMX_M);
    polyvec_reduce(f, PQMX_M);

    // hash t,t0,w

    shake128_init(&state);
    shake128_absorb(&state, rho, 2*PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)t->t0, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)t->tm, 3*PQMX_NV*sizeof(poly));
    shake128_absorb(&state, (uint8_t*) sp->t->t0, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*) sp->t->tm, (4*PQMX_LAMBDA+20*PQMX_NV+3)*sizeof(poly));
    shake128_absorb(&state, (uint8_t*) &sp->h, sizeof(poly));
    shake128_absorb(&state, sp->hash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*) sp->z.e, PQMX_MU*sizeof(poly));
    shake128_absorb(&state, (uint8_t*) sp->z.em, (4*PQMX_LAMBDA+20*PQMX_NV+3)*sizeof(poly));
    shake128_absorb(&state, (uint8_t*) sp->z.s, PQMX_LAMBDA*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    poly alpha;
    poly_uniform_alpha(&alpha, thash, nonce++);
    poly_ntt(&alpha);
    poly *alphapowpis = (poly *) aligned_alloc(32, PQMX_NV*sizeof(poly));
    memset(alphapowpis,0,PQMX_NV*sizeof(poly));
    alphapowpis[0] = alpha;
    for(j=1;j<PQMX_NV;j++){
        poly_basemul_montgomery(&alphapowpis[j], &alphapowpis[j-1], &alpha);
        poly_tomont(&alphapowpis[j]);
    }

    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)(t->tm+3*PQMX_NV), 3*PQMX_NV*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    poly beta, gamma;
    poly_uniform_alpha(&beta, thash, nonce++);
    poly_ntt(&beta);
    poly_uniform_alpha(&gamma, thash, nonce++);
    poly_ntt(&gamma);

    poly P;
    memset(&P, 0, sizeof(poly));
    P.coeffs[0] = 1;
    poly_ntt(&P);
    for(j=0;j<PQMX_NV;j++){
        memset(&tmp, 0, sizeof(poly));
        tmp.coeffs[0] = j+1;
        poly_ntt(&tmp);
        poly_basemul_montgomery(&tmp, &tmp, &beta);
        poly_tomont(&tmp);
        poly_add(&tmp, &tmp, &alphapowpis[j]);
        poly_reduce(&tmp);
        poly_sub(&tmp, &tmp, &gamma);
        poly_reduce(&tmp);
        poly_basemul_montgomery(&P, &P, &tmp);
        poly_tomont(&P);
    }
    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)(t->tm+6*PQMX_NV), 3*PQMX_NV*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)w, PQMX_MU*sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    nonce+= PQMX_MU + PQMX_LAMBDA + PQMX_M;

    poly *epsilon = (poly *) aligned_alloc(32, (4*PQMX_NV+4)*sizeof(poly));
    for(j=0;j<4*PQMX_NV+4;j++){
        poly_uniform(&epsilon[j], thash, nonce++);
        poly_ntt(&epsilon[j]);
    }

    poly v[4], vtmp;
    memset(v, 0, 4*sizeof(poly));

    poly cg;
    poly_basemul_montgomery(&cg,&c,&gamma);
    poly_tomont(&cg);

    for(i=0;i<PQMX_NV;i++){
        poly_basemul_montgomery(&tmp, &beta, &f[2*PQMX_NV+i]);
        poly_tomont(&tmp);
        poly_add(&tmp, &tmp, &f[3*PQMX_NV+i]);
        poly_reduce(&tmp);
        poly_sub(&tmp,&tmp,&f[6*PQMX_NV+i]);
        poly_reduce(&tmp);
        poly_add(&tmp, &tmp, &cg);
        poly_reduce(&tmp);
        poly_basemul_montgomery(&tmp,&tmp, &epsilon[i]);
        poly_tomont(&tmp);
        poly_add(&v[0], &v[0], &tmp);
        poly_reduce(&v[0]);
    }

    for(i=0;i<PQMX_NV;i++){
        poly_basemul_montgomery(&tmp, &f[6*PQMX_NV+i], &f[7*PQMX_NV+i]);
        poly_tomont(&tmp);
        poly_basemul_montgomery(&vtmp, &f[8*PQMX_NV+i], &c);
        poly_tomont(&vtmp);
        poly_add(&vtmp, &vtmp, &tmp);
        poly_reduce(&vtmp);
        poly_basemul_montgomery(&tmp, &vtmp, &epsilon[PQMX_NV+i]);
        poly_tomont(&tmp);
        poly_add(&v[1], &v[1], &tmp);
        poly_reduce(&v[1]);

        poly_basemul_montgomery(&tmp, &f[i], &f[3*PQMX_NV+i]);
        poly_tomont(&tmp);
        poly_basemul_montgomery(&vtmp, &f[4*PQMX_NV+i], &c);
        poly_tomont(&vtmp);
        poly_add(&vtmp, &vtmp, &tmp);
        poly_reduce(&vtmp);
        poly_basemul_montgomery(&tmp, &vtmp, &epsilon[2*PQMX_NV+i]);
        poly_tomont(&tmp);
        poly_add(&v[1], &v[1], &tmp);
        poly_reduce(&v[1]);

        poly_basemul_montgomery(&tmp, &f[PQMX_NV+i], &f[3*PQMX_NV+i]);
        poly_tomont(&tmp);
        poly_basemul_montgomery(&vtmp, &f[5*PQMX_NV+i], &c);
        poly_tomont(&vtmp);
        poly_add(&vtmp, &vtmp, &tmp);
        poly_reduce(&vtmp);
        poly_basemul_montgomery(&tmp, &vtmp, &epsilon[3*PQMX_NV+i]);
        poly_tomont(&tmp);
        poly_add(&v[1], &v[1], &tmp);
        poly_reduce(&v[1]);
    }
    poly_add(&v[1], &v[1], &f[PQMX_M-1]);
    poly_reduce(&v[1]);

    poly M[2];
    memset(M, 0, 2*sizeof(poly));    
    for(i=0;i<PQMX_NV;i++){
        poly_basemul_montgomery(&tmp, &alphapowpis[i], &in[i].u);
        poly_tomont(&tmp);
        poly_add(&M[0], &M[0], &tmp);
        poly_reduce(&M[0]);

        poly_basemul_montgomery(&tmp, &alphapowpis[i], &in[i].v);
        poly_tomont(&tmp);
        poly_add(&M[1], &M[1], &tmp);
        poly_reduce(&M[1]);
    }
    poly_basemul_montgomery(&M[0], &M[0], &c);
    poly_tomont(&M[0]);
    poly_basemul_montgomery(&M[1], &M[1], &c);
    poly_tomont(&M[1]);

    memset(&vtmp, 0, sizeof(poly));
    for(i=0;i<PQMX_NV;i++){
        poly_basemul_montgomery(&tmp, &f[3*PQMX_NV+i], &out[i].u);
        poly_tomont(&tmp);
        poly_add(&vtmp, &vtmp, &tmp);
        poly_reduce(&vtmp);
        poly_sub(&vtmp, &vtmp, &f[4*PQMX_NV+i]);
        poly_reduce(&vtmp);
    }

    poly_add(&vtmp, &vtmp, &M[0]);
    poly_reduce(&vtmp);
    poly_basemul_montgomery(&tmp, &epsilon[4*PQMX_NV], &vtmp);
    poly_tomont(&tmp);
    poly_add(&v[2], &v[2], &tmp);
    poly_reduce(&v[2]);
    
    memset(&vtmp, 0, sizeof(poly));
    for(i=0;i<PQMX_NV;i++){
        poly_basemul_montgomery(&tmp, &f[3*PQMX_NV+i], &out[i].v);
        poly_tomont(&tmp);
        poly_add(&vtmp, &vtmp, &tmp);
        poly_reduce(&vtmp);
        poly_sub(&vtmp, &vtmp, &f[5*PQMX_NV+i]);
        poly_reduce(&vtmp);
    }
    poly_add(&vtmp, &vtmp, &M[1]);
    poly_reduce(&vtmp);
    poly_basemul_montgomery(&tmp, &epsilon[4*PQMX_NV+1], &vtmp);
    poly_tomont(&tmp);
    poly_add(&v[2], &v[2], &tmp);
    poly_reduce(&v[2]);



    poly_basemul_montgomery(&tmp, &P, &c);
    poly_tomont(&tmp);
    poly_add(&tmp, &tmp, &f[9*PQMX_NV-1]);
    poly_reduce(&tmp);
    poly_basemul_montgomery(&tmp, &epsilon[4*PQMX_NV+2], &tmp);
    poly_tomont(&tmp);
    poly_add(&v[3], &v[3], &tmp);

    poly_add(&tmp, &f[7*PQMX_NV], &c);
    poly_reduce(&tmp);
    poly_basemul_montgomery(&tmp, &epsilon[4*PQMX_NV+3], &tmp);
    poly_tomont(&tmp);
    poly_add(&v[3], &v[3], &tmp);
    poly_reduce(&v[3]);

    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t*)v, 4*sizeof(poly));
    shake128_absorb(&state, (uint8_t*)(t->tm+9*PQMX_NV), sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(chash, PQMX_SYMBYTES, &state);
    
    free(f);
    free(epsilon);
    free(alphapowpis);
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