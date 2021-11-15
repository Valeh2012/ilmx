#ifndef PARAMS_H
#define PARAMS_H


//#define PQMX_90S	/* Uncomment this if you want the 90S variant */

#define PQMX_NAMESPACE(s) PQMX_##s

#define PQMX_N 4096
#define PQMX_L 1024
#define PQMX_Q 4611686018427365377
#define PQMX_BETA 4096


#define PQMX_SYMBYTES 512   /* size in bytes of hashes, and seeds */
#define PQMX_SSBYTES  32   /* size in bytes of shared key */

#define PQMX_POLYBYTES		4096*8


#define PQMX_MU 1
#define PQMX_LAMBDA 1
#define PQMX_NV 10

#ifndef SHORT 
#define PQMX_M (9*PQMX_NV+3)
#define PQMX_DELTA (1L<<43)
#define PQMX_POLYCOMPRESSEDBYTES    22528
#else 
#define PQMX_M  (4*PQMX_LAMBDA+20*PQMX_NV+3)  
#define PQMX_DELTA (1L<<47)
#define PQMX_POLYCOMPRESSEDBYTES    24576
#endif



#define PQMX_INDCPA_MSGBYTES       (PQMX_SYMBYTES)
#define PQMX_INDCPA_PUBLICKEYBYTES (PQMX_POLYVECBYTES + PQMX_SYMBYTES)
#define PQMX_INDCPA_SECRETKEYBYTES (PQMX_POLYVECBYTES)
#define PQMX_INDCPA_BYTES          (PQMX_POLYVECCOMPRESSEDBYTES + PQMX_POLYCOMPRESSEDBYTES)


#endif
