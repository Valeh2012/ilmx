#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define MONT 90108 // 2^64 mod q
#define MONT2 8119451664 // 2^128 mod q
#define QINV 5842042167526184961 // q^-1 mod 2^64

#define montgomery_reduce PQMX_NAMESPACE(montgomery_reduce)
int64_t montgomery_reduce(__int128 a);

#define barrett_reduce PQMX_NAMESPACE(barrett_reduce)
int64_t barrett_reduce(int64_t a);

#endif
