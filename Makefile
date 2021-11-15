CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer
NISTFLAGS += -Wno-unused-result -O3 -fomit-frame-pointer
RM = /bin/rm

SOURCES = polyvec.c poly.c ntt.c reduce.c util.c
SOURCESKECCAK = $(SOURCES) fips202.c symmetric-shake.c
HEADERS = params.h polyvec.h poly.h ntt.h reduce.h util.h 
HEADERSKECCAK = $(HEADERS) fips202.h

.PHONY: all speed shared clean

all: \
  test_rlwe \
  test_shuffle


%_short.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -DSHORT=1 -c $< -o $@

%_shuffle.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

test_shuffle: $(SOURCESKECCAK) rlwe.c comm_shuffle.o comm_short.o shortness_proof_short.o mxproof.c test_shuffle.c randombytes.c 
	$(CC) $(CFLAGS) $^ -o test_shuffle -lm

test_shortness: $(SOURCESKECCAK) randombytes.c shortness_proof_short.o comm_short.o comm_shuffle.o rlwe.c test_shortness.c
	$(CC) $(CFLAGS) $^ -o test_shortness -lm

test_rlwe: $(SOURCESKECCAK) randombytes.c rlwe.c test_rlwe.c 
	$(CC) $(CFLAGS) $^ -o test_rlwe -lm

clean:
	-$(RM) -rf *.gcno *.gcda *.gch *.lcov *.o *.so
	-$(RM) -rf test_rlwe
	-$(RM) -rf test_shuffle
	-$(RM) -rf test_shortness
