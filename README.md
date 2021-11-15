# Improved Lattice-based Mix-Net implementation

This repository contains source code of a protocol from the paper *Improved Lattice-Based Mix-Nets for Electronic Voting* accepted to ICISC 2021 by Jan Willemson, Jaan Kristjan Kaasik and Valeh Farzaliyev. 

## Run
There are three test files to test Ring-LWE encryption algorithm, Zero-Knowledge proof of shortness of RLWE params, and Zero-Knowledge proof of shuffle. To compile and run test programs on Linux, run

```sh
make test_rlwe
./test_rlwe

make test_shortness
./test_shortness

make test_shuffle
./test_shuffle
```
WARNING: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.
