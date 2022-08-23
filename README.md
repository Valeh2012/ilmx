# Improved Lattice-based Mix-Net implementation

This repository contains the source code of the protocol from the paper [*Improved Lattice-Based Mix-Nets for Electronic Voting*](https://link.springer.com/chapter/10.1007/978-3-031-08896-4_6) accepted to ICISC 2021 by Jan Willemson, Jaan Kristjan Kaasik and Valeh Farzaliyev. 

UPDATE: The code has been modified partially to match the protocol specification given in the journal version of the paper, [*Improved Lattice-Based Mix-Nets for Electronic Voting*](https://ietresearch.onlinelibrary.wiley.com/doi/full/10.1049/ise2.12089), published at IET Information Security

## Run
There are four test files to test Ring-LWE encryption algorithm `test_rlwe`, Zero-Knowledge proof of shortness of RLWE params `test_shortness`, Zero-Knowledge proof of a shuffle (without shortness) `test_shuffle`, and full Zero-knowledge proof of a shuffle (with shortness) protocol `test_mx`. To compile and run test full protocol on Linux, run

```sh
make test_mx
./text_mx
```
WARNING: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.

## Citation

journal version:
```json
@article{https://doi.org/10.1049/ise2.12089,
author = {Farzaliyev, Valeh and Willemson, Jan and Kaasik, Jaan Kristjan},
title = {Improved lattice-based mix-nets for electronic voting},
journal = {IET Information Security},
volume = {n/a},
number = {n/a},
pages = {},
keywords = {electronic voting, implementation, lattice-based post-quantum cryptography, mix-nets, zero-knowledge proofs},
doi = {https://doi.org/10.1049/ise2.12089},
url = {https://ietresearch.onlinelibrary.wiley.com/doi/abs/10.1049/ise2.12089},
eprint = {https://ietresearch.onlinelibrary.wiley.com/doi/pdf/10.1049/ise2.12089},
}
```

eprint version
```json
@misc{cryptoeprint:2021:1499,
    author       = {Valeh Farzaliyev and
		    Jan Willemson and
		    Jaan Kristjan Kaasik},
    title        = {Improved Lattice-Based Mix-Nets for Electronic Voting},
    howpublished = {Cryptology ePrint Archive, Report 2021/1499},
    year         = {2021},
    note         = {\url{https://ia.cr/2021/1499}},
}
```


conference version:
```json
@InProceedings{10.1007/978-3-031-08896-4_6,
author="Farzaliyev, Valeh
and Willemson, Jan
and Kaasik, Jaan Kristjan",
editor="Park, Jong Hwan
and Seo, Seung-Hyun",
title="Improved Lattice-Based Mix-Nets for Electronic Voting",
booktitle="Information Security and Cryptology -- ICISC 2021",
year="2022",
publisher="Springer International Publishing",
address="Cham",
pages="119--136",
isbn="978-3-031-08896-4"
}
```
