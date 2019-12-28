## QKD

This repo contains code used for implementing post-processing in a QKD system, and borrows ideas from a bunch of research papers. Introductory information is provided below.

## Schematics

1.ldpc.py - Implementation of the hard decoding and sum-product soft decoding LDPC algorithms for error correction. Functions for Gallager construction and Mackay-Neal construction are provided.

2.compressedSensing.py - Perform error correction using signal reconstruction by compressed sensing methods. We minimise the L1 norm for errors using linear programming.

3.localSearch.py - Reduce the error correction problem to a multidimensional subset sum problem, and try to solve using it random-greedy moves. Lots of heuristic tuning required here.

4.base.py - Helper file for other methods

5.hashFunctions.py - Minimal implementations of Toeplitz and Polynomial Hashing for privacy amplification

6.primes.txt - A set of primes of the form 2^a + b , where (a,b) are stored in the file

## Resources

[Lightweight authentication for quantum key distribution](https://arxiv.org/pdf/1903.10237.pdf)

[Post-processing procedure for industrial quantum
key distribution systems](https://iopscience.iop.org/article/10.1088/1742-6596/741/1/012081/pdf)

[Introducing Low-Density Parity-Check Codes
](https://www.researchgate.net/publication/228977165_Introducing_Low-Density_Parity-Check_Codes)

[Quantum Cryptography Notes](http://users.cms.caltech.edu/~vidick/notes/QCryptoX/LN_Week0.pdf)
