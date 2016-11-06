[![Build Status](https://travis-ci.org/currymj/o2lsh.svg?branch=master)](https://travis-ci.org/currymj/o2lsh)

O<sub>2</sub>LSH is a Rust crate for locality-sensitive hashing. It is heavily
inspired by E<sup>2</sup>LSH. Information about the underlying algorithms can be
found at http://www.mit.edu/~andoni/LSH/.

The code is currently licensed under MPL 2.0.

This code remains a work in progress. If you were hoping to use locality sensitive hashing right now in another Rust package, this probably won't get you there at the moment.

Goals:
[x] multiprobe
[x] rayon parallelism
[ ] validate correctness/efficiency of Euclidean LSH implementation
[ ] refactor for ergonomics
[ ] implement other hash functions
[ ] make general enough to implement some really weird hash functions
