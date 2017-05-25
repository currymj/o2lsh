//! This is O2LSH, a crate for locality-sensitive hashing. The idea is to efficiently (linear time) build a hash table 
//! such that similar items collide with high probability, and dissimilar items collide with low probability. One can
//! adjust the precision of the hash function, and increase recall with redundant tables.

pub mod table;
pub mod hashes;
pub mod util;
pub mod multi;
pub mod lsh;
pub mod lshtable;

#[cfg(test)]
mod tests {}
