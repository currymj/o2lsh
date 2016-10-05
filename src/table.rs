extern crate rand;
extern crate num;
use table::rand::Rng;
use std::ops::{Index, IndexMut};
use self::num::pow;
const MASK: u32 = 0x000FFFFF;

use multi;


#[derive(Clone)]
pub struct Bucket {
    pointers: Vec<usize>,
    hash_sig: u32
}

impl Bucket {
    pub fn new(hash_sig: u32) -> Self {
        Bucket {
            pointers: Vec::new(),
            hash_sig: hash_sig
        }

    }
    pub fn push(&mut self, value: usize) {
        self.pointers.push(value);
    }
}

#[derive(Clone, Default)]
struct BucketChain {
    chain: Vec<Bucket>
}

impl BucketChain {
    pub fn new() -> Self {
        BucketChain{
            chain: Vec::new()
        }
    }
    pub fn get(&self, _index: u32) -> Option<&Bucket> {
        for bucket_ref in &self.chain {
            if bucket_ref.hash_sig == _index {
                return Some(bucket_ref);
            }
        }
        None
    }
    pub fn get_mut(&mut self, _index: u32) -> Option<&mut Bucket> {
        for bucket_ref in &mut self.chain {
            if bucket_ref.hash_sig == _index {
                return Some(bucket_ref);
            }
        }
        None
    }
    pub fn create_and_add_bucket(&mut self, bucket_sig: u32, value: usize) {
        let mut new_bucket = Bucket::new(bucket_sig);
        new_bucket.push(value);
        self.chain.push(new_bucket);
    }

    pub fn len(&self) -> usize {
        self.chain.len()
    }
}

const P: u64 = (0xFFFFFFFE as u64 - 4) as u64;

pub struct MaskedSmallPointerSlice<'a, T: 'a> {
    data: &'a [T]
}

impl<'a, T> Index<u32> for MaskedSmallPointerSlice<'a, T> {
    type Output = T;
    fn index(&self, _index: u32) -> &T {
        &(self.data[(_index & MASK)  as usize])
    }
}

/*impl<'a, T> IndexMut<u32> for MaskedSmallPointerSlice<'a, T> {
    fn index_mut(&mut self, _index: u32) -> &mut T {
        &mut (self.data[(_index & MASK) as usize])
    }
}*/

pub struct SmallPointerArray<T> {
    data: Vec<T>
}


impl<T> Index<u32> for SmallPointerArray<T> {
    type Output = T;
    fn index(&self, _index: u32) -> &T {
        &(self.data[_index as usize])
    }
}

impl<T> IndexMut<u32> for SmallPointerArray<T> {
    fn index_mut(&mut self, _index: u32) -> &mut T {
        &mut (self.data[_index as usize])
    }
}

impl<T> SmallPointerArray<T> {
    pub fn with_capacity(length: usize) -> Self {
        SmallPointerArray {
            data: Vec::with_capacity(length)
        }
    }
    pub fn new() -> Self {
        SmallPointerArray {
            data: Vec::new()
        }
    }

    pub fn push(&mut self, x: T) {
        self.data.push(x);
    }
}

impl<'a, T: 'a> MaskedSmallPointerSlice<'a, T> {
    pub fn new(d: &'a [T]) -> Self {
        MaskedSmallPointerSlice {
            data: d
        }
    }

    /*pub fn push(&mut self, x: T) {
        self.data.push(x);
    }*/
}

pub struct LSHTable<'a, T: 'a, Q: 'a+?Sized>  {
    buckets: Vec<BucketChain>,
    data: &'a [T],
    hash_functions: Vec<Box<Q>>,
    ri1: Vec<u32>,
    ri2: Vec<u32>,
    multiprobe_sequence: &'a [Vec<usize>]
}

struct PackedLSHTable<'a, T: 'a, Q: 'a+?Sized> {
    data_refs_and_buckets: SmallPointerArray<u32>,
    bucket_table: Vec<u32>,
    data: MaskedSmallPointerSlice<'a, T>,
    hash_functions: Vec<Box<Q>>,
    ri1: Vec<u32>,
    ri2: Vec<u32>,
    multiprobe_sequence: &'a [Vec<usize>]
}

fn fill_array(buckets: &mut SmallPointerArray<u32>, chain: &BucketChain, start_ind: u32) -> u32 {
    let mut ind = start_ind;
    let num_buckets = chain.chain.len();
    for i in 0..num_buckets-1 {
        buckets[ind] = chain.chain[i].hash_sig;
        ind += 1;
        let current_bucket_len = chain.chain[i].pointers.len();
        if current_bucket_len > (pow(2,11) - 1) {
            // handle long case
            unimplemented!();
        } else {
            // handle normal case
            // make sure high order 12 bits are 0
            // leave highest bit as 0
            // set remaining 11 to bucket_len
            let curr_pointer = chain.chain[i].pointers[0];
            buckets[ind] = ((current_bucket_len as u32) << 20) | (curr_pointer as u32);
            for j in 1..current_bucket_len {
                let curr_pointer = chain.chain[i].pointers[j] as u32;
                assert!(curr_pointer < pow(2,20));
                buckets[ind] = curr_pointer;
                ind += 1;
            }
        }
    }
    // deal with last bucket in chain
    // set highest bit to 1


    let curr_pointer = chain.chain[num_buckets-1].pointers[0];
    let current_bucket_len = chain.chain[num_buckets - 1].pointers.len();
    buckets[ind] = 0x80000000 | (curr_pointer as u32);
    ind += 1;

    for j in 1..current_bucket_len {
        let curr_pointer = chain.chain[num_buckets-1].pointers[j];
        assert!(curr_pointer < pow(2,20));
        buckets[ind] = curr_pointer as u32;
        ind += 1;
    }
    ind
}
impl<'a, T, Q: 'a+?Sized> PackedLSHTable<'a, T, Q> where Q: Fn(&'a T) -> f32 {


    pub fn new_from_table(inner_table: LSHTable<'a, T, Q>) -> Self {
        // unimplemented!();
        // get total number of points, and total number of buckets
        // (slow, requiring iteration over whole list)
        let num_points = inner_table.data.len();
        let num_buckets = inner_table.buckets.iter().map(|x| {x.len()}).sum();

        let mut bucket_array = SmallPointerArray::with_capacity(num_points + num_buckets);

        let mut hash_table = Vec::with_capacity(num_buckets);
        let mut current_write_ind = 0;

        let data_slice = MaskedSmallPointerSlice::new(inner_table.data);

        for (i, chain) in inner_table.buckets.iter().enumerate() {
            hash_table[i] = current_write_ind;
            current_write_ind = fill_array(&mut bucket_array, chain, current_write_ind);
        }



        PackedLSHTable {
            data_refs_and_buckets: bucket_array,
            bucket_table: hash_table,
            data: data_slice,
            hash_functions: inner_table.hash_functions,
            ri1: inner_table.ri1,
            ri2: inner_table.ri2,
            multiprobe_sequence: inner_table.multiprobe_sequence
        }

    }

    pub fn new_build(data: &'a [T], hashes: Vec<Box<Q>>, ms: &'a [Vec<usize>]) -> Self {
        PackedLSHTable::new_from_table(LSHTable::new_build(data, hashes, ms))
    }

    fn get_signature(&self, v: &'a T) -> Vec<f32> {
        self.hash_functions.iter().map(|x| {
            (*x)(v)
        }).collect()
    }

    fn get_quantized_signature(&self, v: &'a T) -> Vec<u32> {
        self.hash_functions.iter().map(|x| {
            (*x)(v) as u32
        }).collect()
    }
    /*pub fn query_vec(&self, v: &'a T) -> Vec<usize> {
        let sig = self.get_quantized_signature(v);
        let sig_ind = hash_func_t1(&sig, &self.ri1, self.buckets.len());
        let chain_ind = hash_func_t2(&sig, &self.ri2);
        // replace this part
        match self.buckets[sig_ind].get(chain_ind) {
           Some(bucket) => bucket.pointers.iter().map(|bucket_ind| {
                *bucket_ind
           }).collect(),
            None => Vec::new()
        }
    }*/

    fn get_all_sigs(&self, sig: &[f32], multiprobe_limit: usize) -> Vec<Vec<u32>> {
        let qsig = self.quantize_signature(sig);
        let pi_map = multi::compute_pi_j(sig, &qsig, 1.0);
        let mut combined_output = Vec::new();
        for perturb_set in self.multiprobe_sequence.iter().take(multiprobe_limit) {
            let mut out_vec = qsig.clone();
            for &val in perturb_set {
                let (i, delta) = pi_map[val];
                if delta == -1 {
                    out_vec[i] -= 1;
                } else {
                    out_vec[i] += delta as u32;
                }
            }
            combined_output.push(out_vec);
        }
        combined_output
    }

    fn quantize_signature(&self, sig: &[f32]) -> Vec<u32> {
        sig.iter().map(|&x| {x as u32}).collect()
    }

    /*pub fn query_multiprobe(&self, v: &'a T, multiprobe_limit: usize) -> Vec<usize> {
        let sig = self.get_signature(v);
        let all_sigs = self.get_all_sigs(&sig, multiprobe_limit);
        let mut output_vec = Vec::new();
        for s in &all_sigs {
            // replace this part
            let sig_ind = hash_func_t1(s, &self.ri1, self.buckets.len());
            let chain_ind = hash_func_t2(s, &self.ri2);
            output_vec.append(
                &mut match self.buckets[sig_ind].get(chain_ind) {
                    Some(bucket) => bucket.pointers.iter().map(|bucket_ind| {
                        *bucket_ind
                    }).collect(),
                    None => Vec::new()
                });
        }
        output_vec
    }*/
}

impl<'a, T, Q: 'a+?Sized> LSHTable<'a, T, Q> where Q: Fn(&'a T) -> f32 {

    pub fn new(data: &'a [T], hashes: Vec<Box<Q>>, ms: &'a [Vec<usize>]) -> Self {
        let length_hash = hashes.len();
        LSHTable {
            buckets: vec![BucketChain::new(); data.len()],
            data: data,
            hash_functions: hashes,
            ri1: rand::thread_rng().gen_iter().take(length_hash).collect(),
            ri2: rand::thread_rng().gen_iter().take(length_hash).collect(),
            multiprobe_sequence: ms
        }
    }

    fn get_signature(&self, v: &'a T) -> Vec<f32> {
        self.hash_functions.iter().map(|x| {
            (*x)(v)
        }).collect()
    }

    fn get_quantized_signature(&self, v: &'a T) -> Vec<u32> {
        self.hash_functions.iter().map(|x| {
            (*x)(v) as u32
        }).collect()
    }

    pub fn new_build(data: &'a [T], hashes: Vec<Box<Q>>, ms: &'a [Vec<usize>]) -> Self {
        let mut x_to_build = LSHTable::new(data, hashes, ms);
        for (i, v) in x_to_build.data.iter().enumerate() {
            let hash_sig = x_to_build.get_quantized_signature(v);
            let bucket_ind = hash_func_t1(&hash_sig, &x_to_build.ri1, x_to_build.buckets.len());
            let mut bucket_chain = &mut (x_to_build.buckets[bucket_ind]);
            let chain_ind = hash_func_t2(&hash_sig, &x_to_build.ri2);
            let mut bad = false;
            {
                let bucket_opt = bucket_chain.get_mut(chain_ind);
                match bucket_opt {
                    Some(bucket_ref) => bucket_ref.push(i),
                    None => bad = true
                };
            }
            if bad {
                bucket_chain.create_and_add_bucket(chain_ind, i);
            }
        }
        x_to_build
    }

    pub fn query_vec(&self, v: &'a T) -> Vec<usize> {
        let sig = self.get_quantized_signature(v);
        let sig_ind = hash_func_t1(&sig, &self.ri1, self.buckets.len());
        let chain_ind = hash_func_t2(&sig, &self.ri2);
        match self.buckets[sig_ind].get(chain_ind) {
           Some(bucket) => bucket.pointers.iter().map(|bucket_ind| {
                *bucket_ind
           }).collect(),
            None => Vec::new()
        }
    }

    fn get_all_sigs(&self, sig: &[f32], multiprobe_limit: usize) -> Vec<Vec<u32>> {
        let qsig = self.quantize_signature(sig);
        let pi_map = multi::compute_pi_j(sig, &qsig, 1.0);
        let mut combined_output = Vec::new();
        for perturb_set in self.multiprobe_sequence.iter().take(multiprobe_limit) {
            let mut out_vec = qsig.clone();
            for &val in perturb_set {
                let (i, delta) = pi_map[val];
                if delta == -1 {
                    out_vec[i] -= 1;
                } else {
                    out_vec[i] += delta as u32;
                }
            }
            combined_output.push(out_vec);
        }
        combined_output
    }

    fn quantize_signature(&self, sig: &[f32]) -> Vec<u32> {
        sig.iter().map(|&x| {x as u32}).collect()
    }

    pub fn query_multiprobe(&self, v: &'a T, multiprobe_limit: usize) -> Vec<usize> {
        let sig = self.get_signature(v);
        let all_sigs = self.get_all_sigs(&sig, multiprobe_limit);
        let mut output_vec = Vec::new();
        for s in &all_sigs {
            let sig_ind = hash_func_t1(s, &self.ri1, self.buckets.len());
            let chain_ind = hash_func_t2(s, &self.ri2);
            output_vec.append(
                &mut match self.buckets[sig_ind].get(chain_ind) {
                    Some(bucket) => bucket.pointers.iter().map(|bucket_ind| {
                        *bucket_ind
                    }).collect(),
                    None => Vec::new()
                });
        }
        output_vec
    }
}


fn hash_func_t1(signature: &[u32], rand_ints: &[u32], num_buckets: usize) -> usize{
    let total = hash_func_t2(signature, rand_ints);
    (total as usize) % num_buckets
}

fn hash_func_inner(signature: &[f64], rand_ints: &[f64], primes: f64) -> u32 {
    let mut counter: f64 = 0.0;
    for element in signature.iter().zip(rand_ints).map(|(a, b)| {(a * b)}) {
        counter = ((counter + element)) % primes;
    }
    counter as u32
}

const HIGH_ORDER: u64 = 0xFFFF0000;

fn hash_func_t2(signature: &[u32], rand_ints: &[u32]) -> u32 {
    let mut total = 0 as u32;
    for (&i, &j) in signature.iter().zip(rand_ints) {
        let i_product = (i * j) as u64;
        let alpha = (i_product & HIGH_ORDER) + 5 * (i_product >> 32);
        if alpha < P {
            total = total + (alpha as u32);
        } else {
            total = (total + ((alpha - P) as u32)) as u32;
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::LSHTable;
    use super::SmallPointerArray;
    use super::MaskedSmallPointerSlice;
    use multi;

    // fn to gen multiprobe sequence
    #[test]
    fn test_init() {

        let test_data = vec![
            vec![1.0,2.0,3.0,4.0,5.0],
            vec![0.0,0.0,0.0,0.0,0.0],
            vec![1.0,2.0,3.0,4.0,5.0]
        ];

        let val = |q: &Vec<f32>| {0.0 as f32};
        let funcs = vec![Box::new(val)];
        let ms = vec![vec![1,2,3]];
        // get actual multiprobe sequence into ms instead
        let x = LSHTable::new(&test_data, funcs, &ms);
    }
    #[test]
    fn test_new_build() {
        let test_data = vec![
            vec![1.0,2.0,3.0,4.0,5.0],
            vec![0.0,0.0,0.0,0.0,0.0],
            vec![1.0,2.0,3.0,4.0,5.0]
        ];

        let val = |q: &Vec<f32>| {0.0 as f32};
        let funcs = vec![Box::new(val)];
        let ms = vec![vec![1,2,3]];
        let x = LSHTable::new_build(&test_data, funcs, &ms);
    }

    #[test]
    fn test_query() {
        let test_data = vec![
            vec![1.0,2.0,3.0,4.0,5.0],
            vec![0.0,0.0,0.0,0.0,0.0],
            vec![1.0,2.0,3.0,4.0,5.0]
        ];

        let val = |q: &Vec<f32>| {0.0 as f32};
        let funcs = vec![Box::new(val)];
        let zjs = multi::get_expected_zj_vals(1,1.0);
        let sets: Vec<multi::PerturbationSet> = multi::gen_perturbation_sets(&zjs)
            .take(5)
            .collect();
        let ms: Vec<Vec<usize>> = sets.into_iter()
            .map(|x| {x.data})
            .collect();
        let x = LSHTable::new_build(&test_data, funcs, &ms);
        x.query_vec(&test_data[0]);
    }

    #[test]
    fn test_small_array() {
        let mut test_arr: SmallPointerArray<i32> = SmallPointerArray::new();
        for i in 1..10 {
            test_arr.push(i);
        }
        println!("{}", test_arr.data.len());
        println!("test!");
        for i in 0u32..9 {
            println!("{}", test_arr[i]);
        };
    }
    #[test]
    fn test_masked_array() {
        let mut test_vec = Vec::new();
            for i in 1..10 {
                test_vec.push(i);
            }
        let test_arr: MaskedSmallPointerSlice<i32> = MaskedSmallPointerSlice::new(&test_vec);
        println!("{}", test_arr.data.len());
        println!("test!");
        for i in 0xFFF00000u32..0xFFF00009 {
            println!("{}", test_arr[i]);
        };
    }
    #[test]
    fn test_query_multiprobe() {
        let test_data = vec![
            vec![1.0,2.0,3.0,4.0,5.0],
            vec![0.0,0.0,0.0,0.0,0.0],
            vec![1.0,2.0,3.0,4.0,5.0]
        ];

        let val = |q: &Vec<f32>| {0.0 as f32};
        let funcs = vec![Box::new(val)];
        let zjs = multi::get_expected_zj_vals(1,1.0);
        let sets: Vec<multi::PerturbationSet> = multi::gen_perturbation_sets(&zjs)
            .take(5)
            .collect();
        let ms: Vec<Vec<usize>> = sets.into_iter()
            .map(|x| {x.data})
            .collect();
        let x = LSHTable::new_build(&test_data, funcs, &ms);
        x.query_multiprobe(&test_data[0], 3);
    }
}
