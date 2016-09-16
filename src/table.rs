extern crate rand;
use table::rand::Rng;

#[derive(Clone)]
pub struct Bucket {
    pointers: Vec<usize>,
    hash_sig: usize
}

impl Bucket {
    pub fn new(hash_sig: usize) -> Self {
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
pub struct BucketChain {
    chain: Vec<Bucket>
}

impl BucketChain {
    pub fn new() -> Self {
        BucketChain{
            chain: Vec::new()
        }
    }
    pub fn get(&self, _index: usize) -> Option<&Bucket> {
        for bucket_ref in &self.chain {
            if bucket_ref.hash_sig == _index {
                return Some(bucket_ref);
            }
        }
        None
    }
    pub fn get_mut(&mut self, _index: usize) -> Option<&mut Bucket> {
        for bucket_ref in &mut self.chain {
            if bucket_ref.hash_sig == _index {
                return Some(bucket_ref);
            }
        }
        None
    }
    pub fn create_and_add_bucket(&mut self, bucket_sig: usize, value: usize) {
        let mut new_bucket = Bucket::new(bucket_sig);
        new_bucket.push(value);
        self.chain.push(new_bucket);
    }
}

const P: f64 = (0xFFFFFFFE as usize - 4) as f64;
pub struct LSHTable<'a, T: 'a, Q: 'a+?Sized>  {
    buckets: Vec<BucketChain>,
    data: &'a [T],
    hash_functions: &'a [Box<Q>],
    ri1: Vec<f64>,
    ri2: Vec<f64>,
    multiprobe_sequence: &'a [Vec<i32>]
}

impl<'a, T, Q: 'a+?Sized> LSHTable<'a, T, Q> where Q: Fn(&'a T) -> f64 {
    pub fn new(data: &'a [T], hashes: &'a [Box<Q>], ms: &'a [Vec<i32>]) -> Self {
        LSHTable {
            buckets: vec![BucketChain::new(); data.len()],
            data: data,
            hash_functions: hashes,
            ri1: rand::thread_rng().gen_iter().take(hashes.len()).collect(),
            ri2: rand::thread_rng().gen_iter().take(hashes.len()).collect(),
            multiprobe_sequence: ms
        }
    }
    fn get_signature(&self, v: &'a T) -> Vec<f64> {
        self.hash_functions.iter().map(|x| {
            (*x)(v)
        }).collect()
    }
    pub fn new_build(data: &'a [T], hashes: &'a [Box<Q>], ms: &'a [Vec<i32>]) -> Self {
        let mut x_to_build = LSHTable::new(data, hashes, ms);
        for (i, v) in x_to_build.data.iter().enumerate() {
            let hash_sig = x_to_build.get_signature(v);
            let bucket_ind = hash_func_t1(&hash_sig, &x_to_build.ri1, P, x_to_build.buckets.len());
            let mut bucket_chain = &mut (x_to_build.buckets[bucket_ind]);
            let chain_ind = hash_func_t2(&hash_sig, &x_to_build.ri2, P );
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

    pub fn query_vec(&self, v: &'a T) -> Vec<&T> {
        let sig = self.get_signature(v);
        let sig_ind = hash_func_t1(&sig, &self.ri1, P, self.buckets.len());
        let chain_ind = hash_func_t2(&sig, &self.ri2, P );
        match self.buckets[sig_ind].get(chain_ind) {
           Some(bucket) => bucket.pointers.iter().map(|bucket_ind| {
                &self.data[*bucket_ind]
           }).collect(),
            None => Vec::new()
        }
    }

    fn perturb_signature(&self, sig: Vec<f64>, v: &T) -> Vec<Vec<f64>> {
        // precompute the multiprobe data...where?
        // don't worry about it, worry about this end for now
        // for j in M_list
        // get pi_j(that) from sig and v
        // add to output vector
        unimplemented!();
    }

    pub fn query_multiprobe(&self, v: &'a T) -> Vec<&T> {
        let sig = self.get_signature(v);
        let all_sigs = self.perturb_signature(sig, v);
        let mut output_vec = Vec::new();
        for s in &all_sigs {
            let sig_ind = hash_func_t1(s, &self.ri1, P, self.buckets.len());
            let chain_ind = hash_func_t2(s, &self.ri2, P);
            output_vec.append(
                &mut match self.buckets[sig_ind].get(chain_ind) {
                    Some(bucket) => bucket.pointers.iter().map(|bucket_ind| {
                        &self.data[*bucket_ind]
                    }).collect(),
                    None => Vec::new()
                });
        }
        output_vec
    }
}

fn hash_func_t1(signature: &[f64], rand_ints: &[f64], primes: f64, num_buckets: usize) -> usize {
    let total: usize = hash_func_t2(signature, rand_ints, primes);
    total % num_buckets
}

fn hash_func_t2(signature: &[f64], rand_ints: &[f64], primes: f64) -> usize {
    let mut counter: f64 = 0.0;
    for element in signature.iter().zip(rand_ints).map(|(a, b)| {(a * b)}) {
        counter = ((counter + element)) % primes;
    }
    counter as usize
}
#[cfg(test)]
mod tests {
    use super::LSHTable;
    #[test]
    fn test_init() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];

        let val = |q: &Vec<i32>| {0.0 as f64};
        let funcs = vec![Box::new(val)];
        let ms = vec![vec![1i32,2,3]];
        let x = LSHTable::new(&test_data, &funcs, &ms);
    }
    #[test]
    fn test_new_build() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];

        let val = |q: &Vec<i32>| {0.0 as f64};
        let funcs = vec![Box::new(val)];
        let ms = vec![vec![1i32,2,3]];
        let x = LSHTable::new_build(&test_data, &funcs, &ms);
    }

    #[test]
    fn test_query() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];

        let val = |q: &Vec<i32>| {0.0 as f64};
        let funcs = vec![Box::new(val)];
        let ms = vec![vec![1i32,2,3]];
        let x = LSHTable::new_build(&test_data, &funcs, &ms);
        x.query_vec(&test_data[0]);
    }
}
