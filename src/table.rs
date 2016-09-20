extern crate rand;
use table::rand::Rng;
use multi;


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
    multiprobe_sequence: &'a [Vec<usize>]
}

impl<'a, T, Q: 'a+?Sized> LSHTable<'a, T, Q> where Q: Fn(&'a T) -> f64 {
    pub fn new(data: &'a [T], hashes: &'a [Box<Q>], ms: &'a [Vec<usize>]) -> Self {
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
    fn get_quantized_signature(&self, v: &'a T) -> Vec<f64> {
        self.hash_functions.iter().map(|x| {
            (*x)(v).floor()
        }).collect()
    }
    pub fn new_build(data: &'a [T], hashes: &'a [Box<Q>], ms: &'a [Vec<usize>]) -> Self {
        let mut x_to_build = LSHTable::new(data, hashes, ms);
        for (i, v) in x_to_build.data.iter().enumerate() {
            let hash_sig = x_to_build.get_quantized_signature(v);
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
        let sig = self.get_quantized_signature(v);
        let sig_ind = hash_func_t1(&sig, &self.ri1, P, self.buckets.len());
        let chain_ind = hash_func_t2(&sig, &self.ri2, P );
        match self.buckets[sig_ind].get(chain_ind) {
           Some(bucket) => bucket.pointers.iter().map(|bucket_ind| {
                &self.data[*bucket_ind]
           }).collect(),
            None => Vec::new()
        }
    }

    fn get_all_sigs(&self, sig: &[f64]) -> Vec<Vec<f64>> {
        let qsig = self.quantize_signature(sig);
        let pi_map = multi::compute_pi_j(sig, &qsig, 1.0);
        let mut combined_output = Vec::new();
        for perturb_set in self.multiprobe_sequence {
            let mut out_vec = qsig.clone();
            for &val in perturb_set {
                let (i, delta) = pi_map[val];
                out_vec[i] += delta as f64;
            }
            combined_output.push(out_vec);
        }
        combined_output
    }
    
    fn quantize_signature(&self, sig: &[f64]) -> Vec<f64> {
        sig.iter().map(|x| {x.floor()}).collect()
    }

    pub fn query_multiprobe(&self, v: &'a T) -> Vec<&T> {
        let sig = self.get_signature(v);
        let all_sigs = self.get_all_sigs(&sig);
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
    use multi;


    // fn to gen multiprobe sequence
    #[test]
    fn test_init() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];


        let val = |q: &Vec<i32>| {0.0 as f64};
        let funcs = vec![Box::new(val)];
        let ms = vec![vec![1,2,3]];
        // get actual multiprobe sequence into ms instead
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
        let ms = vec![vec![1,2,3]];
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
        let zjs = multi::get_expected_zj_vals(1,1.0);
        let sets: Vec<multi::PerturbationSet> = multi::gen_perturbation_sets(&zjs)
            .take(5)
            .collect();
        let ms: Vec<Vec<usize>> = sets.into_iter()
            .map(|x| {x.data})
            .collect();
        let x = LSHTable::new_build(&test_data, &funcs, &ms);
        x.query_vec(&test_data[0]);
    }
}
