extern crate rand;
use table::rand::Rng;
type Bucket = Vec<usize>; // later this will have chaining, for now just blob
 // everything together
const P: i32 = 0xFFFFFFFF - 5;
struct LSHTable<'a, T: 'a, Q: 'a> {
    buckets: Vec<Bucket>,
    data: &'a Vec<T>,
    hash_functions: &'a Vec<&'a Q>,
    ri1: Vec<f64>,
    ri2: Vec<f64>
}

impl<'a, T, Q: 'a> LSHTable<'a, T, Q> where Q: Fn(&'a T) -> f64 {
    pub fn new(data: &'a Vec<T>, hashes: &'a Vec<&'a Q>) -> Self {
        LSHTable {
            buckets: vec![Vec::new(); data.len() ],
            data: data,
            hash_functions: hashes,
            ri1: rand::thread_rng().gen_iter().take(hashes.len()).collect(),
            ri2: rand::thread_rng().gen_iter().take(hashes.len()).collect()
        }
    }
    pub fn new_build(data: &'a Vec<T>, hashes: &'a Vec<&'a Q>) -> Self {
        let mut x_to_build = LSHTable::new(data, hashes);
        for (i, v) in x_to_build.data.iter().enumerate() {
            let hash_sig: Vec<f64> = hashes.iter().map(|x| {
                (*x)(v)
            }).collect();
            let bucket_ind = hash_func_t1(&hash_sig, &x_to_build.ri1, P, x_to_build.buckets.len());
            x_to_build.buckets[bucket_ind].push(i);
        }
        x_to_build
    }
}

fn hash_func_t1(signature: &Vec<f64>, rand_ints: &Vec<f64>, primes: i32, num_buckets: usize) -> usize {
    let total: usize = signature.iter().zip(rand_ints).map(|(a, b)| {(a * b) as usize}).sum();
    (total % (P as usize)) % num_buckets
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
        let funcs = vec![&val];
        let x = LSHTable::new(&test_data, &funcs);
    }
    #[test]
    fn test_new_build() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];

        let val = |q: &Vec<i32>| {0.0 as f64};
        let funcs = vec![&val];
        let x = LSHTable::new_build(&test_data, &funcs);
    }
}
