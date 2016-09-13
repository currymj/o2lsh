extern crate rand;
use table::rand::Rng;
use std::ops::{Index, IndexMut};
type Bucket = Vec<usize>; // later this will have chaining, for now just blob
 // everything together
const P: f64 = (0xFFFFFFFE as usize - 4) as f64;
pub struct LSHTable<'a, T: 'a, Q: 'a+?Sized>  {
    buckets: Vec<Bucket>,
    data: &'a [T],
    hash_functions: &'a [Box<Q>],
    ri1: Vec<f64>
}

pub struct SmallPointerArray<T> {
    data: Vec<T>
}

impl<T> Index<u32> for SmallPointerArray<T> {
    type Output = T;
    fn index<'a>(&'a self, _index: u32) -> &'a T {
        &(self.data[_index as usize])
    }
}

impl<T> IndexMut<u32> for SmallPointerArray<T> {
    fn index_mut<'a>(&'a mut self, _index: u32) -> &'a mut T {
        &mut (self.data[_index as usize])
    }
}

impl<T> SmallPointerArray<T> {
    pub fn new() -> Self {
        SmallPointerArray {
            data: Vec::new()
        }
    }

    pub fn push(&mut self, x: T) {
        self.data.push(x);
    }
}

impl<'a, T, Q: 'a+?Sized> LSHTable<'a, T, Q> where Q: Fn(&'a T) -> f64 {
    pub fn new(data: &'a [T], hashes: &'a [Box<Q>]) -> Self {
        LSHTable {
            buckets: vec![Vec::new(); data.len()],
            data: data,
            hash_functions: hashes,
            ri1: rand::thread_rng().gen_iter().take(hashes.len()).collect(),
        }
    }
    fn get_signature(&self, v: &'a T) -> Vec<f64> {
        self.hash_functions.iter().map(|x| {
            (*x)(v)
        }).collect()
    }
    pub fn new_build(data: &'a [T], hashes: &'a [Box<Q>]) -> Self {
        let mut x_to_build = LSHTable::new(data, hashes);
        for (i, v) in x_to_build.data.iter().enumerate() {
            let hash_sig = x_to_build.get_signature(v);
            let bucket_ind = hash_func_t1(&hash_sig, &x_to_build.ri1, P, x_to_build.buckets.len());
            x_to_build.buckets[bucket_ind].push(i);
        }
        x_to_build
    }

    pub fn query_vec(&self, v: &'a T) -> Vec<&T> {
        let sig = self.get_signature(v);
        let sig_ind = hash_func_t1(&sig, &self.ri1, P, self.buckets.len());
        self.buckets[sig_ind].iter().map(|bucket_ind| {
            &self.data[*bucket_ind]
        }).collect()
    }
}

fn hash_func_t1(signature: &[f64], rand_ints: &[f64], primes: f64, num_buckets: usize) -> usize {
    let total = hash_func_inner(signature, rand_ints, primes);
    total % num_buckets
}
fn hash_func_inner(signature: &[f64], rand_ints: &[f64], primes: f64) -> usize {
    let total: usize = {
        let mut counter: f64 = 0.0;
        for element in signature.iter().zip(rand_ints).map(|(a, b)| {(a * b)}) {
            counter = ((counter + element)) % primes;
        }
        counter as usize
    };
    total
}

#[cfg(test)]
mod tests {
    use super::LSHTable;
    use super::SmallPointerArray;
    #[test]
    fn test_init() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];

        let val = |q: &Vec<i32>| {0.0 as f64};
        let funcs = vec![Box::new(val)];
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
        let funcs = vec![Box::new(val)];
        let x = LSHTable::new_build(&test_data, &funcs);
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
        let x = LSHTable::new_build(&test_data, &funcs);
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
}
