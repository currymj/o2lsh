extern crate rand;
use table::rand::Rng;
use std::ops::{Index, IndexMut};
type Bucket = Vec<usize>; // later this will have chaining, for now just blob
// everything together

const P: f64 = (0xFFFFFFFE as usize - 4) as f64;
const MASK: u32 = 0x000FFFFF;

pub struct LSHTable<'a, T: 'a, Q: 'a+?Sized>  {
    buckets: Vec<Bucket>,
    data: &'a [T],
    hash_functions: &'a [Box<Q>],
    ri1: Vec<f64>
}

pub struct MaskedSmallPointerArray<T> {
    data: Vec<T>
}
impl<T> Index<u32> for MaskedSmallPointerArray<T> {
    type Output = T;
    fn index(&self, _index: u32) -> &T {
        &(self.data[(_index & MASK)  as usize])
    }
}

impl<T> IndexMut<u32> for MaskedSmallPointerArray<T> {
    fn index_mut(&mut self, _index: u32) -> &mut T {
        &mut (self.data[(_index & MASK) as usize])
    }
}
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
    pub fn new() -> Self {
        SmallPointerArray {
            data: Vec::new()
        }
    }

    pub fn push(&mut self, x: T) {
        self.data.push(x);
    }
}

impl<T> MaskedSmallPointerArray<T> {
    pub fn new() -> Self {
        MaskedSmallPointerArray {
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
            let bucket_ind = hash_func_t1(&hash_sig, &x_to_build.ri1, P, x_to_build.buckets.len() as u32) as usize;
            x_to_build.buckets[bucket_ind].push(i);
        }
        x_to_build
    }

    pub fn query_vec(&self, v: &'a T) -> Vec<&T> {
        let sig = self.get_signature(v);
        let sig_ind = hash_func_t1(&sig, &self.ri1, P, (self.buckets.len() as u32)) as usize;
        self.buckets[sig_ind].iter().map(|bucket_ind| {
            &self.data[*bucket_ind]
        }).collect()
    }
}

fn hash_func_t1(signature: &[f64], rand_ints: &[f64], primes: f64, num_buckets: u32) -> u32 {
    let total = hash_func_inner(signature, rand_ints, primes);
    total % num_buckets
}
fn hash_func_inner(signature: &[f64], rand_ints: &[f64], primes: f64) -> u32 {
    let mut counter: f64 = 0.0;
    for element in signature.iter().zip(rand_ints).map(|(a, b)| {(a * b)}) {
        counter = ((counter + element)) % primes;
    }
    counter as u32
}

#[cfg(test)]
mod tests {
    use super::LSHTable;
    use super::SmallPointerArray;
    use super::MaskedSmallPointerArray;
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
    #[test]
    fn test_masked_array() {
        let mut test_arr: MaskedSmallPointerArray<i32> = MaskedSmallPointerArray::new();
        for i in 1..10 {
            test_arr.push(i);
        }
        println!("{}", test_arr.data.len());
        println!("test!");
        for i in 0xFFF00000u32..0xFFF00009 {
            println!("{}", test_arr[i]);
        };
    }
}
