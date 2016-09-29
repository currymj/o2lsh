extern crate rand;
extern crate ndarray;

use self::rand::distributions::{Normal, IndependentSample};
use self::rand::Rng;
use self::ndarray::prelude::*;

/// Returns a random normal vector of size `length`
fn get_random_normal_vector(length: usize) -> Vec<f32> {
    let normal = Normal::new(0.0,1.0);
    let output_vector: Vec<f32> = vec![normal.ind_sample(&mut rand::thread_rng()) as f32; length];
    output_vector
}

/// Returns a boxed closure, which will dot its input of length `length` with a random normal vector, add a random constant, and scale by `w`.
pub fn get_hash_closure(length: usize, w: f32) -> Box<Fn(&Vec<f32>) -> f32> {
    let a = get_random_normal_vector(length);
    let b = rand::thread_rng().gen_range(0.0, w);
    let clos = move |v: &Vec<f32>| -> f32 {
        (aview1(v).dot(&aview1(&a)) + b) / w
    };
    Box::new(clos)
}

#[cfg(test)]
mod tests {
    use super::get_random_normal_vector;
    #[test]
    fn test_get_normal() {
        get_random_normal_vector(10);
    }
}
