extern crate rand;
// p-stable hash function
use self::rand::distributions::{Normal, IndependentSample};
use self::rand::Rng;
fn get_random_normal_vector(length: usize) -> Vec<f64> {
    let normal = Normal::new(0.0,1.0);
    let mut output_vector: Vec<f64> = vec![0.0; length];
    for i in 1..length {
        output_vector[i] = normal.ind_sample(&mut rand::thread_rng());
    }
    output_vector
}

pub fn get_hash_closure<'a>(length: usize, w: f64) -> Box<Fn(&Vec<f64>) -> f64> {
    let a = get_random_normal_vector(length);
    let b = rand::thread_rng().gen_range(0.0, w);
    let clos = move |v: &Vec<f64>| -> f64 {
        let total: f64 = v.iter()
            .zip(a.iter())
            .map(|(x, y)| x * y)
            .sum();
        (total + b) / w
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
