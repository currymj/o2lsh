use std::fs::File;
use std::io::BufReader;
use std::io::Error;
use std::io::prelude::*;
use std::path::Path;

pub fn get_mnist_vector(fname: &str) -> Result<Vec<Vec<f64>>, Error> {
    let path = Path::new(fname);
    match File::open(&path) {
        Ok(file) => {
            let mut new_vec: Vec<Vec<f64>> = Vec::new();
            let reader = BufReader::new(file);
            for line in reader.lines() {
                match line {
                    Ok(s) => {
                        new_vec.push(mnist_test_to_vector(&s));
                    },
                    Err(reason) => {
                        return Err(reason);
                    }
                };
            }
            Ok(new_vec)
        },
        Err(reason) => Err(reason)
    }
}

pub fn mnist_test_to_vector(line: &str) -> Vec<f64> {
    line.trim().split(' ').map(|instr| instr.parse().unwrap()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_read_mnist_file() {
        let fname: &str = "./mnist1k.dts";
        match get_mnist_vector(fname) {
            Ok(v) => println!("{:?}", v),
            Err(reason) => panic!("couldn't open because {}", reason)
        }
    }

    #[test]
    fn test_mnist_line_to_vector() {
        let my_vec = vec![0.0000, 1.0000, 0.00000, 0.20000];
        let my_string = "0.0000 1.0000 0.00000 0.20000";
        assert!(my_vec == mnist_test_to_vector(&my_string));
    }
}
