use std::fs::File;
use std::io::BufReader;
use std::io::Error;
use std::io::prelude::*;
use std::path::Path;

pub fn get_mnist_vector(fname: &str) -> Result<Vec<Vec<f32>>, Error> {
    let path = Path::new(fname);
    match File::open(&path) {
        Ok(file) => {
            let mut new_vec: Vec<Vec<f32>> = Vec::new();
            let reader = BufReader::new(file);
            for line in reader.lines() {
                match line {
                    Ok(s) => {
                        new_vec.push(mnist_test_to_vector(&s));
                    }
                    Err(reason) => {
                        return Err(reason);
                    }
                };
            }
            Ok(new_vec)
        }
        Err(reason) => Err(reason),
    }
}

pub fn mnist_test_to_vector(line: &str) -> Vec<f32> {
    line.trim().split(' ').map(|instr| instr.parse().unwrap()).collect()
}

pub mod xvecs {
    extern crate byteorder;
    use std::fs::File;
    use std::path::Path;
    use std::io::{BufReader, Result};
    use self::byteorder::{ReadBytesExt, LittleEndian};

    pub fn read_fvecs_file(fname: &str) -> Result<Vec<Vec<f32>>> {
        let path = Path::new(fname);
        let f = try!(File::open(path));
        let mut br = BufReader::new(f);
        let mut output_vec = Vec::new();
        while let Ok(i) = br.read_u32::<LittleEndian>() {
            let ind = i as usize;
            let mut line_vec = vec![0.0 as f32; ind];
            for j in 0..ind {
                line_vec[j] = try!(br.read_f32::<LittleEndian>());
            }
            output_vec.push(line_vec);
        }
        Ok(output_vec)
    }
    pub fn read_ivecs_file(fname: &str) -> Result<Vec<Vec<i32>>> {
        let path = Path::new(fname);
        let f = try!(File::open(path));
        let mut br = BufReader::new(f);
        let mut output_vec = Vec::new();
        while let Ok(i) = br.read_u32::<LittleEndian>() {
            let ind = i as usize;
            let mut line_vec = vec![0 as i32; ind];
            for j in 0..ind {
                line_vec[j] = try!(br.read_i32::<LittleEndian>());
            }
            output_vec.push(line_vec);
        }
        Ok(output_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::xvecs::*;
    #[test]
    fn test_read_mnist_file() {
        let fname: &str = "./mnist1k.dts";
        match get_mnist_vector(fname) {
            Ok(v) => println!("{:?}", v),
            Err(reason) => panic!("couldn't open because {}", reason),
        }
    }

    #[test]
    fn test_mnist_line_to_vector() {
        let my_vec = vec![0.0000, 1.0000, 0.00000, 0.20000];
        let my_string = "0.0000 1.0000 0.00000 0.20000";
        assert!(my_vec == mnist_test_to_vector(&my_string));
    }

    #[test]
    fn test_read_fvecs_file() {
        let y = read_fvecs_file("./sift_query.fvecs").unwrap();
        println!("{} vectors, length of first is {}", y.len(), y[0].len());
    }
}
