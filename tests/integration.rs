extern crate o2lsh;

#[test]
fn build_table_from_mnist() {
    let mnist_data = match o2lsh::util::get_mnist_vector("mnist1k.dts") {
        Ok(v) => v,
        Err(reason) => panic!("couldnt open because {}", reason)
    };

    let vec_length = mnist_data[0].len();

    let num_hashes = 10;
    let hash_boxes: Vec<_> = (1..num_hashes).map(|_| o2lsh::hashes::get_hash_closure(vec_length, 1.0)).collect();

    let hash_table: o2lsh::table::LSHTable<Vec<f64>, Fn(&Vec<f64>) -> f64> = o2lsh::table::LSHTable::new_build(&mnist_data, &hash_boxes);
}
