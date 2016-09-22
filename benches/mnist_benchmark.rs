#![feature(test)]
extern crate o2lsh;
extern crate test;
use test::Bencher;


#[bench]
fn build_table_from_mnist_and_query(b: &mut Bencher) {

    let mnist_data = match o2lsh::util::get_mnist_vector("mnist1k.dts") {
        Ok(v) => v,
        Err(reason) => panic!("couldnt open because {}", reason)
    };

    let vec_length = mnist_data[0].len();

    let num_hashes = 10;
    let mnist_q = match o2lsh::util::get_mnist_vector("mnist1k.q") {
        Ok(v) => v,
        Err(reason) => panic!("couldnt open because {}", reason)
    };
    let zjs = o2lsh::multi::get_expected_zj_vals(num_hashes, 1.0);
    let sets: Vec<o2lsh::multi::PerturbationSet> = o2lsh::multi::gen_perturbation_sets(&zjs)
        .take(5)
        .collect();
    let ms: Vec<Vec<usize>> = sets.into_iter()
        .map(|x| {x.data})
        .collect();

    let hash_boxes: Vec<_> = (1..num_hashes).map(|_| o2lsh::hashes::get_hash_closure(vec_length, 1.0)).collect();

    let hash_table: o2lsh::table::LSHTable<Vec<f64>, Fn(&Vec<f64>) -> f64> = o2lsh::table::LSHTable::new_build(&mnist_data, hash_boxes, &ms);

    let run_test = || {
        for vec in mnist_q.iter() {
            hash_table.query_multiprobe(&vec);
        }
    };

    b.iter(run_test);
}
#[bench]
fn build_many_tables_from_mnist_and_time_query(b: &mut Bencher) {
    let mnist_data = match o2lsh::util::get_mnist_vector("mnist1k.dts") {
        Ok(v) => v,
        Err(reason) => panic!("couldnt open because {}", reason)
    };

    let vec_length = mnist_data[0].len();

    let num_hashes = 10;
    let mnist_q = match o2lsh::util::get_mnist_vector("mnist1k.q") {
        Ok(v) => v,
        Err(reason) => panic!("couldnt open because {}", reason)
    };
    let zjs = o2lsh::multi::get_expected_zj_vals(num_hashes, 1.0);
    let sets: Vec<o2lsh::multi::PerturbationSet> = o2lsh::multi::gen_perturbation_sets(&zjs)
        .take(5)
        .collect();
    let ms: Vec<Vec<usize>> = sets.into_iter()
        .map(|x| {x.data})
        .collect();

    let mut all_tables = o2lsh::lsh::LSHLookup::new(&mnist_data);
    for _ in 1..10 {
        let hash_boxes: Vec<_> = (1..num_hashes).map(|_| o2lsh::hashes::get_hash_closure(vec_length, 1.0)).collect();

        let hash_table: o2lsh::table::LSHTable<Vec<f64>, Fn(&Vec<f64>) -> f64> = o2lsh::table::LSHTable::new_build(&mnist_data, hash_boxes, &ms);
        all_tables.add_table(hash_table);

    }
    let run_test = || {
        for vec in mnist_q.iter() {
            all_tables.query_vec(&vec);
        }
    };

    b.iter(run_test);
}
