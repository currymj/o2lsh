extern crate o2lsh;


fn main() {
    let sift_data = o2lsh::util::xvecs::read_fvecs_file("sift/sift_base.fvecs")
        .expect("base file failed to open");
    let vec_length = sift_data[0].len();

    let num_hashes = 10;
    let sift_q = o2lsh::util::xvecs::read_fvecs_file("sift/sift_query.fvecs")
        .expect("query file failed to open");
    let zjs = o2lsh::multi::get_expected_zj_vals(num_hashes, 1.0);
    let sets: Vec<o2lsh::multi::PerturbationSet> = o2lsh::multi::gen_perturbation_sets(&zjs)
        .take(5)
        .collect();
    let ms: Vec<Vec<usize>> = sets.into_iter()
        .map(|x| {x.data})
        .collect();

    let mut all_tables  = o2lsh::lsh::LSHLookup::new();
    let hash_boxes: Vec<Box<Fn(&Vec<f32>) -> f32 + Send + Sync>> = (1..num_hashes).map(|_| o2lsh::hashes::get_hash_closure(vec_length, 10.0)).collect();
    all_tables.add_table(  o2lsh::table::StandardLSHTable::new_build(&sift_data, hash_boxes, &ms));

    let ground_truth = o2lsh::util::xvecs::read_ivecs_file("sift/sift_groundtruth.ivecs")
        .expect("failed to read ground truth");

    for (i, vec) in sift_q.iter().enumerate() {
        let vecvec = all_tables.query_vec(vec, 5);
        if (vecvec.len() == 100) {
            println!("100");
        }
    }
}
