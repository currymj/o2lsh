extern crate rayon;
use std::collections::BTreeSet;
use self::rayon::prelude::*;
use super::lshtable::LSHTable;

use std::marker::PhantomData;

/// Represents a collection of many redundant LSH tables.
pub struct LSHLookup<'a, T: 'a, O:, L: LSHTable<'a, T, O>> {
    _m1: PhantomData<&'a T>,
    _m2: PhantomData<O>,
    tables: Vec<L>
}

impl<'a, T: Sync, L: LSHTable<'a, T, f32>>  LSHLookup<'a, T, f32, L>  {

    /// Add a new table to the collection.
    pub fn add_table(&mut self, new_table: L) {
        self.tables.push(new_table);
    }

    pub fn new()-> Self {
        LSHLookup {
            _m1: Default::default(),
            _m2: Default::default(),
            tables: Vec::new()
        }
    }


    /// Query all tables, returning deduplicated elements.
    pub fn query_vec(&self, v: &'a T, multiprobe_limit: usize) -> Vec<usize> {
        let mut output_set = BTreeSet::new();
        for table in &self.tables {
            for &reference in &table.query_multiprobe(v, multiprobe_limit) {
                output_set.insert(reference);
            }
        }
        let mut result = Vec::new();
        result.extend(output_set.into_iter());
        result
    }

    /// Query in parallel using Rayon.
    pub fn pquery_vec(&self, v: &'a T, multiprobe_limit: usize) -> Vec<usize> {
        let all_result_sets = self.tables.par_iter()
            .map(|table| {
                let mut output_set = BTreeSet::new();
                for &reference in &table.query_multiprobe(v, multiprobe_limit) {
                    output_set.insert(reference);
                }
                output_set
            })
            .reduce(|| {BTreeSet::new()}, |mut set1, set2| {
                set1.extend(set2);
                set1
            });


        let mut result = Vec::new();
        result.extend(all_result_sets.into_iter());
        result
    }

}
#[cfg(test)]
mod tests {
use super::*;
    use multi;
    use super::super::table::StandardLSHTable;
    #[test]
    fn test_init_tables() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];
        let zjs = multi::get_expected_zj_vals(1,1.0);
        let sets: Vec<multi::PerturbationSet> = multi::gen_perturbation_sets(&zjs)
            .take(5)
            .collect();
        let ms: Vec<Vec<usize>> = sets.into_iter()
            .map(|x| {x.data})
            .collect();
        let mut mylookup = LSHLookup::new();
        for _ in 1..10 {
            let val = |_: &Vec<i32>| {0.0 as f32};
            let funcs: Vec<Box<Fn(&Vec<i32>) -> f32 + Send + Sync>> = vec![Box::new(val)];
            let new_single_table = StandardLSHTable::new(&test_data, funcs, &ms);
            mylookup.add_table(new_single_table);
        }
    }
    fn is_send<T: Send>(_: &T) {}

    #[test]
    fn test_query_tables() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];
        let zjs = multi::get_expected_zj_vals(1,1.0);
        let sets: Vec<multi::PerturbationSet> = multi::gen_perturbation_sets(&zjs)
            .take(5)
            .collect();
        let ms: Vec<Vec<usize>> = sets.into_iter()
            .map(|x| {x.data})
            .collect();
        let mut mylookup = LSHLookup::new();
        for _ in 1..10 {
            let val = |_: &Vec<i32>| {0.0 as f32};
            let funcs: Vec<Box<Fn(&Vec<i32>) -> f32 + Send + Sync>> = vec![Box::new(val)];
            let new_single_table = StandardLSHTable::new(&test_data, funcs, &ms);
            is_send(&new_single_table);
            mylookup.add_table(new_single_table);
        }

        for data in &test_data {
            println!("{:?}", mylookup.query_vec(data,3));
        }
    }
}

