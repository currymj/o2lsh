use table::LSHTable;
use std::collections::BTreeSet;
// we want many lsh table

pub struct LSHLookup<'a, T: 'a, Q: 'a+?Sized> {
    tables: Vec<LSHTable<'a, T, Q>>,
}

impl<'a, T, Q: 'a+?Sized> LSHLookup<'a, T, Q> where Q: Fn(&'a T) -> f32 {
    pub fn add_table(&mut self, new_table: LSHTable<'a, T, Q>) {
        self.tables.push(new_table);
    }

    pub fn new()-> Self {
        LSHLookup {
            tables: Vec::new()
        }
    }


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
}
#[cfg(test)]
mod tests {
use super::*;
    use table::LSHTable;
    use multi;
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
            let funcs = vec![Box::new(val)];
            let new_single_table = LSHTable::new(&test_data, funcs, &ms);
            mylookup.add_table(new_single_table);
        }
    }

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
            let funcs = vec![Box::new(val)];
            let new_single_table = LSHTable::new(&test_data, funcs, &ms);
            mylookup.add_table(new_single_table);
        }

        for data in &test_data {
            println!("{:?}", mylookup.query_vec(data,3));
        }
    }
}
