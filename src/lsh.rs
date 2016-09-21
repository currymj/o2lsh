use table;
use table::LSHTable;

// we want many lsh table

pub struct LSHLookup<'a, T: 'a, Q: 'a+?Sized> {
    tables: Vec<LSHTable<'a, T, Q>>,
    data: &'a [T]
}

impl<'a, T, Q> LSHLookup<'a, T, Q> {
    pub fn add_table(&mut self, new_table: LSHTable<'a, T, Q>) {
        self.tables.push(new_table);
    }

    pub fn new(data: &'a [T]) -> Self {
        LSHLookup {
            tables: Vec::new(),
            data: data
        }
    }
}
#[cfg(test)]
mod tests {
use super::*;
    use table::LSHTable;
    #[test]
    fn test_init_tables() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];
        let ms = vec![vec![1,2,3]];
        let mut mylookup = LSHLookup::new(&test_data);
        for i in 1..10 {
            let val = |q: &Vec<i32>| {0.0 as f64};
            let funcs = vec![Box::new(val)];
            let new_single_table = LSHTable::new(&test_data, funcs, &ms);
            mylookup.add_table(new_single_table);
        }
    }
}
