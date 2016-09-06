// vector of buckets, each bucket being a linked list
// two constant/defined hash functions for taking signatures to bucket chains, and getting the bucket out of the chain
// pass in closures as signature functions (pass them here or just compute elsewhere? field of the data structure?)

// bucket should just be a Vector<usize>

type Bucket = Vec<usize>;
struct LSHTable<T> {
    buckets_array: Vec<Bucket>,
    data: Vec<T> // the LSH table is going to own its data. we should, though, make it possible to borrow refs
}

impl<T> LSHTable<T> {
    pub fn new(in_data: Vec<T>) -> Self {
        LSHTable {
            buckets_array: vec![vec![0 as usize; 1]; in_data.len()],
            data: in_data
        }
    }
}


#[cfg(test)]
mod tests {
    use super::LSHTable;
    #[test]
    fn test_it_works() {
        let test_data = vec![
            vec![1,2,3,4,5],
            vec![0,0,0,0,0],
            vec![1,2,3,4,5]
        ];

        let x = LSHTable::new(test_data);
    }
}

