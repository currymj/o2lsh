pub trait LSHTable<'a, T, O> : Sync + Send {
    fn query_multiprobe(&self, v: &'a T, multiprobe_limit: usize) -> Vec<usize>;
    fn query_vec(&self, v: &'a T) -> Vec<usize>;
}
