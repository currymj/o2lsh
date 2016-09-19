extern crate revord;
use std::cmp::Ordering;
use self::revord::RevOrd;
use std::collections::BinaryHeap;

fn bucket_distance(fi: f64, hi: f64, delta: i32, W: f64) -> f64 {
    if delta == -1 {
        W*(fi - hi)
    } else if delta == 1 {
        W - W*(fi - hi)
    } else {
        0 as f64
    }
}
fn compute_pi_j<T>(q: &[T], f_sig: &[f64], h_sig: &[f64], W: f64) -> Vec<(usize, i32)> {
    let mut intermediate_vec: Vec<((usize,i32), f64)> = f_sig.iter()
        .zip(h_sig.iter())
        .enumerate()
        .flat_map(|(i, (fi, hi))| {
            vec![((i, 1), bucket_distance(*fi, *hi, 1, W)),
                 ((i, -1), bucket_distance(*fi, *hi, -1, W))].into_iter()
        }).collect();
    intermediate_vec.sort_by(|a, b| {
        if a.1 > b.1 {
            Ordering::Greater
        } else if a.1 < b.1 {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    });
    intermediate_vec.iter().map(|a| {a.0}).collect()
}

#[test]
fn sorted_delta_test() {
    let test_q = vec![1.0,2.0,3.0,4.0,5.0];
    let f_sig = vec![1.5,1.2,2.2];
    let h_sig = vec![1.0,1.0,2.0];
    let W = 10.0;
    compute_pi_j(&test_q, &f_sig, &h_sig, W);
}

fn score_set(perturbation_set: &[usize], square_zj_list: &[f64]) -> f64 {
    perturbation_set.iter().map(|ind| {square_zj_list[*ind]}).sum()
}

fn gen_perturbation_sets<'b>(zj_l: &'b[f64]) -> PerturbationIterator<'b> {
    let mut perturb_return: PerturbationIterator<'b> = PerturbationIterator {
        heap: BinaryHeap::new(),
        zj_list: zj_l
    };
    let zero_vec = vec![0];
    let a0: PerturbationSet<'b> = PerturbationSet {
        data: zero_vec,
        zj_list: zj_l,
        max_m: perturb_return.zj_list.len() / 2
    };
    perturb_return.heap.push(RevOrd(a0));
    perturb_return
}

struct PerturbationIterator<'a> {
    heap: BinaryHeap<RevOrd<PerturbationSet<'a>>>,
    zj_list: &'a[f64]
}

#[test]
fn test_perturbation_iterator() {
    let zjl = vec![1.0,2.0,3.0,4.0];
    let ii = gen_perturbation_sets(&zjl);
    let z: Vec<_> = ii.take(3).collect();
}
impl<'a> Iterator for PerturbationIterator<'a> {
    type Item = PerturbationSet<'a>;
    fn next(&mut self) -> Option<PerturbationSet<'a>> {
        loop {
            match self.heap.pop() {
                Some(revord_val) => {
                    let a_i = revord_val.0;
                    let a_s = a_i.shift();
                    if a_s.valid() {
                        self.heap.push(RevOrd(a_s));
                    }
                    let a_e = a_i.expand();
                    if a_e.valid() {
                        self.heap.push(RevOrd(a_e));
                    }

                    if a_i.valid() {
                        return Some(a_i);
                    }
                },
                None => return None
            };
        }
    }
}


#[derive(PartialEq, Debug)]
struct PerturbationSet<'a> {
    data: Vec<usize>,
    zj_list: &'a [f64],
    max_m:  usize
}
impl<'a> Eq for PerturbationSet<'a> {}

impl<'a> Ord for PerturbationSet<'a> {
    fn cmp(&self, other: &PerturbationSet) -> Ordering {
        let self_score = score_set(&(self.data), self.zj_list);
        let other_score = score_set(&(other.data), other.zj_list);
        if (self_score - other_score) <= 1e-6 {
            Ordering::Equal
        } else if self_score > other_score {
            Ordering::Greater
        } else {
            Ordering::Less
        }
    }
}

impl<'a> PartialOrd for PerturbationSet<'a> {
    fn partial_cmp(&self, other: &PerturbationSet) -> Option<Ordering> {
        let self_score = score_set(&(self.data), self.zj_list);
        let other_score = score_set(&(other.data), other.zj_list);
        self_score.partial_cmp(&other_score)
    }
}
impl<'a> PerturbationSet<'a> {
    pub fn len(&self) -> usize {
        self.data.len()
    }
    pub fn shift<'b>(&'b self) -> PerturbationSet<'a> {
        let mut new_data = Vec::new();
        let max_val = self.data.iter().max().unwrap();
        let newlist = self.zj_list;
        for &x in &self.data {
            let y = x;
            if y == *max_val {
                new_data.push(y+1);
            } else {
                new_data.push(y);
            }
        }
        PerturbationSet {
            data: new_data,
            zj_list: newlist,
            max_m: self.max_m
        }
    }
    pub fn valid(&self) -> bool {
        let &max_val = self.data.iter().max().unwrap();
        if max_val >= self.max_m {
            return false;
        }
        for &val in &self.data {
            let other_should_be_missing = 2 * self.max_m + 1 - val;
            if let Some(_) = self.data.iter().position(|&y| { y == other_should_be_missing}) {
                return false;
            }
        }
        return true;
    }
    pub fn expand<'b>(&'b self) -> PerturbationSet<'a> {
        let mut new_data = Vec::new();
        let max_val = self.data.iter().max().unwrap();
        for &x in &self.data {
            new_data.push(x);
        }
        new_data.push(*max_val + 1);
        PerturbationSet {
            data: new_data,
            zj_list: self.zj_list,
            max_m: self.max_m
        }
    }
}

fn expected_zj_squared(j: usize, M: usize, W: f64) -> f64 {
    if j <= M {
        (((j * (j + 1)) as f64)
            / ((4 * (M + 1) * (M + 2)) as f64)) * (W * W)
    } else {
        (W*W) * (1.0 - ((2 * M + 1 - j) as f64) / ((M + 1) as f64) +
                 (((2 * M + 1 - j) * (2 * M + 2 - j)) as f64) / ((4 * (M + 1) * (M + 2)) as f64))
    }
}

/*fn get_perturbation_iterator<'b>(j: usize, M: usize, W: f64) -> PerturbationIterator<'b> {
    let mut zj_vals = Vec::new();
    for j in 1..20 {
        zj_vals.push(expected_zj_squared(j, 10, 1.0));
    }
    gen_perturbation_sets(&zj_vals)
}*/

#[test]
fn test_expected_zj_squared() {
    expected_zj_squared(1,10,1.0);
}

#[test]
fn test_perturbation_iterator_zj() {
    let mut zj_vals = Vec::new();
    for j in 1..20 {
        zj_vals.push(expected_zj_squared(j, 10, 1.0));
    }
    let perturbation_iterator = gen_perturbation_sets(&zj_vals);
    for pert in perturbation_iterator {
        println!("{}", pert.len());
    }
}
