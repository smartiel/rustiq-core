use crate::routines::f2_linalg::*;
use crate::structures::{CliffordCircuit, CliffordGate, PauliSet};

pub fn build_table(pauli_set: &PauliSet) -> (Vec<Vec<bool>>, Vec<Vec<bool>>) {
    let n = pauli_set.n;
    let mut table_z = vec![vec![false; pauli_set.len()]; n];
    let mut table_x = vec![vec![false; pauli_set.len()]; n];
    for i in 0..pauli_set.len() {
        let (_, vec) = pauli_set.get_as_vec_bool(i);
        for j in 0..n {
            table_z[j][i] = vec[j + n];
            table_x[j][i] = vec[j];
        }
    }
    (table_z, table_x)
}

fn swap_rows(table: &mut Vec<Vec<bool>>, rows: &[usize]) {
    let mut new_table = Vec::new();
    for j in rows.iter() {
        new_table.push(table[*j].clone());
    }
    *table = new_table;
}

type FullRank = (
    CliffordCircuit,
    Vec<usize>,
    usize,
    (Vec<Vec<bool>>, Vec<Vec<bool>>),
);

pub fn make_full_rank(pauli_set: &PauliSet) -> FullRank {
    let (mut t_z, mut t_x) = build_table(pauli_set);
    let mut circuit = CliffordCircuit::new(pauli_set.n);
    for i in 0..pauli_set.n {
        let rk = f2_rank(&t_x);
        std::mem::swap(t_z.get_mut(i).unwrap(), t_x.get_mut(i).unwrap());
        let new_rk = f2_rank(&t_x);
        if new_rk > rk {
            circuit.gates.push(CliffordGate::H(i));
        } else {
            std::mem::swap(t_z.get_mut(i).unwrap(), t_x.get_mut(i).unwrap());
        }
    }
    let mut fr_rows = Vec::new();
    let mut wit = Vec::new();
    let mut rk = 0;
    for (i, item) in t_x.iter().enumerate().take(pauli_set.n) {
        wit.push(item.clone());
        if f2_rank(&wit) == rk + 1 {
            fr_rows.push(i);
            rk += 1;
        }
    }
    let rk = f2_rank(&t_x);

    for i in 0..pauli_set.n {
        if !fr_rows.contains(&i) {
            fr_rows.push(i);
        }
    }
    swap_rows(&mut t_x, &fr_rows);
    swap_rows(&mut t_z, &fr_rows);
    diagonalize(&mut t_x, &mut t_z, rk);
    (circuit, fr_rows, rk, (t_z, t_x))
}
pub fn permute_circuit(circuit: &CliffordCircuit, permutation: &[usize]) -> CliffordCircuit {
    let mut permuted_circuit = CliffordCircuit::new(circuit.nqbits);
    for gate in circuit.gates.iter() {
        match gate {
            CliffordGate::CNOT(i, j) => permuted_circuit
                .gates
                .push(CliffordGate::CNOT(permutation[*i], permutation[*j])),
            CliffordGate::CZ(i, j) => permuted_circuit
                .gates
                .push(CliffordGate::CZ(permutation[*i], permutation[*j])),
            CliffordGate::S(i) => permuted_circuit
                .gates
                .push(CliffordGate::S(permutation[*i])),
            CliffordGate::Sd(i) => permuted_circuit
                .gates
                .push(CliffordGate::Sd(permutation[*i])),
            _ => panic!("This should never happen: encountered gate {:?}", gate),
        }
    }
    permuted_circuit
}

#[cfg(test)]
mod common_codiag_tests {

    use super::*;
    use crate::structures::pauli_like::PauliLike;
    #[test]
    fn test_full_r() {
        let mut pset = PauliSet::new(4);
        pset.insert("ZZII", false);
        pset.insert("IIXI", false);
        pset.insert("IIIZ", false);
        let (c, m, _, _) = make_full_rank(&pset);
        println!("{:?}", m);
        println!("{:?}", c);
    }
    #[test]
    fn test_build_table() {
        let mut pset = PauliSet::new(4);
        pset.insert("ZZII", false);
        pset.insert("IXXI", false);
        pset.insert("IIZY", false);
        let (z, x) = build_table(&pset);
        assert_eq!(z[0], &[true, false, false]);
        assert_eq!(z[1], &[true, false, false]);
        assert_eq!(z[2], &[false, false, true]);
        assert_eq!(z[3], &[false, false, true]);

        assert_eq!(x[0], &[false, false, false]);
        assert_eq!(x[1], &[false, true, false]);
        assert_eq!(x[2], &[false, true, false]);
        assert_eq!(x[3], &[false, false, true]);
    }
    use rand::Rng;
    fn random_instance(n: usize, m: usize) -> PauliSet {
        let mut rng = rand::rng();
        let mut pset = PauliSet::new(n);
        for _ in 0..m {
            let mut vec: Vec<bool> = vec![false; 2 * n];
            for b in vec.iter_mut().take(n) {
                *b = rng.random::<bool>();
            }
            pset.insert_vec_bool(&vec, false);
        }
        for _ in 0..n * n {
            let i = rng.random::<u64>() as usize % n;
            loop {
                let j = rng.random::<u64>() as usize % n;
                if j != i {
                    pset.cnot(i, j);
                    let g2 = rng.random::<bool>();
                    if g2 {
                        pset.h(j);
                    } else {
                        pset.s(j);
                    }
                    break;
                }
            }
            let g1 = rng.random::<bool>();
            if g1 {
                pset.h(i);
            } else {
                pset.s(i);
            }
        }
        pset
    }
    #[test]
    fn test_random_thin() {
        for _ in 0..10 {
            let pset = random_instance(50, 20);
            let (_, _, rk, (_, x)) = make_full_rank(&pset);
            for i in 0..rk {
                let mut v = vec![false; pset.len()];
                v[i] = true;
                assert_eq!(x[i], v);
            }
        }
    }
    #[test]
    fn test_random_thick() {
        for _ in 0..10 {
            let pset = random_instance(20, 50);
            let (_, _, rk, (_, x)) = make_full_rank(&pset);
            for i in 0..rk {
                let mut v = vec![false; pset.len()];
                v[i] = true;
                assert_eq!(x[i], v);
            }
        }
    }
}
