use super::measurement_like::MeasurementLike;
use super::pauli_set::PauliSet;
use super::{pauli_like::PauliLike, Pauli};
use crate::routines::f2_linalg::{row_echelon, Matrix};
use itertools::Itertools;
use rand::Rng;
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct IsometryTableau {
    pub n: usize,
    pub k: usize,
    pub logicals: PauliSet,
    pub stabilizers: PauliSet,
}

impl IsometryTableau {
    pub fn new(n: usize, k: usize) -> Self {
        let mut logicals = PauliSet::new(n + k);
        for i in 0..n {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[i] = true;
            logicals.insert_vec_bool(&vecbool, false);
        }
        for i in 0..n {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[i + n + k] = true;
            logicals.insert_vec_bool(&vecbool, false);
        }
        let mut stabilizers = PauliSet::new(n + k);
        for i in 0..k {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[n + k + n + i] = true;
            stabilizers.insert_vec_bool(&vecbool, false);
        }
        IsometryTableau {
            n,
            k,
            logicals,
            stabilizers,
        }
    }
    pub fn random(n: usize, k: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut iso = Self::new(n, k);
        for _ in 0..(n + k) * (n + k) {
            let i = rng.gen::<usize>() % (n + k);
            loop {
                let j = rng.gen::<usize>() % (n + k);
                if i == j {
                    continue;
                }
                iso.cnot(i, j);
                break;
            }
            for _ in 0..(n + k) {
                let gate_i = rng.gen::<u8>() % 3;
                if gate_i == 1 {
                    let q = rng.gen::<usize>() % (n + k);
                    iso.h(q);
                }
                if gate_i == 2 {
                    let q = rng.gen::<usize>() % (n + k);
                    iso.s(q);
                }
            }
        }
        iso
    }
    /// Put the full Tableau in column echelon form
    /// Warning: this method scratches the phases
    pub fn normalize_inplace(&mut self) {
        let mut table: Matrix = Vec::new();
        // The first k rows are the stabilizers
        for i in 0..self.k {
            let (_, vec) = self.stabilizers.get_as_vec_bool(i);
            table.push(vec);
        }
        // The next 2n rows are the logical operators
        for i in 0..2 * self.n {
            let (_, vec) = self.logicals.get_as_vec_bool(i);
            table.push(vec);
        }
        row_echelon(&mut table, self.k);
        let mut stabs = PauliSet::new(self.n + self.k);
        for row in table.iter().take(self.k) {
            stabs.insert_vec_bool(row, false);
        }
        let mut logicals = PauliSet::new(self.n + self.k);
        for i in 0..2 * self.n {
            logicals.insert_vec_bool(&table[self.k + i], false);
        }
        self.logicals = logicals;
        self.stabilizers = stabs;
    }
}

impl PauliLike for IsometryTableau {
    fn h(&mut self, i: usize) {
        self.logicals.h(i);
        self.stabilizers.h(i);
    }

    fn s(&mut self, i: usize) {
        self.logicals.s(i);
        self.stabilizers.s(i);
    }

    fn sd(&mut self, i: usize) {
        self.logicals.sd(i);
        self.stabilizers.sd(i);
    }

    fn sqrt_x(&mut self, i: usize) {
        self.logicals.sqrt_x(i);
        self.stabilizers.sqrt_x(i);
    }

    fn sqrt_xd(&mut self, i: usize) {
        self.logicals.sqrt_xd(i);
        self.stabilizers.sqrt_xd(i);
    }

    fn cnot(&mut self, i: usize, j: usize) {
        self.logicals.cnot(i, j);
        self.stabilizers.cnot(i, j);
    }
}

impl MeasurementLike for IsometryTableau {
    fn measure(&mut self, meas_string: Pauli) {
        // Find an anti-commuting stabilizer, if there isn't, you measured data!
        // Todo: generalize to allow measurements of data
        let (corr_i, _) = (0..self.stabilizers.len())
            .find_position(|i| !self.stabilizers.get_as_pauli(*i).commutes(&meas_string))
            .unwrap();

        let correction = self.stabilizers.get_as_pauli(corr_i);
        for j in corr_i + 1..self.stabilizers.len() {
            let other = self.stabilizers.get_as_pauli(j);
            if !other.commutes(&meas_string) {
                self.stabilizers.mul_pauli(j, &correction);
            }
        }

        // Now correct logicals
        for j in 0..self.logicals.len() {
            let other = self.logicals.get_as_pauli(j);
            if !other.commutes(&meas_string) {
                self.logicals.mul_pauli(j, &correction);
            }
        }
    }
}

impl fmt::Display for IsometryTableau {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Logicals:\n{}Stabilizers:\n{}",
            self.logicals, self.stabilizers
        )?;
        fmt::Result::Ok(())
    }
}
