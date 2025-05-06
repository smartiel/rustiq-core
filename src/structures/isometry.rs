use super::measurement_like::MeasurementLike;
use super::pauli_set::PauliSet;
use super::{pauli_like::PauliLike, Pauli};
use crate::routines::f2_linalg::{row_echelon, Matrix};
use itertools::Itertools;
use rand::Rng;
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct IsometryTableau {
    pub n: usize, // Number of logical operators
    pub k: usize, // Number of stabilizers
    pub logicals: PauliSet,
    pub stabilizers: PauliSet,
}

impl IsometryTableau {
    // Construct an IsometryTableau with n logical qubits and k ancilla in |0>.
    pub fn new(n: usize, k: usize) -> Self {
        let mut logicals = PauliSet::new(n + k);
        for i in 0..n {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[i] = true; // Logical operator X_i -> X_i
            logicals.insert_vec_bool(&vecbool, false);
        }
        for i in 0..n {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[i + n + k] = true; // Logical operator Z_i -> Z_i
            logicals.insert_vec_bool(&vecbool, false);
        }
        let mut stabilizers = PauliSet::new(n + k);
        for i in 0..k {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[n + k + n + i] = true; // stabilizer Z_{n+i}
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
    fn measure(&mut self, basis: &Pauli) {
        // Find an anti-commuting stabilizer, if there isn't, you measured data!
        // Todo: generalize to allow measurements of data
        let corr_opt = (0..self.stabilizers.len())
            .find_position(|i| !self.stabilizers.get_as_pauli(*i).commutes(basis));

        if let Some((corr_i, _)) = corr_opt {
            let correction = self.stabilizers.get_as_pauli(corr_i);
            self.stabilizers.set_pauli(corr_i, basis);
            for j in corr_i + 1..self.stabilizers.len() {
                let other = self.stabilizers.get_as_pauli(j);
                if !other.commutes(basis) {
                    self.stabilizers.mul_pauli(j, &correction);
                }
            }

            // Now correct logicals
            for j in 0..self.logicals.len() {
                let other = self.logicals.get_as_pauli(j);
                if !other.commutes(basis) {
                    self.logicals.mul_pauli(j, &correction);
                }
            }
        } else {
            panic!("No anti-commuting stabilizer, so you are measuring data. This has not been implemented.")
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_measure_x2() {
        let mut tableau = IsometryTableau::new(1, 1);
        let no_qubit = vec![0];
        let qubit2 = vec![2];
        let x2 = Pauli {
            n: 2,
            x_paulis: qubit2.clone(),
            z_paulis: no_qubit.clone(),
            sign: false,
        };
        let z2 = Pauli {
            n: 2,
            x_paulis: no_qubit,
            z_paulis: qubit2,
            sign: false,
        };
        tableau.measure(&x2);

        let mut ref_tableau = IsometryTableau::new(1, 1);
        ref_tableau.stabilizers.set_pauli(0, &x2);
        assert_eq!(tableau, ref_tableau);
    }

    #[test]
    fn test_measure_x1x2() {
        let mut tableau = IsometryTableau::new(1, 1);
        let no_qubit = vec![0];
        let qubit12 = vec![3];
        let x12 = Pauli {
            n: 2,
            x_paulis: qubit12.clone(),
            z_paulis: no_qubit.clone(),
            sign: false,
        };
        let z12 = Pauli {
            n: 2,
            x_paulis: no_qubit,
            z_paulis: qubit12,
            sign: false,
        };
        tableau.measure(&x12);

        let mut ref_tableau = IsometryTableau::new(1, 1);
        ref_tableau.stabilizers.set_pauli(0, &x12);
        ref_tableau.logicals.set_pauli(1, &z12);
        assert_eq!(tableau, ref_tableau);
    }
}
