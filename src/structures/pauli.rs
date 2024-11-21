use super::pauli_like::PauliLike;
use std::ops;

#[derive(Debug, PartialEq, Clone)]
pub struct Pauli {
    pub n: usize,
    pub x_paulis: Vec<u64>,
    pub z_paulis: Vec<u64>,
    pub sign: bool,
}

impl Pauli {
    pub fn new(n: usize) -> Self {
        let nr_blocks = n.div_ceil(64);
        Pauli {
            n,
            x_paulis: vec![0; nr_blocks],
            z_paulis: vec![0; nr_blocks],
            sign: false,
        }
    }

    pub fn commutes(&self, other: &Pauli) -> bool {
        if self.n != other.n {
            panic!("Can't compare two Paulis on different number of qubits");
        }

        let self_xz = self.x_paulis.iter().chain(self.z_paulis.iter());
        let other_zx = other.z_paulis.iter().chain(other.x_paulis.iter());
        let added = self_xz
            .zip(other_zx)
            .map(|(x1, z2)| x1 & z2)
            // This is equivalent to x.count_ones() + y.count_ones() % 2
            .reduce(|x, y| x ^ y)
            .unwrap_or(0);
        added.count_ones() % 2 == 0
    }

    pub fn x_bit(&self, idx: usize) -> bool {
        self.x_paulis[Pauli::stride(idx)] >> Pauli::offset(idx) & 1 != 0
    }

    pub fn z_bit(&self, idx: usize) -> bool {
        self.z_paulis[Pauli::stride(idx)] >> Pauli::offset(idx) & 1 != 0
    }

    fn stride(index: usize) -> usize {
        index / 64
    }
    fn offset(index: usize) -> usize {
        index % 64
    }
}

impl ops::Mul<&Pauli> for &Pauli {
    type Output = Pauli;

    fn mul(self, rhs: &Pauli) -> Self::Output {
        assert_eq!(self.n, rhs.n);
        let mut output = Pauli::new(self.n);
        output.sign = self.sign ^ rhs.sign;

        // X^x Z^z X^x' Z^z' = (-1)^(zx') X^(x+x') Z^(z+z')
        // Compute whether a -1 sign is applied
        let commute_phase = self
            .z_paulis
            .iter()
            .zip(rhs.x_paulis.iter())
            .map(|(z1, x2)| z1 & x2)
            .reduce(|x, y| x ^ y)
            .unwrap_or(0)
            .count_ones()
            % 2;
        output.sign ^= commute_phase != 0;

        output.x_paulis = self
            .x_paulis
            .iter()
            .zip(rhs.x_paulis.iter())
            .map(|(x1, x2)| x1 ^ x2)
            .collect();

        output.z_paulis = self
            .z_paulis
            .iter()
            .zip(rhs.z_paulis.iter())
            .map(|(x1, x2)| x1 ^ x2)
            .collect();
        output
    }
}

impl From<(&[bool], bool)> for Pauli {
    fn from((data, sign): (&[bool], bool)) -> Self {
        let n = data.len() / 2;
        let mut me = Pauli::new(n);
        me.sign = sign;

        let (x_bools, z_bools) = data.split_at(n);
        me.x_paulis = x_bools
            .chunks(64)
            .map(|chunk| {
                let mut sum = 0;
                // Encode in LSB
                for b in chunk.iter().rev() {
                    sum = (sum << 1) | (*b as u64);
                }
                sum
            })
            .collect();

        me.z_paulis = z_bools
            .chunks(64)
            .map(|chunk| {
                let mut sum = 0;
                // Encode in LSB
                for b in chunk.iter().rev() {
                    sum = (sum << 1) | (*b as u64);
                }
                sum
            })
            .collect();

        me
    }
}

// TODO
// impl PauliLike for Pauli {
//     fn h(&mut self, i: usize) {
//         self.data.swap(i, i + self.n);
//     }

//     fn s(&mut self, i: usize) {
//         self.data[i + self.n] ^= self.data[i];
//     }

//     fn sd(&mut self, i: usize) {
//         self.data[i + self.n] ^= self.data[i];
//     }

//     fn sqrt_x(&mut self, i: usize) {
//         self.data[i] ^= self.data[i + self.n];
//     }

//     fn sqrt_xd(&mut self, i: usize) {
//         self.data[i] ^= self.data[i + self.n];
//     }

//     fn cnot(&mut self, i: usize, j: usize) {
//         self.data[i + self.n] ^= self.data[j + self.n];
//         self.data[j] ^= self.data[i];
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_bool_vec() {
        let mut x = Pauli::new(2);
        x.x_paulis[0] = 3;
        let vector = vec![true, true, false, false];
        let p = Pauli::from((vector.as_slice(), false));

        assert_eq!(x, p)
    }

    #[test]
    fn test_get_bits() {
        let mut x = Pauli::new(2);
        x.x_paulis[0] = 3;
        x.z_paulis[0] = 1;
        assert!(x.x_bit(0));
        assert!(x.x_bit(1));
        assert!(x.z_bit(0));
        assert!(!x.z_bit(1));
    }
}
