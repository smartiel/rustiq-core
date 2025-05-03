use super::pauli::Pauli;
use super::pauli_like::PauliLike;
use crate::synthesis::pauli_network::chunks::CHUNK_CONJUGATION_SCORE;
use itertools::izip;
use std::cmp::max;
use std::fmt;

const WIDTH: usize = 64;

fn get_stride(index: usize) -> usize {
    index / WIDTH
}

fn get_offset(index: usize) -> usize {
    index % WIDTH
}

/// A set of Pauli operators (module global phase)
/// Conjugation by Clifford gates are vectorized
#[derive(Clone, Debug, PartialEq)]
pub struct PauliSet {
    pub n: usize,
    nstrides: usize,
    noperators: usize,
    start_offset: usize,
    /// The X and Z parts of the Pauli operators (in row major)
    /// The X part spans the first `n` rows and the Z part spans the last `n` rows
    data_array: Vec<Vec<u64>>,
    phases: Vec<u64>,
}

impl PauliSet {
    /// Allocate an empty set of n-qubit Pauli operators
    pub fn new(n: usize) -> Self {
        Self {
            n,
            nstrides: 0,
            noperators: 0,
            start_offset: 0,
            data_array: vec![Vec::new(); 2 * n],
            phases: Vec::new(),
        }
    }
    /// Allocate a set of m n-qubit Pauli operators set to the identity
    pub fn new_empty(n: usize, m: usize) -> Self {
        let nstrides = get_stride(m) + 1;
        Self {
            n,
            nstrides,
            noperators: m,
            start_offset: 0,
            data_array: vec![vec![0; nstrides]; 2 * n],
            phases: vec![0; nstrides],
        }
    }
    // Construction from a list of operators
    pub fn from_slice(data: &[String]) -> Self {
        if data.is_empty() {
            return Self::new(0);
        }
        let n = data.first().unwrap().len();
        let mut pset = Self::new(n);
        for piece in data {
            pset.insert(piece, false);
        }
        pset
    }
    /// Returns the number of operators stored in the set
    pub fn len(&self) -> usize {
        self.noperators
    }

    pub fn is_empty(&self) -> bool {
        self.noperators == 0
    }

    /// Inserts a new Pauli operator in the set and returns its index
    pub fn insert(&mut self, axis: &str, phase: bool) -> usize {
        let stride = get_stride(self.noperators + self.start_offset);
        let offset = get_offset(self.noperators + self.start_offset);
        if stride == self.nstrides {
            self.nstrides += 1;
            self.data_array.iter_mut().for_each(|row| row.push(0));
            self.phases.push(0);
        }
        // Setting the phase
        if phase {
            self.phases[stride] |= 1 << offset;
        }
        // Setting the operator
        for (index, pauli) in axis.chars().enumerate() {
            match pauli {
                'Z' => self.data_array[index + self.n][stride] |= 1 << offset,
                'X' => self.data_array[index][stride] |= 1 << offset,
                'Y' => {
                    self.data_array[index][stride] |= 1 << offset;
                    self.data_array[index + self.n][stride] |= 1 << offset
                }
                _ => {}
            }
        }
        self.noperators += 1;
        self.noperators - 1
    }

    /// Inserts a new Pauli operator described as a vector of bool in the set and returns its index
    pub fn insert_vec_bool(&mut self, axis: &[bool], phase: bool) -> usize {
        let stride = get_stride(self.noperators + self.start_offset);
        let offset = get_offset(self.noperators + self.start_offset);
        if stride == self.nstrides {
            self.nstrides += 1;
            self.data_array.iter_mut().for_each(|row| row.push(0));
            self.phases.push(0);
        }
        if phase {
            self.phases[stride] |= 1 << offset;
        }
        for (index, value) in axis.iter().enumerate() {
            if *value {
                self.data_array[index][stride] |= 1 << offset;
            }
        }
        self.noperators += 1;
        self.noperators - 1
    }
    pub fn insert_pauli(&mut self, pauli: &Pauli) -> usize {
        self.insert_vec_bool(&pauli.data, pauli.phase == 2)
    }
    pub fn set_phase(&mut self, col: usize, phase: bool) {
        let stride = get_stride(col);
        let offset = get_offset(col);
        if phase != (((self.phases[stride] >> offset) & 1) != 0) {
            self.phases[stride] ^= 1 << offset;
        }
    }

    pub fn set_entry(&mut self, operator_index: usize, qbit: usize, x_part: bool, z_part: bool) {
        let stride = get_stride(operator_index + self.start_offset);
        let offset = get_offset(operator_index + self.start_offset);
        if x_part != (1 == (self.data_array[qbit][stride] >> offset) & 1) {
            self.data_array[qbit][stride] ^= 1 << offset;
        }
        if z_part != (1 == (self.data_array[qbit + self.n][stride] >> offset) & 1) {
            self.data_array[qbit + self.n][stride] ^= 1 << offset;
        }
    }
    pub fn set_raw_entry(&mut self, row: usize, col: usize, value: bool) {
        let stride = get_stride(col);
        let offset = get_offset(col);
        if value != (1 == (self.data_array[row][stride] >> offset) & 1) {
            self.data_array[row][stride] ^= 1 << offset;
        }
    }

    /// Clears the data of the Pauli set
    pub fn clear(&mut self) {
        for j in 0..self.nstrides {
            for i in 0..2 * self.n {
                self.data_array[i][j] = 0;
            }
            self.phases[j] = 0;
        }
        self.noperators = 0;
        self.start_offset = 0;
    }
    /// Pops the first rotation in the set
    pub fn pop(&mut self) {
        let stride = get_stride(self.start_offset);
        let offset = get_offset(self.start_offset);
        for i in 0..2 * self.n {
            self.data_array[i][stride] &= !(1 << offset);
        }
        self.phases[stride] &= !(1 << offset);
        self.start_offset += 1;
        self.noperators -= 1;
    }
    /// Pops the last rotation in the set
    pub fn pop_last(&mut self) {
        let stride = get_stride(self.start_offset + self.noperators - 1);
        let offset = get_offset(self.start_offset + self.noperators - 1);
        for i in 0..2 * self.n {
            self.data_array[i][stride] &= !(1 << offset);
        }
        self.phases[stride] &= !(1 << offset);
        self.noperators -= 1;
    }
    /// Set some operator to identity (because popping in the middle is expensive :O)
    pub fn set_to_identity(&mut self, operator_index: usize) {
        // set_entry(&mut self, operator_index: usize, qbit: usize, x_part: bool, z_part: bool)
        for i in 0..self.n {
            self.set_entry(operator_index, i, false, false);
        }
    }
    /// Get the operator at index `operator_index` as a pair (phase, string)
    pub fn get(&self, operator_index: usize) -> (bool, String) {
        let operator_index = operator_index + self.start_offset;
        let mut output = String::new();
        let stride = get_stride(operator_index);
        let offset = get_offset(operator_index);
        for i in 0..self.n {
            match (
                (self.data_array[i][stride] >> offset) & 1,
                (self.data_array[i + self.n][stride] >> offset) & 1,
            ) {
                (1, 0) => {
                    output += "X";
                }
                (0, 1) => {
                    output += "Z";
                }
                (1, 1) => {
                    output += "Y";
                }
                _ => {
                    output += "I";
                }
            }
        }
        (((self.phases[stride] >> offset) & 1 != 0), output)
    }

    /// Get the operator at index `operator_index` as a pair `(bool, Vec<bool>)`
    pub fn get_as_vec_bool(&self, operator_index: usize) -> (bool, Vec<bool>) {
        let operator_index = operator_index + self.start_offset;
        let mut output = Vec::new();
        let stride = get_stride(operator_index);
        let offset = get_offset(operator_index);
        for i in 0..2 * self.n {
            output.push(((self.data_array[i][stride] >> offset) & 1) != 0);
        }
        (((self.phases[stride] >> offset) & 1 != 0), output)
    }

    /// Get the operator at index `operator_index` as a `Pauli` object
    pub fn get_as_pauli(&self, operator_index: usize) -> Pauli {
        let (phase, data) = self.get_as_vec_bool(operator_index);
        Pauli::from_vec_bool(data, if phase { 2 } else { 0 })
    }
    /// Get a single entry of the PauliSet
    pub fn get_entry(&self, row: usize, col: usize) -> bool {
        let col = col + self.start_offset;
        let stride = get_stride(col);
        let offset = get_offset(col);
        ((self.data_array[row][stride] >> offset) & 1) != 0
    }
    pub fn get_phase(&self, col: usize) -> bool {
        let col = col + self.start_offset;
        let stride = get_stride(col);
        let offset = get_offset(col);
        ((self.phases[stride] >> offset) & 1) != 0
    }

    pub fn get_i_factors(&self) -> Vec<u8> {
        let mut output = Vec::new();
        for i in 0..self.len() {
            let mut ifact: u8 = 0;
            for j in 0..self.n {
                if self.get_entry(j, i) & self.get_entry(j + self.n, i) {
                    ifact += 1;
                }
            }
            output.push(ifact % 4);
        }
        output
    }
    pub fn get_i_factors_single_col(&self, col: usize) -> u8 {
        let mut ifact: u8 = 0;
        for j in 0..self.n {
            if self.get_entry(j, col) & self.get_entry(j + self.n, col) {
                ifact += 1;
            }
        }
        ifact % 4
    }
    /// Get the inverse Z output of the tableau (assuming the PauliSet is a Tableau, i.e. has exactly 2n operators storing X1...Xn Z1...Zn images)
    pub fn get_inverse_z(&self, qbit: usize) -> (bool, String) {
        let mut pstring = String::new();
        for i in 0..self.n {
            let x_bit = self.get_entry(qbit, i + self.n);
            let z_bit = self.get_entry(qbit, i);
            match (x_bit, z_bit) {
                (false, false) => {
                    pstring.push('I');
                }
                (true, false) => {
                    pstring.push('X');
                }
                (false, true) => {
                    pstring.push('Z');
                }
                (true, true) => {
                    pstring.push('Y');
                }
            }
        }
        (self.get_phase(qbit + self.n), pstring)
    }
    /// Get the inverse X output of the tableau (assuming the PauliSet is a Tableau, i.e. has exactly 2n operators storing X1...Xn Z1...Zn images)
    pub fn get_inverse_x(&self, qbit: usize) -> (bool, String) {
        let mut pstring = String::new();
        let mut cy = 0;
        for i in 0..self.n {
            let x_bit = self.get_entry(qbit + self.n, i + self.n);
            let z_bit = self.get_entry(qbit + self.n, i);
            match (x_bit, z_bit) {
                (false, false) => {
                    pstring.push('I');
                }
                (true, false) => {
                    pstring.push('X');
                }
                (false, true) => {
                    pstring.push('Z');
                }
                (true, true) => {
                    pstring.push('Y');
                    cy += 1;
                }
            }
        }
        ((cy % 2 != 0), pstring)
    }

    /// Returns the sum mod 2 of the logical AND of a row with an external vector of booleans
    pub fn and_row_acc(&self, row: usize, vec: &[bool]) -> bool {
        let mut output = false;
        for (i, item) in vec.iter().enumerate().take(2 * self.n) {
            output ^= self.get_entry(row, i) & item;
        }
        output
    }

    /// Check equality between two operators
    pub fn equals(&self, i: usize, j: usize) -> bool {
        let (_, vec1) = self.get_as_vec_bool(i);
        let (_, vec2) = self.get_as_vec_bool(j);
        vec1 == vec2
    }

    /// Checks if two operators in the set commute
    pub fn commute(&self, i: usize, j: usize) -> bool {
        let mut count_diff = 0;
        for k in 0..self.n {
            if (self.get_entry(k, i) & self.get_entry(k + self.n, j))
                ^ (self.get_entry(k + self.n, i) & self.get_entry(k, j))
            {
                count_diff += 1;
            }
        }
        (count_diff % 2) == 0
    }

    // Returns the support of the operator (i.e. the list of qbit indices on which the operator acts non trivially)
    pub fn get_support(&self, index: usize) -> Vec<usize> {
        let mut support = Vec::new();
        let index = index + self.start_offset;
        let stride = get_stride(index);
        let offset = get_offset(index);
        for i in 0..self.n {
            if (((self.data_array[i][stride] | self.data_array[i + self.n][stride]) >> offset) & 1)
                != 0
            {
                support.push(i);
            }
        }
        support
    }

    /// Returns the support size of the operator (i.e. the number of non-I Pauli term in the operator)
    pub fn support_size(&self, index: usize) -> usize {
        let index = index + self.start_offset;
        let mut count = 0;
        let stride = get_stride(index);
        let offset = get_offset(index);
        for i in 0..self.n {
            if (((self.data_array[i][stride] | self.data_array[i + self.n][stride]) >> offset) & 1)
                != 0
            {
                count += 1;
            }
        }
        count
    }

    /*
           Internal methods
    */

    /// XORs row `i` into row `j`
    fn row_op(&mut self, i: usize, j: usize) {
        let (left, right) = self.data_array.split_at_mut(max(i, j));
        let (target_row, source_row) = if i < j {
            (right.get_mut(0).unwrap(), left.get(i).unwrap())
        } else {
            (left.get_mut(j).unwrap(), right.first().unwrap())
        };

        for (v1, v2) in source_row.iter().zip(target_row.iter_mut()) {
            *v2 ^= *v1;
        }
    }

    pub fn swap_qbits(&mut self, i: usize, j: usize) {
        self.data_array.swap(i, j);
        self.data_array.swap(self.n + i, self.n + j);
    }

    /// Offset the phases by the logical bitwise and of two target rows
    fn update_phase_and(&mut self, i: usize, j: usize) {
        for (v1, v2, phase) in izip!(
            self.data_array[i].iter(),
            self.data_array[j].iter(),
            self.phases.iter_mut()
        ) {
            *phase ^= *v1 & *v2;
        }
    }

    /// Same thing as `update_phase_and` but computes the and of 4 rows instead
    fn update_phase_and_many(&mut self, i: usize, j: usize, k: usize, l: usize) {
        for (v1, v2, v3, v4, phase) in izip!(
            self.data_array[i].iter(),
            self.data_array[j].iter(),
            self.data_array[k].iter(),
            self.data_array[l].iter(),
            self.phases.iter_mut()
        ) {
            *phase ^= *v1 & *v2 & *v3 & *v4;
        }
    }

    /*
       Gate conjugation
    */

    /*
       Metrics for synthesis algorithms
    */
    pub fn count_id(&self, qbit: usize) -> usize {
        let mut count: usize = 0;
        let nstart_stride = get_stride(self.start_offset);
        for ns in nstart_stride..self.nstrides {
            let r_x = self.data_array[qbit][ns];
            let r_z = self.data_array[qbit + self.n][ns];
            let value = r_x | r_z;
            count += value.trailing_zeros() as usize;
            if ns == nstart_stride {
                count -= self.start_offset;
            }
            if ns == self.nstrides - 1 && value == 0 {
                count -= 64 - get_offset(self.start_offset + self.noperators);
            }
            if value != 0 {
                return count;
            }
        }
        count
    }

    /// Returns `true` if pauli `col` for qubit `q` is `I`
    #[inline]
    pub fn is_i(&self, qbit: usize, col: usize) -> bool {
        !self.get_entry(qbit, col) && !self.get_entry(qbit + self.n, col)
    }

    #[inline]
    pub fn count_leading_i(&self, qbit: usize, order: &[usize]) -> usize {
        order
            .iter()
            .take_while(|&&col| self.is_i(qbit, col))
            .count()
    }

    /// Returns (the index of) the Pauli pair over the qubits `i` and `j`
    /// for the Pauli operator in column `col`.
    #[inline]
    pub fn pauli_pair_index(&self, i: usize, j: usize, col: usize) -> usize {
        let n = self.n;
        let s0 = self.get_entry(i, col);
        let d0 = self.get_entry(i + n, col);
        let s1 = self.get_entry(j, col);
        let d1 = self.get_entry(j + n, col);

        ((s0 as usize) << 3) | ((d0 as usize) << 2) | ((s1 as usize) << 1) | (d1 as usize)
    }

    /// Computes the score of conjugating the Pauli pair over the qubits `i` and `j`
    /// by chunk `c`. This is equivalent to the scoring function described in the
    /// paper but uses the precomputed table lookup instead of performing conjugation.
    #[inline]
    pub fn count_leading_i_conjugation(
        &self,
        i: usize,
        j: usize,
        q: usize,
        c: usize,
        order: &[usize],
    ) -> usize {
        order
            .iter()
            .take_while(|&&col| CHUNK_CONJUGATION_SCORE[c][q][self.pauli_pair_index(i, j, col)] > 0)
            .count()
    }

    /// Sorts the set by support size
    pub fn support_size_sort(&mut self) {
        // We first build the "transpose" of data_array (cheaper this way)
        let mut transposed: Vec<(bool, Vec<bool>)> = (0..self.noperators)
            .map(|i| self.get_as_vec_bool(i))
            .collect();
        // We sort the operators by support size
        transposed.sort_by_key(|(_, vec)| {
            (0..self.n)
                .map(|i| if vec[i] | vec[i + self.n] { 1 } else { 0 })
                .sum::<i32>()
        });
        // We clear the set
        self.clear();
        // And insert everybody back in place
        for (phase, axis) in transposed {
            self.insert_vec_bool(&axis, phase);
        }
    }
}

impl fmt::Display for PauliSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.len() {
            let (phase, string) = self.get(i);
            writeln!(f, "{}{}", if phase { "-" } else { "+" }, string)?;
        }
        writeln!(f)
    }
}

impl PauliLike for PauliSet {
    /// Conjugate the set of rotations via a H gate
    fn h(&mut self, i: usize) {
        self.data_array.swap(i, i + self.n);
        self.update_phase_and(i, i + self.n);
    }
    /// Conjugate the set of rotations via a S gate
    fn s(&mut self, i: usize) {
        self.update_phase_and(i, i + self.n);
        self.row_op(i, i + self.n);
    }
    /// Conjugate the set of rotations via a S dagger gate
    fn sd(&mut self, i: usize) {
        self.row_op(i, i + self.n);
        self.update_phase_and(i, i + self.n);
    }
    /// Conjugate the set of rotations via a SQRT_X gate
    fn sqrt_x(&mut self, i: usize) {
        self.row_op(i + self.n, i);
        self.update_phase_and(i, i + self.n);
    }
    /// Conjugate the set of rotations via a SQRT_X dagger gate
    fn sqrt_xd(&mut self, i: usize) {
        self.update_phase_and(i, i + self.n);
        self.row_op(i + self.n, i);
    }
    /// Conjugate the set of rotations via a CNOT gate
    fn cnot(&mut self, i: usize, j: usize) {
        self.update_phase_and_many(i, j, i + self.n, j + self.n);
        self.row_op(j + self.n, i + self.n);
        self.row_op(i, j);
        self.update_phase_and_many(i, j, i + self.n, j + self.n);
    }
}

#[cfg(test)]
mod pauli_set_tests {
    use super::*;

    #[test]
    fn construction() {
        let pset = PauliSet::new(10);
        assert_eq!(pset.data_array.len(), 20);
        assert_eq!(pset.n, 10);
        assert_eq!(pset.nstrides, 0);
        assert_eq!(pset.noperators, 0);
    }

    #[test]
    fn insertion() {
        let mut pset = PauliSet::new(4);
        pset.insert("XYZI", false);
        assert_eq!(pset.data_array.len(), 8);
        assert_eq!(pset.n, 4);
        assert_eq!(pset.nstrides, 1);
        assert_eq!(pset.noperators, 1);
        assert_eq!(pset.data_array[0][0], 1);
        assert_eq!(pset.data_array[1][0], 1);
        assert_eq!(pset.data_array[2][0], 0);
        assert_eq!(pset.data_array[3][0], 0);

        assert_eq!(pset.data_array[4][0], 0);
        assert_eq!(pset.data_array[5][0], 1);
        assert_eq!(pset.data_array[6][0], 1);
        assert_eq!(pset.data_array[7][0], 0);
    }

    #[test]
    fn get() {
        let mut pset = PauliSet::new(4);
        pset.insert("XYZI", false);
        assert_eq!(pset.get(0), (false, "XYZI".to_owned()));
    }

    #[test]
    fn h_test() {
        let mut pset = PauliSet::new(1);
        pset.insert("X", false);
        pset.insert("Z", false);
        pset.insert("Y", false);
        pset.insert("I", false);

        pset.h(0);

        assert_eq!(pset.get(0), (false, "Z".to_owned()));
        assert_eq!(pset.get(1), (false, "X".to_owned()));
        assert_eq!(pset.get(2), (true, "Y".to_owned()));
        assert_eq!(pset.get(3), (false, "I".to_owned()));
    }

    #[test]
    fn s_test() {
        let mut pset = PauliSet::new(1);
        pset.insert("X", false);
        pset.insert("Z", false);
        pset.insert("Y", false);
        pset.insert("I", false);

        pset.s(0);

        assert_eq!(pset.get(0), (false, "Y".to_owned()));
        assert_eq!(pset.get(1), (false, "Z".to_owned()));
        assert_eq!(pset.get(2), (true, "X".to_owned()));
        assert_eq!(pset.get(3), (false, "I".to_owned()));
    }

    #[test]
    fn sqrt_x_test() {
        let mut pset = PauliSet::new(1);
        pset.insert("X", false);
        pset.insert("Z", false);
        pset.insert("Y", false);
        pset.insert("I", false);

        pset.sqrt_x(0);

        assert_eq!(pset.get(0), (false, "X".to_owned()));
        assert_eq!(pset.get(1), (true, "Y".to_owned()));
        assert_eq!(pset.get(2), (false, "Z".to_owned()));
        assert_eq!(pset.get(3), (false, "I".to_owned()));
    }
    // TODO:write this test :'(
    #[test]
    fn cnot_test() {
        const INPUTS: [&str; 16] = [
            "II", "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX",
            "ZY", "ZZ",
        ];
        const OUTPUTS: [(bool, &str); 16] = [
            (false, "II"),
            (false, "IX"),
            (false, "ZY"),
            (false, "ZZ"),
            (false, "XX"),
            (false, "XI"),
            (false, "YZ"),
            (true, "YY"),
            (false, "YX"),
            (false, "YI"),
            (true, "XZ"),
            (false, "XY"),
            (false, "ZI"),
            (false, "ZX"),
            (false, "IY"),
            (false, "IZ"),
        ];
        for (ins, (phase, outs)) in INPUTS.iter().zip(OUTPUTS.iter()) {
            let mut pset = PauliSet::new(2);
            pset.insert(ins, false);
            pset.cnot(0, 1);
            assert_eq!(pset.get(0), (*phase, (*outs).to_owned()));
        }
    }

    #[test]
    fn support_size_test() {
        let mut pset = PauliSet::new(4);
        pset.insert("XYIZ", false);
        pset.insert("XYII", false);
        pset.insert("IYIZ", false);
        pset.insert("IIII", false);
        assert_eq!(pset.support_size(0), 3);
        assert_eq!(pset.support_size(1), 2);
        assert_eq!(pset.support_size(2), 2);
        assert_eq!(pset.support_size(3), 0);
    }
    #[test]
    fn count_id() {
        let mut pset = PauliSet::new(5);
        pset.insert("IIIII", false);
        pset.insert("XIIII", false);
        pset.insert("XXIII", false);
        pset.insert("XXXII", false);
        pset.insert("XXXXI", false);
        for i in 0..5 {
            assert_eq!(pset.count_id(i), i + 1);
        }
    }
    #[test]
    fn sort_test() {
        let mut pset = PauliSet::new(4);
        pset.insert("IIII", false);
        pset.insert("XXII", false);
        pset.insert("XXXX", false);
        pset.insert("XIII", false);
        pset.insert("XXXI", false);
        pset.support_size_sort();
        assert_eq!(pset.get(0), (false, "IIII".to_owned()));
        assert_eq!(pset.get(1), (false, "XIII".to_owned()));
        assert_eq!(pset.get(2), (false, "XXII".to_owned()));
        assert_eq!(pset.get(3), (false, "XXXI".to_owned()));
        assert_eq!(pset.get(4), (false, "XXXX".to_owned()));
    }
    #[test]
    fn pop_test() {
        let mut pset = PauliSet::new(1);
        pset.insert("I", false);
        pset.insert("X", false);
        assert_eq!(pset.noperators, 2);
        pset.pop();
        assert_eq!(pset.noperators, 1);
        assert_eq!(pset.start_offset, 1);
        assert_eq!(pset.get(0), (false, "X".to_owned()));
    }
    #[test]
    fn commute_test() {
        let mut pset = PauliSet::new(2);
        pset.insert("ZI", false);
        pset.insert("XI", false);
        pset.insert("ZZ", false);
        pset.insert("XX", false);
        pset.insert("YY", false);
        assert!(pset.commute(0, 2));
        assert!(!pset.commute(0, 1));
        assert!(pset.commute(2, 3));
        assert!(pset.commute(2, 4));
        assert!(pset.commute(3, 4));
        assert!(pset.commute(1, 3));
    }
}
