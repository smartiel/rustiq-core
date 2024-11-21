use super::Pauli;

/// This trait should be implemented by any struct that can be measured by a Pauli operator
pub trait MeasurementLike {
    // Measure in the given basis
    fn measure(&mut self, basis: &Pauli);
}
