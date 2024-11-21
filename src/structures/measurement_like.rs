use super::Pauli;

/// This trait should be implemented by any struct that can be measured by a Pauli operator
pub trait MeasurementLike {
    fn measure(&mut self, p: Pauli);
}
