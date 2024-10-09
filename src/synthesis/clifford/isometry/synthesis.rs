use crate::structures::CliffordCircuit;
use crate::structures::IsometryTableau;
use crate::structures::Metric;

use super::count::isometry_count_synthesis;
use super::depth::isometry_depth_synthesis;

use super::common::fix_phases;

pub fn isometry_synthesis(
    isometry: &IsometryTableau,
    metric: &Metric,
    niter: usize,
) -> CliffordCircuit {
    let mut result = match metric {
        Metric::COUNT => isometry_count_synthesis(isometry, niter),
        Metric::DEPTH => isometry_depth_synthesis(isometry),
    };
    fix_phases(isometry, &mut result);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        routines::f2_linalg::Matrix,
        structures::{PauliLike, Tableau},
    };

    fn print_matrix(matrix: &Matrix) {
        for row in matrix.iter() {
            for elem in row.iter() {
                if *elem {
                    print!("1");
                } else {
                    print!("0");
                }
            }
            println!("");
        }
    }
    #[test]
    fn test_phases_clifford_count() {
        for _ in 0..10 {
            let n = 10;
            let tableau = Tableau::random(n).to_isometry();
            let circuit = isometry_synthesis(&tableau, &Metric::COUNT, 1);
            let mut simulated = IsometryTableau::new(n, 0);
            simulated.conjugate_with_circuit(&circuit);
            assert_eq!(simulated, tableau);
        }
    }

    #[test]
    fn test_phases_clifford_depth() {
        for _ in 0..10 {
            let n = 10;
            let tableau = Tableau::random(n).to_isometry();
            let circuit = isometry_synthesis(&tableau, &Metric::DEPTH, 1);
            let mut simulated = IsometryTableau::new(n, 0);
            simulated.conjugate_with_circuit(&circuit);
            assert_eq!(simulated, tableau);
        }
    }
}
