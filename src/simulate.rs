use super::matrix::Matrix;
use rand::seq::IndexedRandom;
use rand::{self, Rng};

#[derive(Debug, thiserror::Error)]
pub enum SimulationError {
    #[error("invalid num rows of frequency matrix")]
    FreqDimensionError,
    #[error("invalid length of genetic distance vector")]
    CDimensionError,
    #[error("invalid length of IBD state vector")]
    StatesDimensionError,
    #[error("invalid length of genotype vector")]
    GtDimensionError,
    // #[error("Invalid value of frequency matrix")]
    // FreqValueError,
    #[error("weights slice used in sampling is empty")]
    RngWeightError,
}

type Result<T> = std::result::Result<T, SimulationError>;

/// this function is adapted from paneljudge package
/// See: https://github.com/aimeertaylor/paneljudge/blob/b489f6ab4669d0d3c6ff6e3bdfa9666375889a31/R/simulate_Ys.R#L60C1-L60C12
pub fn simulate_genotype_for_pair_from_k_and_r(
    f_mat1: &Matrix<f64>,
    f_mat2: &Matrix<f64>,
    c_vec: &[f64],
    k: f64,
    r: f64,
    epsilon: f64,
    states: &mut [bool],
    gt1: &mut [u8],
    gt2: &mut [u8],
) -> Result<()> {
    let num_marker = f_mat1.get_nrows();
    let num_alleles_max = f_mat1.get_ncols();

    if num_marker == 0 {
        return Err(SimulationError::FreqDimensionError);
    }
    if c_vec.len() != num_marker - 1 {
        return Err(SimulationError::CDimensionError);
    }
    if states.len() != num_marker {
        return Err(SimulationError::StatesDimensionError);
    }
    if (gt1.len() != num_marker) || (gt2.len() != num_marker) {
        return Err(SimulationError::GtDimensionError);
    }
    let mut rng = rand::rng();

    let mut weights1: Vec<(usize, f64)> = Vec::with_capacity(num_alleles_max);
    let mut weights2: Vec<(usize, f64)> = Vec::with_capacity(num_alleles_max);

    for m in 0..num_marker {
        // transmission ====================================
        if m == 0 {
            // initialize
            states[0] = rng.random::<f64>() < r;
        } else {
            if states[m - 1] {
                states[m] =
                    rng.random::<f64>() < 1.0 - (1.0 - r) * (1.0 - (-k * c_vec[m - 1]).exp());
            } else {
                states[m] = rng.random::<f64>() < r * (1.0 - (-k * c_vec[m - 1]).exp());
            }
        }

        // emission ================================
        let f_vec1 = f_mat1.get_row_raw_slice(m);
        let f_vec2 = f_mat2.get_row_raw_slice(m);
        let gamma = f_vec1
            .iter()
            .zip(f_vec2.iter())
            .filter(|(&f1, &f2)| (f1 > 1e-10) || (f2 > 1e-10))
            .count();

        // populate weights for each site
        weights1.clear();
        f_vec1
            .iter()
            .enumerate()
            .filter(|(_idx, &f)| f > 1e-10)
            .for_each(|(idx, &f)| weights1.push((idx, f)));
        weights2.clear();
        f_vec2
            .iter()
            .enumerate()
            .filter(|(_idx, &f)| f > 1e-10)
            .for_each(|(idx, &f)| weights2.push((idx, f)));

        let true_gt1 = weights1
            .choose_weighted(&mut rng, |x| x.1)
            .map_err(|_| SimulationError::RngWeightError)?
            .0;

        let true_gt2 = if states[m] {
            // if ibd
            true_gt1
        } else {
            // if not ibd, independent draw
            weights2
                .choose_weighted(&mut rng, |x| x.1)
                .map_err(|_| SimulationError::RngWeightError)?
                .0
        };

        use rand::seq::IteratorRandom;
        // simulate genotping error
        gt1[m] = if rng.random::<f64>() < (gamma - 1) as f64 * epsilon {
            // uniform draw from alleles other than true gt1
            weights1
                .iter()
                .filter(|(i, _)| *i != true_gt1)
                .choose(&mut rng)
                .ok_or(SimulationError::RngWeightError)?
                .0 as u8
        } else {
            // not error
            true_gt1 as u8
        };

        gt2[m] = if rng.random::<f64>() < (gamma - 1) as f64 * epsilon {
            // uniform draw from alleles other than true gt2
            weights2
                .iter()
                .filter(|(i, _)| *i != true_gt2)
                .choose(&mut rng)
                .ok_or(SimulationError::RngWeightError)?
                .0 as u8
        } else {
            // not error
            true_gt2 as u8
        };
    }
    Ok(())
}

pub fn simulate_genotype_for_pair_from_states(
    f_mat1: &Matrix<f64>,
    f_mat2: &Matrix<f64>,
    c_vec: &[f64],
    epsilon: f64,
    states: &[bool],
    gt1: &mut [u8],
    gt2: &mut [u8],
) -> Result<()> {
    let num_marker = f_mat1.get_nrows();
    let num_alleles_max = f_mat1.get_ncols();

    if num_marker == 0 {
        return Err(SimulationError::FreqDimensionError);
    }
    if c_vec.len() != num_marker - 1 {
        return Err(SimulationError::CDimensionError);
    }
    if states.len() != num_marker {
        return Err(SimulationError::StatesDimensionError);
    }
    if (gt1.len() != num_marker) || (gt2.len() != num_marker) {
        return Err(SimulationError::GtDimensionError);
    }
    let mut rng = rand::rng();

    let mut weights1: Vec<(usize, f64)> = Vec::with_capacity(num_alleles_max);
    let mut weights2: Vec<(usize, f64)> = Vec::with_capacity(num_alleles_max);

    for m in 0..num_marker {
        // emission ================================
        let f_vec1 = f_mat1.get_row_raw_slice(m);
        let f_vec2 = f_mat2.get_row_raw_slice(m);
        let gamma = f_vec1
            .iter()
            .zip(f_vec2.iter())
            .filter(|(&f1, &f2)| (f1 > 1e-10) || (f2 > 1e-10))
            .count();

        // populate weights for each site
        weights1.clear();
        f_vec1
            .iter()
            .enumerate()
            .filter(|(_idx, &f)| f > 1e-10)
            .for_each(|(idx, &f)| weights1.push((idx, f)));
        weights2.clear();
        f_vec2
            .iter()
            .enumerate()
            .filter(|(_idx, &f)| f > 1e-10)
            .for_each(|(idx, &f)| weights2.push((idx, f)));

        let true_gt1 = weights1
            .choose_weighted(&mut rng, |x| x.1)
            .map_err(|_| SimulationError::RngWeightError)?
            .0;

        let true_gt2 = if states[m] {
            // if ibd
            true_gt1
        } else {
            // if not ibd, independent draw
            weights2
                .choose_weighted(&mut rng, |x| x.1)
                .map_err(|_| SimulationError::RngWeightError)?
                .0
        };

        use rand::seq::IteratorRandom;
        // simulate genotping error
        gt1[m] = if rng.random::<f64>() < (gamma - 1) as f64 * epsilon {
            // uniform draw from alleles other than true gt1
            weights1
                .iter()
                .filter(|(i, _)| *i != true_gt1)
                .choose(&mut rng)
                .ok_or(SimulationError::RngWeightError)?
                .0 as u8
        } else {
            // not error
            true_gt1 as u8
        };

        gt2[m] = if rng.random::<f64>() < (gamma - 1) as f64 * epsilon {
            // uniform draw from alleles other than true gt2
            weights2
                .iter()
                .filter(|(i, _)| *i != true_gt2)
                .choose(&mut rng)
                .ok_or(SimulationError::RngWeightError)?
                .0 as u8
        } else {
            // not error
            true_gt2 as u8
        };
    }
    Ok(())
}
