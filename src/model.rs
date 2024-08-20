use crate::{args::Arguments, data::InputData, matrix::*};

pub struct ModelParamState<'a> {
    pub model: ModelParameters<'a>,
    pub last_model: ModelParameters<'a>,
    pub finish_fit: bool,
    pub iiter: usize,
    pub num_info_site: usize,
    pub num_discord: usize,
}

impl<'a> ModelParamState<'a> {
    pub fn new(data: &'a InputData) -> Self {
        let model = ModelParameters::new(&data.args);
        let last_model = model;
        let finish_fit = false;
        Self {
            model,
            last_model,
            finish_fit,
            iiter: 0,
            num_info_site: 0,
            num_discord: 0,
        }
    }
    pub fn is_fit_finished(&self) -> bool {
        self.finish_fit
    }
}

#[derive(Clone, Copy)]
pub struct ModelParameters<'a> {
    /// logarithm of the probability of the final set of state assignments and
    /// the set of observations (calculated with the Viterbi algorithm).
    pub max_phi: f64,
    /// HMM initial state probability
    pub pi: [f64; 2],
    /// number of generations of recombination (1 of 2 free parameters in fit).
    pub k_rec: f64,
    /// cmd arguments
    pub args: &'a Arguments,
}

impl<'a> ModelParameters<'a> {
    pub fn new(args: &'a Arguments) -> Self {
        Self {
            args,
            max_phi: 0.0,
            pi: [0.5, 0.5],
            k_rec: 1.0,
        }
    }
}

pub struct PerChrModelVariables {
    /// Rabiner's emission probabily matrix: 1 x nsites
    /// This is not the full matrix, as it is calcualted per pair per site
    /// the observation symbol is already known
    // pub b: Matrix<f64>,
    /// Forward/alpha variable matrix: nstates (2) x nsites
    pub alpha: Matrix<f64>,
    /// Backward/alpha variable matrix: nstates (2) x nsites
    pub beta: Matrix<f64>,
    /// Rabiner's Delta variable matrix (best score):
    pub phi: Matrix<f64>,
    pub psi: Matrix<u8>,
    pub traj: Vec<u8>,
    pub scale: Vec<f64>,
    // non-missing sites per chromosome per sample pair
    pub nsites: usize,
}

impl Default for PerChrModelVariables {
    fn default() -> Self {
        Self::new()
    }
}
impl PerChrModelVariables {
    pub fn new() -> Self {
        Self {
            alpha: Matrix::from_shape(2, 0, 0.0),
            beta: Matrix::from_shape(2, 0, 0.0),
            phi: Matrix::from_shape(2, 0, 0.0),
            psi: Matrix::from_shape(2, 0, 0),
            traj: vec![],
            scale: vec![],
            nsites: 0,
        }
    }

    pub fn resize_and_clear(&mut self, data: &InputData, chrid: usize) {
        let r = data.sites.get_chrom_pos_idx_ranges(chrid);
        let n = r.1 - r.0;
        self.alpha.resize_and_clear(2, n, 0.0);
        self.beta.resize_and_clear(2, n + 1, 0.0); // TODO figure out plus 1
        self.phi.resize_and_clear(2, n, 0.0);
        self.psi.resize_and_clear(2, n, 0);
        self.traj.clear();
        self.traj.resize(n, 0);
        self.scale.clear();
        self.scale.resize(n, 0.0);
    }
}

#[derive(Clone, Copy, Default)]
pub struct PerSnpModelVariables {
    /// Rabiner's xi matrix but are restricted at a given sites. nstates (2) x nstates (x)
    pub xi: [[f64; 2]; 2],
    /// Rabiner's gamma matrix but are restricted a given sites: nsites (2)
    pub gamma: [f64; 2],
}
impl PerSnpModelVariables {
    pub fn new() -> Self {
        Self::default()
    }
}
