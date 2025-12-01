#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0:?}")]
    Io(#[from] std::io::Error),
    #[error("Too many fields")]
    TooManyFields,
    #[error("Too few fields")]
    TooFewFields,
    #[error("{0:?}")]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error("{0:?}")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("Invalid pop id, should be 0 or 1")]
    InvalidPopId,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct SimulationParams {
    pub r: f64,
    pub k: f64,
    pub pop_id1: u8,
    pub pop_id2: u8,
}

pub fn read_params_file(param_file_path: &str) -> Result<(Vec<SimulationParams>, bool), Error> {
    use std::io::BufRead;
    let mut use_2nd_freq_file = false;
    let reader = std::fs::File::open(param_file_path).map(std::io::BufReader::new)?;
    let mut params_vec = Vec::<SimulationParams>::new();
    for line_res in reader.lines() {
        let mut params = SimulationParams::default();
        let mut num_fiels = 0;
        line_res?.split('\t').enumerate().try_for_each(
            |(which, field): (usize, &str)| -> Result<(), Error> {
                num_fiels += 1;
                match which {
                    0 => params.r = field.parse()?,
                    1 => params.k = field.parse()?,
                    2 => params.pop_id1 = field.parse()?,
                    3 => params.pop_id2 = field.parse()?,
                    _ => return Err(Error::TooManyFields),
                }
                Ok(())
            },
        )?;
        if num_fiels < 4 {
            return Err(Error::TooFewFields);
        }
        if (params.pop_id1 > 1) || (params.pop_id2 > 1) {
            return Err(Error::InvalidPopId);
        }
        if (!use_2nd_freq_file) && ((params.pop_id1 == 1) || (params.pop_id2 == 1)) {
            use_2nd_freq_file = true;
        }
        params_vec.push(params);
    }

    Ok((params_vec, use_2nd_freq_file))
}

// #[test]
// pub fn test_read_params_file() -> Result<(), Error> {
//     let p = "tmp/params.txt";
//     let v = read_params_file(p)?;
//     dbg!(&v[..3]);
//     Ok(())
// }
