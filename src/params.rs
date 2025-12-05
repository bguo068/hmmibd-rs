use crate::samples;
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
    #[error("invalid sample id. {0} is not in the pop file")]
    SampleIdNotInPopFile(String),
}

#[derive(Debug, Clone, Copy, Default)]
pub struct SimulationParams {
    pub id1: u32,
    pub id2: u32,
    pub r: f64,
    pub k: f64,
}

pub fn read_params_file(
    param_file_path: &str,
    samples: &samples::Samples,
) -> Result<Vec<SimulationParams>, Error> {
    let m = samples.m();
    use std::io::BufRead;
    let reader = std::fs::File::open(param_file_path).map(std::io::BufReader::new)?;
    let mut params_vec = Vec::<SimulationParams>::new();
    for line_res in reader.lines() {
        let mut params = SimulationParams::default();
        let mut num_fiels = 0;
        line_res?.split('\t').enumerate().try_for_each(
            |(which, field): (usize, &str)| -> Result<(), Error> {
                num_fiels += 1;
                match which {
                    0 => {
                        params.id1 = *m
                            .get(field)
                            .ok_or(Error::SampleIdNotInPopFile(field.to_string()))?
                    }
                    1 => {
                        params.id2 = *m
                            .get(field)
                            .ok_or(Error::SampleIdNotInPopFile(field.to_string()))?
                    }
                    2 => params.r = field.parse()?,
                    3 => params.k = field.parse()?,
                    _ => return Err(Error::TooManyFields),
                }
                Ok(())
            },
        )?;
        if num_fiels < 4 {
            return Err(Error::TooFewFields);
        }
        params_vec.push(params);
    }

    Ok(params_vec)
}

// #[test]
// pub fn test_read_params_file() -> Result<(), Error> {
//     let p = "tmp/params.txt";
//     let v = read_params_file(p)?;
//     dbg!(&v[..3]);
//     Ok(())
// }
