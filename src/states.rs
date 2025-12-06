use crate::{genome::Genome, samples, sites::Sites};
use slice_group_by::GroupBy;
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
    #[error("simulation genotype from states assumes one chromosome")]
    TooManyChromosome,
    #[error("Site should be sorted")]
    SiteNotSorted,
    #[error("When a pair is marked as non-IBD pair, there should one record for this pair")]
    TooManyStateRecordForAMarkerPair,
    #[error("segments in states file is overlapping for some pair")]
    OverlappingSegments,
    #[error("some segment has start position larger than or equal to end position ")]
    StartLargerThanEnd,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct StateRecord {
    pub id1: u32,
    pub id2: u32,
    // for pairs with not ibd segment sharing, should use id1, id2, 0, 0
    pub start: u32,
    pub end: u32,
}

pub struct States {
    states: Vec<StateRecord>,
    positions: Vec<u32>,
    state_pair: Vec<bool>,
    current_idx: usize,
}

impl States {
    pub fn from_state_file(
        state_file_path: &str,
        samples: &samples::Samples,
        sites: &Sites,
        genome: &Genome,
    ) -> Result<Self, Error> {
        let m = samples.m();
        use std::io::BufRead;
        let reader = std::fs::File::open(state_file_path).map(std::io::BufReader::new)?;
        let mut states = Vec::<StateRecord>::new();
        for line_res in reader.lines() {
            let mut params = StateRecord::default();
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
                        2 => params.start = field.parse()?,
                        3 => params.end = field.parse()?,
                        _ => return Err(Error::TooManyFields),
                    }
                    Ok(())
                },
            )?;
            if num_fiels < 4 {
                return Err(Error::TooFewFields);
            }
            if params.id1 > params.id2 {
                std::mem::swap(&mut params.id1, &mut params.id2);
            }
            states.push(params);
        }

        states.sort();
        let positions = sites.get_pos_slice().to_vec();
        if !positions
            .iter()
            .zip(positions.iter().skip(1))
            .all(|(a, b)| a.le(b))
        {
            return Err(Error::SiteNotSorted);
        }

        if genome.get_nchrom() > 1 {
            return Err(Error::TooManyChromosome);
        }
        let nposition = positions.len();

        let s = Self {
            states,
            positions,
            state_pair: vec![false; nposition],
            current_idx: 0,
        };
        s.check()?;

        Ok(s)
    }

    fn check(&self) -> Result<(), Error> {
        for blk in self.states.linear_group_by_key(|e| (e.id1, e.id2)) {
            if (blk.len() > 1) && blk.iter().any(|e| e.end == 0) {
                return Err(Error::TooManyStateRecordForAMarkerPair);
            }
            if (blk[0].end != 0) && (!blk.iter().all(|e| e.start < e.end)) {
                return Err(Error::StartLargerThanEnd);
            }
            if !blk
                .iter()
                .zip(blk.iter().skip(1))
                .all(|(a, b)| a.end <= b.start)
            {
                return Err(Error::OverlappingSegments);
            }
        }
        Ok(())
    }

    pub fn next_pair_and_states(&mut self) -> Option<((u32, u32), &[bool])> {
        let npos = self.positions.len();
        let nrec = self.states.len();
        if self.current_idx >= nrec {
            return None;
        }
        let blk = self.states[self.current_idx..]
            .linear_group_by_key(|e| (e.id1, e.id2))
            .next()
            .unwrap_or(&[]);
        self.current_idx += blk.len();
        // reset the vector
        self.state_pair.clear();
        self.state_pair.resize(npos, false);

        if blk[0].end != 0 {
            for e in blk {
                let i = self.positions.partition_point(|x| *x < e.start);
                let j = self.positions.partition_point(|x| *x <= e.end);
                self.state_pair[i..j].iter_mut().for_each(|s| {
                    *s = true;
                });
            }
        }
        Some(((blk[0].id1, blk[0].id2), &self.state_pair[..]))
    }

    pub fn get_num_unique_pairs(&self) -> u32 {
        self.states.linear_group_by_key(|e| (e.id1, e.id2)).count() as u32
    }
}
