use serde::Deserialize;
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("IO error, source: {source:?}, path: {path:?}")]
    Io {
        source: std::io::Error,
        path: std::path::PathBuf,
    },
    #[error("toml parsing error: {0:?}")]
    TomlParsingError(#[from] toml::de::Error),

    #[error("GenomeError: {0}")]
    GenomeValidationError(#[from] GenomeValidationError),
}

#[derive(Debug, thiserror::Error)]
pub enum GenomeValidationError {
    #[error("chromsize vector and chromname vector have different length")]
    InconsistencyChromsizeVsChromnames,
    #[error("chromsize vector and gmap vector have different length")]
    InconsistencyChromsizeVsGmaps,
    #[error("not all chromnames are in idx hashmap")]
    ChromnamesNotInIdx,
    #[error("Gwstart vector is empty")]
    GwstartsIsEmpty,
}

#[derive(Deserialize, Debug)]
struct GenomeFile {
    name: String,
    chromsize: Vec<u32>,
    idx: HashMap<String, usize>,
    chromnames: Vec<String>,
    gmaps: Vec<String>,
}

impl GenomeFile {
    fn check(&self) -> Result<()> {
        if self.chromsize.len() != self.chromnames.len() {
            Err(GenomeValidationError::InconsistencyChromsizeVsChromnames.into())
        } else if self.chromsize.len() != self.gmaps.len() {
            Err(GenomeValidationError::InconsistencyChromsizeVsGmaps.into())
        } else if !self.chromnames.iter().all(|x| self.idx.contains_key(x)) {
            Err(GenomeValidationError::ChromnamesNotInIdx.into())
        } else {
            Ok(())
        }
    }
}

#[test]
fn test_genome() {
    let mut s = String::new();
    std::fs::File::open("sim_data/genome.toml")
        .unwrap()
        .read_to_string(&mut s)
        .unwrap();
    let _: GenomeFile = toml::from_str(&s).unwrap();
}

#[derive(Default)]
pub struct GenomeInfo {
    pub name: String,
    pub chromsize: Vec<u32>,
    pub chromnames: Vec<String>,
    pub idx: HashMap<String, usize>,
    pub gwstarts: Vec<u32>,
    pub gmaps: Vec<String>,
}

impl GenomeInfo {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_toml_file<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let mut s = String::new();
        let p: &Path = path.as_ref();
        std::fs::File::open(p)
            .map_err(|e| Error::Io {
                source: e,
                path: p.to_owned(),
            })?
            .read_to_string(&mut s)
            .map_err(|e| Error::Io {
                source: e,
                path: p.to_owned(),
            })?;
        let mut gfile: GenomeFile = toml::from_str(&s)?;
        gfile.check()?;
        let mut ginfo = Self::new();
        use std::mem::swap;
        swap(&mut gfile.name, &mut ginfo.name);
        swap(&mut gfile.chromsize, &mut ginfo.chromsize);
        swap(&mut gfile.chromnames, &mut ginfo.chromnames);
        swap(&mut gfile.idx, &mut ginfo.idx);
        swap(&mut gfile.gmaps, &mut ginfo.gmaps);
        ginfo.gwstarts.push(0u32);
        for chrlen in &ginfo.chromsize[0..(ginfo.chromsize.len() - 1)] {
            let last = ginfo
                .gwstarts
                .last()
                .ok_or(GenomeValidationError::GwstartsIsEmpty)?;
            ginfo.gwstarts.push(*chrlen + *last);
        }
        Ok(ginfo)
    }
    pub fn to_gw_pos(&self, chrid: usize, pos: u32) -> u32 {
        self.gwstarts[chrid] + pos
    }
    pub fn to_chr_pos(&self, gw_pos: u32) -> (usize, &str, u32) {
        let chrid = self.gwstarts.partition_point(|x| *x <= (gw_pos)) - 1;
        let pos = gw_pos - self.gwstarts[chrid];
        let chrname = &self.chromnames[chrid];
        (chrid, chrname, pos)
    }
    pub fn get_total_len_bp(&self) -> u32 {
        self.chromsize.iter().sum()
    }

    // pub fn split_chromosomes_by_regions(&self, regions: &Intervals<u32>) -> Self {
    //     let name = format!("{}_rmpeaks", &self.name);
    //     let mut chromsize = vec![];
    //     let mut idx = HashMap::new();
    //     let mut gwstarts = vec![];
    //     let mut chromnames = vec![];
    //     let gmaps = vec![];
    //     let gw_tot_len_bp = self.get_total_len_bp();

    //     // new gwstart are from old gwstarts plus region boundaries
    //     gwstarts.extend(self.gwstarts.iter());
    //     regions
    //         .iter()
    //         .for_each(|r| gwstarts.extend(&[r.start, r.end]));
    //     // dedup
    //     gwstarts.sort();
    //     gwstarts.dedup();
    //     // filter boundaries beyond the total size
    //     gwstarts.retain(|x| *x < gw_tot_len_bp);

    //     for (i, gwstart) in gwstarts.iter().enumerate() {
    //         let (_chrid, chrname, chrpos) = self.to_chr_pos(*gwstart);
    //         let new_chrname = format!("{}_{}", chrname, chrpos);
    //         chromnames.push(new_chrname.clone());
    //         idx.insert(new_chrname, i);
    //     }

    //     // use gwstarts to get chrommsome sizes
    //     gwstarts.push(gw_tot_len_bp);
    //     chromsize.extend(
    //         gwstarts
    //             .iter()
    //             .zip(gwstarts.iter().skip(1))
    //             .map(|(a, b)| *b - *a),
    //     );
    //     gwstarts.pop();

    //     Self {
    //         name,
    //         chromnames,
    //         chromsize,
    //         idx,
    //         gwstarts,
    //         gmaps,
    //     }
    // }
}

#[derive(Clone)]
pub struct Genome {
    chromnames: Vec<String>,
    // chromsizes: Vec<u32>,
    gwchrstarts: Vec<u32>,
    idx: HashMap<String, u32>,
}

impl Genome {
    pub fn to_gw_pos(&self, chrname: &str, chr_pos: u32) -> u32 {
        self.gwchrstarts[self.idx[chrname] as usize] + chr_pos
    }
    pub fn to_chr_pos(&self, gw_pos: u32) -> (u32, &str, u32) {
        // x[0], x[1], x[2], x[3]
        //     x[i] <= u  < x[i+1]
        let chrid = self.gwchrstarts.partition_point(|x| *x <= gw_pos) - 1;
        let chrname = &self.chromnames[chrid];
        let chr_pos = gw_pos - self.gwchrstarts[chrid];
        (chrid as u32, chrname, chr_pos)
    }
    pub fn get_gwchrstarts(&self) -> &[u32] {
        &self.gwchrstarts[..]
    }

    pub fn is_identical(&self, other: &Self) -> bool {
        (self.chromnames == other.chromnames) && (self.gwchrstarts == other.gwchrstarts)
    }
    pub fn get_nchrom(&self) -> u32 {
        self.get_gwchrstarts().len() as u32
    }
    pub fn get_chrname(&self, idx: usize) -> &str {
        &self.chromnames[idx]
    }
    pub fn from_genome_info(ginfo: &GenomeInfo) -> Self {
        let chromnames = ginfo.chromnames.clone();
        let gwchrstarts = ginfo.gwstarts.clone();
        let idx = ginfo
            .idx
            .iter()
            .map(|(k, v)| (k.clone(), *v as u32))
            .collect();
        Self {
            chromnames,
            gwchrstarts,
            idx,
        }
    }
}

/// It is simple gnome builder, useful when no genetic map and genome infor is provided
#[derive(Default)]
pub struct GenomeBuilder {
    // chromsome and its index
    gwchrstarts: Vec<u32>,
    chromnames: Vec<String>,
    idx: HashMap<String, (u32, u32)>,
}

impl GenomeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// convert string to idx and inform GENOME about the chromosme length
    pub fn encode_pos(&mut self, chrname: &str, pos: u32) -> Result<(u32, u32)> {
        if self.idx.contains_key(chrname) {
            let (chrid, max_pos) = self.idx[chrname];
            assert!(max_pos < pos);
            *self
                .idx
                .get_mut(chrname)
                .ok_or(GenomeValidationError::ChromnamesNotInIdx)? = (chrid, pos);
            let gw_chr_start = self
                .gwchrstarts
                .last()
                .ok_or(GenomeValidationError::GwstartsIsEmpty)?;
            Ok((chrid, pos + gw_chr_start))
        } else {
            let last_max_pos = match self.chromnames.last() {
                None => 0,
                Some(chrname) => self.idx[chrname].1 + 1,
            };
            let last_gw_start = match self.gwchrstarts.last() {
                Some(x) => *x,
                None => 0,
            };
            // println!("last+max+pos-{last_max_pos} lat_gw_start={last_gw_start}");
            let this_gw_start = last_max_pos + last_gw_start;
            self.gwchrstarts.push(this_gw_start);

            let chrid = self.idx.len() as u32;
            self.chromnames.push(chrname.to_owned());
            self.idx.insert(chrname.to_owned(), (chrid, pos));

            Ok((chrid, pos + this_gw_start))
        }
    }

    pub fn finish(&mut self) -> Genome {
        let idx: HashMap<_, _> = std::mem::take(&mut self.idx)
            .into_iter()
            .map(|(chrname, (chrid, _max_pos))| (chrname, chrid))
            .collect();
        let chromnames = std::mem::take(&mut self.chromnames);
        let gwchrstarts = std::mem::take(&mut self.gwchrstarts);

        Genome {
            gwchrstarts,
            idx,
            chromnames,
        }
    }
}
