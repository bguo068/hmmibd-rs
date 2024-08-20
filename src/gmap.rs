pub struct GeneticMap(Vec<(u32, f64)>);
use crate::genome::GenomeInfo;
use csv;
use std::{
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("io error, source:{source:?}, path: {path:?}")]
    Io {
        source: std::io::Error,
        path: std::path::PathBuf,
    },
    #[error("error in converting Path to String: {0:?}")]
    PathToStringError(std::path::PathBuf),
    #[error("csv error: {0:?}")]
    CsvError(#[from] csv::Error),
    #[error("csv not enough column error: expect {expect} columns, found {actual} columns")]
    CsvNotEnoughColumns { expect: usize, actual: usize },
    #[error("{0:?}")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("{0:?}")]
    ParseIntError(#[from] std::num::ParseIntError),
}

impl FromIterator<(u32, f64)> for GeneticMap {
    fn from_iter<T: IntoIterator<Item = (u32, f64)>>(iter: T) -> Self {
        let v: Vec<_> = iter.into_iter().collect();
        Self(v)
    }
}

impl GeneticMap {
    // pub fn from_iter(it: impl Iterator<Item = (u32, f64)>) -> Self {
    //     let v = it.collect();
    //     Self(v)
    // }

    pub fn from_genome_info(ginfo: &GenomeInfo) -> Result<Self, Error> {
        let mut gw_chr_start_bp = 0u32;
        let mut gw_chr_start_cm = 0.0f64;
        let mut v = Vec::new();

        for (chrlen, plinkmap_fn) in ginfo.chromsize.iter().zip(ginfo.gmaps.iter()) {
            let mut chrmap = GeneticMap::from_plink_map(plinkmap_fn, *chrlen)?;
            let chrlen_cm = chrmap.get_size_cm();

            chrmap.update_to_genome_wide_coords(gw_chr_start_bp, gw_chr_start_cm);

            v.extend(chrmap.0);

            gw_chr_start_bp += chrlen;
            gw_chr_start_cm += chrlen_cm;
        }
        Ok(Self(v))
    }

    pub fn get_gw_chr_start_cm_from_chrid(&self, chrid: usize, ginfo: &GenomeInfo) -> f64 {
        let gw_ch_start_bp = ginfo.gwstarts[chrid];
        self.get_cm(gw_ch_start_bp)
    }

    pub fn get_gw_chr_start_cm_vec(&self, ginfo: &GenomeInfo) -> Vec<f64> {
        let mut v = Vec::new();
        for chrname in ginfo.chromnames.iter() {
            let chrid = ginfo.idx[chrname];
            let gw_ch_start_bp = ginfo.gwstarts[chrid];
            let gw_chr_start_cm = self.get_cm(gw_ch_start_bp);
            v.push(gw_chr_start_cm);
        }

        v
    }

    pub fn from_plink_map(p: impl AsRef<Path>, chrlen: u32) -> Result<Self, Error> {
        let mut v = vec![(0, 0.0)];
        let mut record = csv::StringRecord::new();

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b' ')
            .from_path(&p)?;

        while reader.read_record(&mut record)? {
            if record.len() < 3 {
                return Err(Error::CsvNotEnoughColumns {
                    expect: 3,
                    actual: record.len(),
                });
            }
            // println!("{:?}", record);
            let cm = record[2].parse::<f64>()?;
            // use 0-based position
            let bp: u32 = record[3].parse::<u32>()? - 1;
            // println!("record: cm={}, bp={}", cm, bp);
            if bp == 0 {
                continue;
            }
            assert!(
                (bp, cm) > v[v.len() - 1],
                "genetic map should be ordered by position"
            );
            v.push((bp, cm));
        }

        let (x1, y1) = v[v.len() - 2];
        let (x2, y2) = v[v.len() - 1];
        if x2 != chrlen {
            let slope = (y2 - y1) / ((x2 - x1) as f64);
            // println!("slope={slope}");
            let mut end_cm = ((chrlen - x2) as f64) * slope + y2;
            if end_cm < y2 {
                end_cm = y2;
            }
            v.push((chrlen, end_cm));
        }
        // println!("the third last pair: {:?}", v[v.len() - 3]);
        // println!("the second last pair: {:?}", v[v.len() - 2]);
        // println!("last pair: {:?}", v.last().unwrap());
        Ok(Self(v))
    }
    pub fn get_cm(&self, bp: u32) -> f64 {
        let idx = self.0.partition_point(|e| e.0 <= bp) - 1;
        let (x1, y1) = self.0[idx];
        let (x2, y2) = self.0[idx + 1];
        let slope = (y2 - y1) / (x2 - x1) as f64;
        let mut cm = (bp - x1) as f64 * slope + y1;
        if cm < y1 {
            cm = y1;
        } else if cm > y2 {
            cm = y2;
        }
        cm
    }

    pub fn get_cm_len(&self, s: u32, e: u32) -> f64 {
        self.get_cm(e) - self.get_cm(s)
    }

    pub fn get_bp(&self, cm: f64) -> u32 {
        let idx = self.0.partition_point(|e| e.1 <= cm) - 1;
        let (x1, y1) = self.0[idx];
        let (x2, y2) = self.0[idx + 1];
        let slope = (x2 - x1) as f64 / (y2 - y1);
        let mut bp = ((cm - y1) * slope) as u32 + x1;
        if bp < x1 {
            bp = x1;
        } else if bp > x2 {
            bp = x2;
        }
        bp
    }

    pub fn get_size_cm(&self) -> f64 {
        let n = self.0.len();
        self.0[n - 1].1 - self.0[0].1
    }

    pub fn update_to_genome_wide_coords(&mut self, gw_chr_start_bp: u32, gw_chr_start_cm: f64) {
        self.0.iter_mut().for_each(|(bp, cm)| {
            *bp += gw_chr_start_bp;
            *cm += gw_chr_start_cm;
        });
    }

    pub fn to_plink_map_files(
        &self,
        ginfo: &GenomeInfo,
        prefix: impl AsRef<Path>,
    ) -> Result<(), Error> {
        // make folder
        let parent = match prefix.as_ref().parent() {
            Some(parent) => parent,
            None => match prefix.as_ref() {
                path if path == Path::new("/") => Path::new("/"),
                _ => Path::new(""),
            },
        };
        if !parent.exists() {
            std::fs::create_dir_all(parent).map_err(|e| Error::Io {
                source: e,
                path: parent.to_path_buf(),
            })?;
        }
        let filename = prefix
            .as_ref()
            .file_name()
            .ok_or(Error::PathToStringError(prefix.as_ref().to_path_buf()))?
            .to_str()
            .ok_or(Error::PathToStringError(prefix.as_ref().to_path_buf()))?;
        for (i, chrname) in ginfo.chromnames.iter().enumerate() {
            // make a file name per chromosome
            let path = parent.with_file_name(format!("{}_{}.map", filename, chrname));
            let mut f = std::fs::File::create(&path)
                .map(BufWriter::new)
                .map_err(|e| Error::Io {
                    source: e,
                    path: path.to_path_buf(),
                })?;
            // find the first end record for each chromosome
            let gwstart = ginfo.gwstarts[i];
            let gwend = match ginfo.gwstarts.get(i + 1) {
                None => ginfo.get_total_len_bp(),
                Some(gwend) => *gwend,
            };
            let s = self.0.partition_point(|x| x.0 < gwstart);
            let e = self.0.partition_point(|x| x.0 < gwend);
            // write records for each chromosomes
            let (pos_offset, cm_offset) = self.0[s];
            for (pos, cm) in &self.0[s..e] {
                let pos = *pos - pos_offset + 1; // 1-based position
                let cm = *cm - cm_offset;

                writeln!(f, "{} . {} {}", chrname, cm, pos).map_err(|e| Error::Io {
                    source: e,
                    path: path.to_path_buf(),
                })?;
            }
        }
        Ok(())
    }
}
