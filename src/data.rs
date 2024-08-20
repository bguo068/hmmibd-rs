use itertools::{EitherOrBoth, Itertools};
use smallvec::SmallVec;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read},
    path::Path,
    sync::{Arc, Mutex},
};

use crate::{
    args::Arguments,
    bcf::{self, DominantGenotype, DominantGenotypeArgs},
    genome::Genome,
    matrix::*,
    samples::{self, Samples},
    sites::{self, SiteInfoRaw, Sites},
};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("error in processing bcf file: {0:#?}")]
    Bcf(#[from] crate::bcf::Error),
    #[error("io error, source {source:?}, file: {file:?}")]
    Io {
        source: std::io::Error,
        file: Option<String>,
    },
    // #[error("custom: {0:?}")]
    // Custom(String),
    #[error("{0:?}")]
    ParseLineError(#[from] ParseLineError),
    #[error("irow={irow}, iallele={iallele}, value={value:?}")]
    MissingFrequency {
        irow: usize,
        iallele: usize,
        value: Option<f64>,
    },
    #[error("EmptyIterator")]
    EmptyIterator,

    #[error("lockerror: {0}")]
    LockError(&'static str),

    #[error("sites error: {0:?}")]
    Site(#[from] sites::Error),

    #[error("sample error: {0:?}")]
    Sample(#[from] samples::Error),

    #[error("sample error: {0:?}")]
    PartialCmpIsNone(&'static str),

    #[error("toml serialize error: {0:?}")]
    TomlSerializeError(#[from] toml::ser::Error),
}

#[derive(Debug, thiserror::Error)]
pub enum ParseLineError {
    #[error("Cannot read column {0}")]
    ReadColumnError(&'static str),
    #[error("Cannot parse column {0}")]
    ParseColumnError(&'static str),
    #[error("Cannot read line")]
    ReadLineError,
}

pub struct InputData {
    /// arguments
    pub args: Arguments,
    /// genotype matrix for population 1
    pub geno: Matrix<u8>,
    /// an array of number of unique alleles, length number sites
    pub nall: Vec<u8>,
    pub majall: Vec<u8>,
    /// genotype matrix for population 2
    pub freq1: Matrix<f64>,
    /// allele frequency matrix for population 2
    pub freq2: Option<Matrix<f64>>,
    pub sites: Sites,
    pub genome: Genome,
    pub samples: Samples,
    pub pairs: Vec<(u32, u32)>,
}

impl InputData {
    /// The first element of return value is always a chunk index in the samples in pop1
    /// The 2nd element of return value is a chunk index that can be in pop1 or pop2.
    /// the chunk index for pop2 is always no less than max chunk index of pop1 which can be used
    /// to determine which population the chunk is from
    pub fn get_chunk_pairs(&self) -> Vec<(u32, u32)> {
        let mut chunk_pairs = vec![];
        let mut num_chunk1 = self.samples.pop1_nsam() / self.args.par_chunk_size;
        if self.samples.pop1_nsam() % self.args.par_chunk_size > 0 {
            num_chunk1 += 1;
        }
        if self.freq2.is_some() {
            // two different populations
            assert!(self.samples.pop2_nsam() > 0);
            let mut num_chunk2 = self.samples.pop2_nsam() / self.args.par_chunk_size;
            if self.samples.pop2_nsam() % self.args.par_chunk_size > 0 {
                num_chunk2 += 1;
            }
            for i in 0..num_chunk1 {
                for j in 0..num_chunk2 {
                    let j = num_chunk1 + j;
                    chunk_pairs.push((i, j));
                }
            }
        } else {
            // same population
            for i in 0..num_chunk1 {
                for j in i..num_chunk1 {
                    chunk_pairs.push((i, j));
                }
            }
        }
        chunk_pairs
    }

    pub fn clone_inputdata_for_chunkpair(&self, chunkpair: (u32, u32)) -> Result<Self, Error> {
        let (ichunk, jchunk) = chunkpair;

        // num of chunks for samples.s1()
        let num_chunk1 = {
            let a = self.samples.pop1_nsam();
            let b = self.args.par_chunk_size;
            let mut c = a / b;
            let d = a % b;
            if d > 0 {
                c += 1;
            }
            c
        };

        // chunk1 start,end rows/samples, len
        let (ichunk_start, ichunk_end, ichunk_len) = {
            let ichunk_start = ichunk * self.args.par_chunk_size;
            let mut ichunk_end = (1 + ichunk) * self.args.par_chunk_size;
            if ichunk_end > self.samples.pop1_nsam() {
                ichunk_end = self.samples.pop1_nsam();
            }
            let ichunk_len = ichunk_end - ichunk_start;
            (ichunk_start, ichunk_end, ichunk_len)
        };

        // chunk2  start, end rows/samples, len
        let (jchunk_start, jchunk_end, jchunk_len) = {
            if jchunk < num_chunk1 {
                let jchunk_start = jchunk * self.args.par_chunk_size;
                let mut jchunk_end = (1 + jchunk) * self.args.par_chunk_size;
                if jchunk_end > self.samples.pop1_nsam() {
                    jchunk_end = self.samples.pop1_nsam();
                }
                let jchunk_len = jchunk_end - jchunk_start;
                (jchunk_start, jchunk_end, jchunk_len)
            } else {
                assert!(self.samples.pop2_nsam() > 0);
                let jjchunk = jchunk - num_chunk1;
                let jchunk_start = self.samples.pop1_nsam() + jjchunk * self.args.par_chunk_size;
                let mut jchunk_end =
                    self.samples.pop1_nsam() + (1 + jjchunk) * self.args.par_chunk_size;
                if jchunk_end > self.samples.pop1_nsam() + self.samples.pop2_nsam() {
                    jchunk_end = self.samples.pop1_nsam() + self.samples.pop2_nsam();
                }
                (jchunk_start, jchunk_end, jchunk_end - jchunk_start)
            }
        };
        let nsites = self.geno.get_ncols();

        let args = self.args.clone();
        let nall = self.nall.clone();
        let majall = self.majall.clone();
        let freq1 = self.freq1.clone();
        let geno = {
            let mut nrows = ichunk_len as usize;
            let ncols = nsites;
            let mut s = ichunk_start as usize * nsites;
            let mut e = ichunk_end as usize * nsites;
            let mut data = self.geno.as_slice()[s..e].to_vec();
            if ichunk != jchunk {
                s = jchunk_start as usize * nsites;
                e = jchunk_end as usize * nsites;
                data.extend_from_slice(&self.geno.as_slice()[s..e]);
                nrows += jchunk_len as usize;
            }
            Matrix::<u8>::from_shape_vec(nrows, ncols, data)
        };
        let freq2 = match (ichunk == jchunk, jchunk < num_chunk1) {
            (true, _) => None,
            (false, true) => Some(self.freq1.clone()),
            (false, false) => self.freq2.clone(),
        };
        let samples = {
            let pop1_slice = &self.samples.v()[(ichunk_start as usize)..(ichunk_end as usize)];
            if ichunk == jchunk {
                Samples::from_slice(pop1_slice, &[])
            } else {
                let pop2_slice = &self.samples.v()[(jchunk_start as usize)..(jchunk_end as usize)];
                Samples::from_slice(pop1_slice, pop2_slice)
            }
        };
        let sites = self.sites.clone();
        let genome = self.genome.clone();
        let pairs = Self::get_valid_pair_file(&args, &samples)?;
        Ok(Self {
            args,
            geno,
            nall,
            majall,
            freq1,
            freq2,
            sites,
            genome,
            samples,
            pairs,
        })
    }

    pub fn from_args(args: &Arguments) -> Result<Self, Error> {
        let dgt = if args.from_bcf {
            let dga = match args.dom_gt_config.as_ref() {
                Some(dom_gt_config_path) => {
                    DominantGenotypeArgs::new_from_toml_file(dom_gt_config_path)
                }
                None => {
                    let config = DominantGenotypeArgs::new_from_builtin()?;
                    std::fs::write("tmp_dom_gt_config.toml", toml::to_string(&config)?).map_err(
                        |e| Error::Io {
                            source: e,
                            file: args.dom_gt_config.to_owned(),
                        },
                    )?;
                    eprintln!(concat!(
                        "WARN: --dom-gt-config not specified, a builtin configuration is used",
                        " and is written to 'tmp_dom_gt_config.toml'",
                    ));
                    Ok(config)
                }
            }?;
            let dgt = DominantGenotype::new_from_processing_bcf(&dga, &args.data_file1)?;
            Some(dgt)
        } else {
            None
        };

        let valid_samples = Samples::from_args(args, dgt.as_ref())?;
        let min_snp_sep = args.min_snp_sep;

        let (mut geno1, geno2, sitesinfo) = match args.from_bcf {
            true => {
                let dgt = dgt.ok_or(bcf::Error::DominantGenotypeEmpty {
                    file: file!(),
                    line: line!(),
                })?;
                let (geno1, sitesinfo) =
                    Self::read_data_dominant_genotype(dgt, &valid_samples, min_snp_sep);
                (geno1, None, sitesinfo)
            }
            false => {
                let (geno1, sitesinfo) =
                    Self::read_data_hmmibd_format(&args.data_file1, &valid_samples, min_snp_sep)?;
                let geno2 = match args.data_file2.as_ref() {
                    Some(data_file2) => {
                        let (geno2, sites2) =
                            Self::read_data_hmmibd_format(data_file2, &valid_samples, min_snp_sep)?;
                        assert_eq!(&sitesinfo, &sites2);
                        Some(geno2)
                    }
                    None => None,
                };
                (geno1, geno2, sitesinfo)
            }
        };

        // create the genome object and site object
        let (sites, genome) = sitesinfo.into_sites_and_genome(&args.rec_args)?;

        let nsam1_valid = geno1.get_ncols();
        let nsam2_valid = geno2.as_ref().map(|geno2| geno2.get_ncols());
        let freq1 = if let Some(freq_file1) = args.freq_file1.as_ref() {
            Self::read_freq_file(freq_file1, args, &genome, &sites)?
        } else {
            Self::infer_freq_from_data(&geno1, &sites, args)
        };

        let freq2 = match (geno2.as_ref(), args.freq_file2.as_ref()) {
            (Some(_), Some(freq_file2)) => {
                Some(Self::read_freq_file(freq_file2, args, &genome, &sites)?)
            }
            (Some(geno2), None) => Some(Self::infer_freq_from_data(geno2, &sites, args)),
            _ => None,
        };
        let pairs = Self::get_valid_pair_file(args, &valid_samples)?;
        let nall = Self::get_nall(&freq1, freq2.as_ref())?;
        let majall = Self::get_major_all(&freq1, freq2.as_ref(), nsam1_valid, nsam2_valid)?;

        // sample oriented
        geno1.transpose();
        if let Some(mut geno2) = geno2 {
            geno2.transpose();
            geno1.merge(&geno2);
        }

        Ok(Self {
            args: args.clone(),
            geno: geno1,
            nall,
            majall,
            freq1,
            freq2,
            samples: valid_samples,
            sites,
            genome,
            pairs,
        })
    }
    // fn read_data_file_with_ginfo_and_gmap(
    //     data_file: impl AsRef<Path>,
    //     valid_samples: &Samples,
    //     min_snp_sep: u32,
    //     ginfo: &GenomeInfo,
    //     gmap: &GeneticMap,
    //     genome: &Genome,
    // ) -> (Matrix<u8>, Sites) {
    //     let mut sites = Sites::new();

    //     let mut line = String::with_capacity(100000);
    //     let mut f = std::fs::File::open(data_file.as_ref())
    //         .map(BufReader::new)
    //         .unwrap();

    //     // get all sample names
    //     f.read_line(&mut line).unwrap();
    //     let samples: Vec<_> = line
    //         .trim()
    //         .split('\t')
    //         .skip(2)
    //         .map(|x| x.to_owned())
    //         .collect();
    //     line.clear();

    //     let n_valid_samples = samples
    //         .iter()
    //         .filter(|s| valid_samples.m().contains_key(*s))
    //         .count();

    //     let mut geno1 = MatrixBuilder::<u8>::new(n_valid_samples);
    //     let mut last_chrname = String::new();
    //     let mut last_pos = 0;
    //     while f.read_line(&mut line).unwrap() != 0 {
    //         let mut fields = line.trim().split("\t");
    //         let chrname = fields.next().unwrap();
    //         let pos: u32 = fields.next().unwrap().parse().unwrap();

    //         // println!("pos: {pos}, last_pos: {last_pos}");
    //         if (chrname == last_chrname) && (last_pos + min_snp_sep > pos) {
    //             line.clear();
    //             continue;
    //         } else {
    //             last_chrname.clear();
    //             last_chrname.push_str(chrname);
    //             last_pos = pos;
    //         }

    //         let chrid = ginfo.idx[chrname];
    //         let gw_pos = ginfo.to_gw_pos(chrid, pos);
    //         let gw_pos_cm = gmap.get_cm(gw_pos);

    //         sites.add(gw_pos, gw_pos_cm);

    //         let mut cnt = 0;
    //         fields.enumerate().for_each(|(i, field)| {
    //             if valid_samples.m().contains_key(&samples[i]) {
    //                 let allel = match field.parse::<i8>().unwrap() {
    //                     -1 => None,
    //                     x => Some(x as u8),
    //                 };
    //                 geno1.push(allel);
    //                 cnt += 1;
    //             }
    //         });
    //         assert_eq!(cnt, n_valid_samples);
    //         line.clear();
    //     }

    //     let geno1 = geno1.finish();
    //     sites.finish(&genome);

    //     (geno1, sites)
    // }

    fn read_data_dominant_genotype(
        dgt: DominantGenotype,
        valid_samples: &Samples,
        min_snp_sep: u32,
    ) -> (Matrix<u8>, SiteInfoRaw) {
        // check seleted_samples are correct
        dgt.into_genotype_siteinfo(valid_samples, min_snp_sep)
    }

    fn read_data_hmmibd_format(
        data_file: impl AsRef<Path>,
        valid_samples: &Samples,
        min_snp_sep: u32,
    ) -> Result<(Matrix<u8>, SiteInfoRaw), Error> {
        let mut siteinfo = SiteInfoRaw::new();

        let mut line = String::with_capacity(100000);
        let mut f = std::fs::File::open(data_file.as_ref())
            .map(BufReader::new)
            .map_err(|source| Error::Io {
                source,
                file: Some(data_file.as_ref().to_string_lossy().to_string()),
            })?;

        // get all sample names
        f.read_line(&mut line).map_err(|e| Error::Io {
            source: e,
            file: None,
        })?;
        let samples: Vec<_> = line
            .trim()
            .split('\t')
            .skip(2)
            .map(|x| x.to_owned())
            .collect();
        line.clear();

        let valid_sample_col: Vec<usize> = samples
            .iter()
            .enumerate()
            .filter(|(_i, s)| valid_samples.m().contains_key(*s))
            .map(|(i, _)| i)
            .collect();
        let n_valid_samples = valid_sample_col.len();

        let mut geno = MatrixBuilder::<u8>::new(n_valid_samples);
        let mut last_chrname = String::new();
        let mut last_pos = 0;

        while f
            .read_line(&mut line)
            .map_err(|source| Error::Io { source, file: None })?
            != 0
        {
            let mut fields = line.trim().split("\t");
            let chrname = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("chrname"))?;
            let pos: u32 = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("pos"))?
                .parse()
                .map_err(|_| ParseLineError::ParseColumnError("pos"))?;

            // println!("pos: {pos}, last_pos: {last_pos}");
            if chrname == last_chrname {
                if last_pos + min_snp_sep > pos {
                    line.clear();
                    continue;
                } else {
                    last_chrname.clear();
                    last_chrname.push_str(chrname);
                    last_pos = pos;
                }
            } else {
                last_chrname.clear();
                last_chrname.push_str(chrname);
                last_pos = pos;
                siteinfo.add_chr_name(chrname);
            }

            siteinfo.add_chr_idx(chrname);
            siteinfo.add_chr_pos(pos);

            let mut cnt = 0;

            // use merge_join_by to only parse columns that contain selected samples
            fields
                .enumerate()
                .merge_join_by(valid_sample_col.iter(), |a, b| a.0.cmp(b))
                .try_for_each(|mergeby_res| -> Result<(), Error> {
                    if let EitherOrBoth::Both((_, field), _) = mergeby_res {
                        let allel = match field
                            .parse::<i8>()
                            .map_err(|_| ParseLineError::ParseColumnError("alleles"))?
                        {
                            -1 => None,
                            x => Some(x as u8),
                        };
                        geno.push(allel);
                        cnt += 1;
                    }
                    Ok(())
                })?;
            assert_eq!(cnt, n_valid_samples);
            line.clear();
        }

        let geno = geno.finish();
        Ok((geno, siteinfo))
    }

    // fn read_data_file(
    //     data_file: impl AsRef<Path>,
    //     valid_samples: &Samples,
    //     min_snp_sep: u32,
    //     rec_rate: f64,
    // ) -> (Genome, Matrix<u8>, Sites) {
    //     let mut sites = Sites::new();
    //     let mut gbuilder = GenomeBuilder::new();

    //     let mut line = String::with_capacity(100000);
    //     let mut f = std::fs::File::open(data_file.as_ref())
    //         .map(BufReader::new)
    //         .unwrap();

    //     // get all sample names
    //     f.read_line(&mut line).unwrap();
    //     let samples: Vec<_> = line
    //         .trim()
    //         .split('\t')
    //         .skip(2)
    //         .map(|x| x.to_owned())
    //         .collect();
    //     line.clear();

    //     let n_valid_samples = samples
    //         .iter()
    //         .filter(|s| valid_samples.m().contains_key(*s))
    //         .count();

    //     let mut geno1 = MatrixBuilder::<u8>::new(n_valid_samples);
    //     let mut last_chrname = String::new();
    //     let mut last_pos = 0;
    //     while f.read_line(&mut line).unwrap() != 0 {
    //         let mut fields = line.trim().split("\t");
    //         let chrname = fields.next().unwrap();
    //         let pos: u32 = fields.next().unwrap().parse().unwrap();

    //         // println!("pos: {pos}, last_pos: {last_pos}");
    //         if (chrname == last_chrname) && (last_pos + min_snp_sep > pos) {
    //             line.clear();
    //             continue;
    //         } else {
    //             last_chrname.clear();
    //             last_chrname.push_str(chrname);
    //             last_pos = pos;
    //         }

    //         let (_, gw_pos) = gbuilder.encode_pos(chrname, pos);

    //         sites.add(gw_pos, gw_pos as f64 * rec_rate * 100.0);

    //         let mut cnt = 0;
    //         fields.enumerate().for_each(|(i, field)| {
    //             if valid_samples.m().contains_key(&samples[i]) {
    //                 let allel = match field.parse::<i8>().unwrap() {
    //                     -1 => None,
    //                     x => Some(x as u8),
    //                 };
    //                 geno1.push(allel);
    //                 cnt += 1;
    //             }
    //         });
    //         assert_eq!(cnt, n_valid_samples);
    //         line.clear();
    //     }

    //     let genome = gbuilder.finish();
    //     let geno1 = geno1.finish();
    //     sites.finish(&genome);

    //     (genome, geno1, sites)
    // }
    fn read_freq_file(
        freq_file: impl AsRef<Path>,
        args: &Arguments,
        genome: &Genome,
        sites: &Sites,
    ) -> Result<Matrix<f64>, Error> {
        let mut f = std::fs::File::open(freq_file.as_ref())
            .map(BufReader::new)
            .map_err(|source| Error::Io {
                source,
                file: Some(freq_file.as_ref().to_string_lossy().to_string()),
            })?;

        let mut freq = MatrixBuilder::<f64>::new(args.max_all as usize);
        let mut v = Vec::new();

        let mut line = String::with_capacity(100000);
        line.clear();

        let mut last_chrname = String::new();
        let mut last_pos = 0;
        while f
            .read_line(&mut line)
            .map_err(|_| ParseLineError::ReadLineError)?
            != 0
        {
            let mut fields = line.trim().split("\t");
            let chrname = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("chrname"))?;
            let pos: u32 = fields
                .next()
                .ok_or(ParseLineError::ReadColumnError("pos"))?
                .parse()
                .map_err(|_| ParseLineError::ParseColumnError("pos"))?;
            if (chrname == last_chrname) && (last_pos + args.min_snp_sep > pos) {
                line.clear();
                continue;
            } else {
                last_chrname.clear();
                last_chrname.push_str(chrname);
                last_pos = pos;
            }

            let gw_pos = genome.to_gw_pos(chrname, pos);
            assert!(sites.has_gw_pos(gw_pos));

            v.clear();
            v.resize(args.max_all as usize, 0.0);
            fields
                .enumerate()
                .try_for_each(|(i, field)| -> Result<(), Error> {
                    let af = field
                        .parse::<f64>()
                        .map_err(|_| ParseLineError::ParseColumnError("freq"))?;
                    v[i] = af;
                    Ok(())
                })?;
            for af in v.iter() {
                freq.push(Some(*af));
            }
            line.clear();
        }

        let freq = freq.finish();
        Ok(freq)
    }
    pub fn infer_freq_from_data(geno: &Matrix<u8>, sites: &Sites, args: &Arguments) -> Matrix<f64> {
        // assert gentoeyps is still site oriented
        assert_eq!(geno.get_nrows(), sites.get_pos_slice().len());

        let mut freq = MatrixBuilder::<f64>::new(args.max_all as usize);
        let mut cnts = vec![0u32; args.max_all as usize];
        let mut total;

        for row in 0..geno.get_nrows() {
            cnts.clear();
            cnts.resize(args.max_all as usize, 0);
            total = 0;
            for allele in geno.get_row_iter(row).flatten() {
                total += 1;
                cnts[allele as usize] += 1;
            }
            for each in cnts.iter() {
                let af = Some(*each as f64 / total as f64);
                freq.push(af);
            }
        }

        freq.finish()
    }

    fn get_valid_pair_file(
        args: &Arguments,
        valid_samples: &Samples,
    ) -> Result<Vec<(u32, u32)>, Error> {
        let mut v = vec![];
        if let Some(good_file) = args.good_file.as_ref() {
            let mut buf = String::new();
            std::fs::File::open(good_file)
                .map(BufReader::new)
                .map_err(|e| Error::Io {
                    source: e,
                    file: Some(good_file.to_owned()),
                })?
                .read_to_string(&mut buf)
                .map_err(|e| Error::Io {
                    source: e,
                    file: Some(good_file.to_owned()),
                })?;
            let m = valid_samples.m();
            for line in buf.trim().split("\n") {
                let mut fields = line.split("\t");
                let s1 = fields.next().ok_or(ParseLineError::ReadColumnError("s1"))?;
                let s2 = fields.next().ok_or(ParseLineError::ReadColumnError("s2"))?;
                if m.contains_key(s1) && m.contains_key(s2) {
                    v.push((m[s1], m[s2]));
                }
            }
        } else {
            // two populatoin
            if valid_samples.pop2_nsam() > 0 {
                for i in 0..valid_samples.pop1_nsam() {
                    for j in valid_samples.pop1_nsam()..valid_samples.v().len() as u32 {
                        v.push((i, j))
                    }
                }
            }
            // one population
            else {
                let s1 = valid_samples.v();
                let n = s1.len() as u32;
                for i in 0..(n - 1) {
                    for j in (i + 1)..n {
                        v.push((i, j))
                    }
                }
            }
        }
        Ok(v)
    }

    /// get number of different alleles for all valide sites
    fn get_nall(freq1: &Matrix<f64>, freq2: Option<&Matrix<f64>>) -> Result<Vec<u8>, Error> {
        let mut v = vec![];
        // println!("{:?}", freq1.get_row_raw_slice(0));
        // println!("{:?}", freq2.unwrap().get_row_raw_slice(0));
        for i in 0..freq1.get_nrows() {
            let mut it1 =
                freq1
                    .get_row_iter(i)
                    .enumerate()
                    .map(|(iallele, x)| -> Result<bool, Error> {
                        Ok(x.ok_or(Error::MissingFrequency {
                            irow: i,
                            iallele,
                            value: x,
                        })? > 0.0)
                    });
            let n = match freq2 {
                Some(geno2) => {
                    let it2 = geno2.get_row_iter(i).enumerate().map(
                        |(iallele, x)| -> Result<bool, Error> {
                            Ok(x.ok_or(Error::MissingFrequency {
                                irow: i,
                                iallele,
                                value: x,
                            })? > 0.0)
                        },
                    );
                    it1.zip(it2)
                        .try_fold(0usize, |acc, (x, y)| -> Result<usize, Error> {
                            Ok(acc + (x? || y?) as usize)
                        })?
                }
                None => it1.try_fold(0usize, |acc, x| -> Result<usize, Error> {
                    Ok(acc + x? as usize)
                })?,
            };
            v.push(n as u8);
        }
        Ok(v)
    }

    /// Get the major allele for each site
    ///
    /// When one site has zero alleles, `Err(Error)` will be returned.
    fn get_major_all(
        freq1: &Matrix<f64>,
        freq2: Option<&Matrix<f64>>,
        nsam1_valid: usize,
        nsam2_valid: Option<usize>,
    ) -> Result<Vec<u8>, Error> {
        let mut v = vec![];
        let freq2 = freq2.unwrap_or(freq1);
        let nsam2_valid = nsam2_valid.unwrap_or(nsam1_valid) as f64;
        let nsam1_valid = nsam1_valid as f64;

        for ipos in 0..freq1.get_nrows() {
            let it = freq1.get_row_iter(ipos).map(|x| x.unwrap_or(0.0));
            let it2 = freq2.get_row_iter(ipos).map(|x| x.unwrap_or(0.0));

            let major_all = it
                .zip(it2)
                .map(|(a, b)| a * nsam1_valid + b * nsam2_valid)
                .enumerate()
                .try_fold(
                    (0, f64::MIN),
                    |(max_idx, max_val), (this_idx, this_val)| -> Result<(usize, f64), Error> {
                        if let std::cmp::Ordering::Greater = this_val
                            .partial_cmp(&max_val)
                            .ok_or(Error::PartialCmpIsNone("get_major_all"))?
                        {
                            Ok((this_idx, this_val))
                        } else {
                            Ok((max_idx, max_val))
                        }
                    },
                )?
                .0;

            // let major_all = it
            //     .zip(it2)
            //     .map(|(a, b)| a * nsam1_valid + b * nsam2_valid)
            //     .enumerate()
            //     .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            //     .ok_or(Error::EmptyIterator)?
            //     .0;
            v.push(major_all as u8);
            // if major_all > 1 {
            //     println!("{major_all}");
            // }
        }
        // println!("freq1.shape={},{}", freq1.get_nrows(), freq1.get_ncols());
        // println!("freq2.shape={},{}", freq2.get_nrows(), freq2.get_ncols());
        // println!("nsam1_valid={}, nsam2_valid={}", nsam1_valid, nsam2_valid);
        // println!("v.len={}", v.len());
        Ok(v)
    }
}

pub struct OutputFiles {
    pub seg_file: Arc<Mutex<BufWriter<File>>>,
    pub frac_file: Arc<Mutex<BufWriter<File>>>,
}

impl OutputFiles {
    pub fn new_from_args(
        args: &Arguments,
        buffer_size_segments: Option<usize>,
        buffer_size_frac: Option<usize>,
    ) -> Result<Self, Error> {
        let prefix = match args.output.as_ref() {
            Some(output) => output,
            None => &args.data_file1,
        };
        let seg_fn = format!("{prefix}.hmm.txt");
        let frac_fn = format!("{prefix}.hmm_fract.txt");
        use std::io::Write;

        let mut seg_file = match buffer_size_segments {
            Some(bfsz) => {
                std::fs::File::create(&seg_fn).map(|seg_fn| BufWriter::with_capacity(bfsz, seg_fn))
            }
            None => std::fs::File::create(&seg_fn).map(BufWriter::new),
        }
        .map_err(|e| Error::Io {
            source: e,
            file: Some(seg_fn.clone()),
        })?;

        let mut frac_file = match buffer_size_frac {
            Some(bfsz) => std::fs::File::create(&frac_fn)
                .map(|frac_fn| BufWriter::with_capacity(bfsz, frac_fn)),
            None => std::fs::File::create(&frac_fn).map(BufWriter::new),
        }
        .map_err(|e| Error::Io {
            source: e,
            file: Some(seg_fn.clone()),
        })?;

        // write header:
        writeln!(
            &mut frac_file,
            "sample1\tsample2\tN_informative_sites\tdiscordance\tlog_p\tN_fit_iteration\tN_generation\tN_state_transition\tseq_shared_best_traj\tfract_sites_IBD\tfract_vit_sites_IBD"
        )
        .map_err(|e| Error::Io {
            source: e,
            file: Some(frac_fn.clone()),
        })?;

        writeln!(
            &mut seg_file,
            "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp"
        )
        .map_err(|e| Error::Io {
            source: e,
            file: Some(seg_fn.clone()),
        })?;

        Ok(Self {
            seg_file: Arc::new(Mutex::new(seg_file)),
            frac_file: Arc::new(Mutex::new(frac_file)),
        })
    }
}

#[test]
fn read_inputdata() {
    let mut args = Arguments::new_for_test();
    args.freq_file1 = None;
    let _input = InputData::from_args(&args);
}

pub struct FracRecord<'a> {
    pub sample1: &'a str,
    pub sample2: &'a str,
    pub sum: usize,
    pub discord: f64,
    pub max_phi: f64,
    pub iter: usize,
    pub k_rec: f64,
    pub ntrans: usize,
    pub seq_ibd_ratio: f64,
    pub count_ibd_fb_ratio: f64,
    pub count_ibd_vit_ratio: f64,
}

pub struct SegRecord<'a> {
    pub sample1: &'a str,
    pub sample2: &'a str,
    pub chrname: &'a str,
    pub start_pos: u32,
    pub end_pos: u32,
    pub ibd: u8,
    pub n_snp: usize,
}

pub struct OutputBuffer<'a> {
    seg_file: Arc<Mutex<BufWriter<File>>>,
    frac_file: Arc<Mutex<BufWriter<File>>>,
    segs: SmallVec<[SegRecord<'a>; 5]>,
    fracs: SmallVec<[FracRecord<'a>; 1]>,
}

impl<'a> OutputBuffer<'a> {
    pub fn new(out: &OutputFiles, segs_capacity: usize, fracs_capacity: usize) -> Self {
        Self {
            seg_file: Arc::clone(&out.seg_file),
            frac_file: Arc::clone(&out.frac_file),
            segs: SmallVec::<[SegRecord<'a>; 5]>::with_capacity(segs_capacity),
            fracs: SmallVec::<[FracRecord<'a>; 1]>::with_capacity(fracs_capacity),
        }
    }

    pub fn add_seg(&mut self, seg: SegRecord<'a>) -> Result<(), Error> {
        if self.segs.len() == self.segs.capacity() {
            self.flush_segs()?;
        }
        self.segs.push(seg);
        Ok(())
    }

    pub fn add_frac(&mut self, frac: FracRecord<'a>) -> Result<(), Error> {
        if self.fracs.len() == self.fracs.len() {
            self.flush_frac()?;
        }
        self.fracs.push(frac);
        Ok(())
    }

    pub fn flush_segs(&mut self) -> Result<(), Error> {
        use std::io::Write;
        let mut file = self
            .seg_file
            .lock()
            .map_err(|_| Error::LockError("flush_segs"))?;
        for seg in self.segs.iter() {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                seg.sample1,
                seg.sample2,
                seg.chrname,
                seg.start_pos,
                seg.end_pos,
                seg.ibd,
                seg.n_snp,
            )
            .unwrap();
        }
        self.segs.clear();
        Ok(())
    }

    pub fn flush_frac(&mut self) -> Result<(), Error> {
        use std::io::Write;
        let mut file = self
            .frac_file
            .lock()
            .map_err(|_| Error::LockError("flush_frac"))?;
        for frac in self.fracs.iter() {
            writeln!(
                file,
                "{}\t{}\t{}\t{:.4}\t{:0.5e}\t{}\t{:.3}\t{}\t{:.5}\t{:.5}\t{:.5}",
                frac.sample1,
                frac.sample2,
                frac.sum,
                frac.discord,
                frac.max_phi,
                frac.iter,
                frac.k_rec,
                frac.ntrans,
                frac.seq_ibd_ratio,
                frac.count_ibd_fb_ratio,
                frac.count_ibd_vit_ratio,
            )
            .unwrap();
        }
        self.fracs.clear();
        Ok(())
    }
}
