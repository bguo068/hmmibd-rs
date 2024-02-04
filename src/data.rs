use smallvec::SmallVec;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read},
    path::Path,
    sync::{Arc, Mutex},
};

use crate::{
    args::Arguments, genome::*, genome2::*, gmap::*, matrix::*, samples::Samples, sites::Sites,
};

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
                    chunk_pairs.push((i as u32, j as u32));
                }
            }
        } else {
            // same population
            for i in 0..num_chunk1 {
                for j in i..num_chunk1 {
                    chunk_pairs.push((i as u32, j as u32));
                }
            }
        }
        chunk_pairs
    }

    pub fn clone_inputdata_for_chunkpair(&self, chunkpair: (u32, u32)) -> Self {
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
        let pairs = Self::get_valid_pair_file(&args, &samples);
        Self {
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
        }
    }

    pub fn from_args(args: &Arguments) -> Self {
        let valid_samples = Samples::from_args(&args);
        let rec_rate = args.rec_args.rec_rate;
        let min_snp_sep = args.min_snp_sep;

        let (genome, mut geno1, geno2, sites) = match args.rec_args.genome.as_ref() {
            Some(genome_toml_path) => {
                let ginfo = GenomeInfo::from_toml_file(genome_toml_path);
                let genome = Genome::from_genome_info(&ginfo);
                let gmap = GeneticMap::from_genome_info(&ginfo);
                let (geno1, sites) = Self::read_data_file_with_ginfo_and_gmap(
                    &args.data_file1,
                    &valid_samples,
                    min_snp_sep,
                    &ginfo,
                    &gmap,
                    &genome,
                );
                let geno2 = args.data_file2.as_ref().map(|data_file2| {
                    let (geno2, sites2) = Self::read_data_file_with_ginfo_and_gmap(
                        data_file2,
                        &valid_samples,
                        min_snp_sep,
                        &ginfo,
                        &gmap,
                        &genome,
                    );
                    assert!(sites.is_identical(&sites2));
                    geno2
                });
                (genome, geno1, geno2, sites)
            }
            None => {
                let (genome, geno1, sites) =
                    Self::read_data_file(&args.data_file1, &valid_samples, min_snp_sep, rec_rate);
                let geno2 = args.data_file2.as_ref().map(|data_file2| {
                    let (genome2, geno2, sites2) =
                        Self::read_data_file(data_file2, &valid_samples, min_snp_sep, rec_rate);
                    assert!(genome.is_identical(&genome2));
                    assert!(sites.is_identical(&sites2));
                    geno2
                });
                (genome, geno1, geno2, sites)
            }
        };
        let nsam1_valid = geno1.get_ncols();
        let nsam2_valid = geno2.as_ref().map(|geno2| geno2.get_ncols());
        let freq1 = if let Some(freq_file1) = args.freq_file1.as_ref() {
            Self::read_freq_file(freq_file1, &args, &genome, &sites)
        } else {
            Self::infer_freq_from_data(&geno1, &sites, &args)
        };

        let freq2 = match (geno2.as_ref(), args.freq_file2.as_ref()) {
            (Some(_), Some(freq_file2)) => {
                Some(Self::read_freq_file(freq_file2, &args, &genome, &sites))
            }
            (Some(geno2), None) => Some(Self::infer_freq_from_data(&geno2, &sites, &args)),
            _ => None,
        };
        let pairs = Self::get_valid_pair_file(&args, &valid_samples);
        let nall = Self::get_nall(&freq1, freq2.as_ref());
        let majall = Self::get_major_all(&freq1, freq2.as_ref(), nsam1_valid, nsam2_valid);

        // sample oriented
        geno1.transpose();
        if let Some(mut geno2) = geno2 {
            geno2.transpose();
            geno1.merge(&geno2);
        }

        Self {
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
        }
    }
    fn read_data_file_with_ginfo_and_gmap(
        data_file: impl AsRef<Path>,
        valid_samples: &Samples,
        min_snp_sep: u32,
        ginfo: &GenomeInfo,
        gmap: &GeneticMap,
        genome: &Genome,
    ) -> (Matrix<u8>, Sites) {
        let mut sites = Sites::new();

        let mut line = String::with_capacity(100000);
        let mut f = std::fs::File::open(data_file.as_ref())
            .map(BufReader::new)
            .unwrap();

        // get all sample names
        f.read_line(&mut line).unwrap();
        let samples: Vec<_> = line
            .trim()
            .split('\t')
            .skip(2)
            .map(|x| x.to_owned())
            .collect();
        line.clear();

        let n_valid_samples = samples
            .iter()
            .filter(|s| valid_samples.m().contains_key(*s))
            .count();

        let mut geno1 = MatrixBuilder::<u8>::new(n_valid_samples);
        let mut last_chrname = String::new();
        let mut last_pos = 0;
        while f.read_line(&mut line).unwrap() != 0 {
            let mut fields = line.trim().split("\t");
            let chrname = fields.next().unwrap();
            let pos: u32 = fields.next().unwrap().parse().unwrap();

            // println!("pos: {pos}, last_pos: {last_pos}");
            if (chrname == last_chrname) && (last_pos + min_snp_sep > pos) {
                line.clear();
                continue;
            } else {
                last_chrname.clear();
                last_chrname.push_str(chrname);
                last_pos = pos;
            }

            let chrid = ginfo.idx[chrname];
            let gw_pos = ginfo.to_gw_pos(chrid, pos);
            let gw_pos_cm = gmap.get_cm(gw_pos);

            sites.add(gw_pos, gw_pos_cm);

            let mut cnt = 0;
            fields.enumerate().for_each(|(i, field)| {
                if valid_samples.m().contains_key(&samples[i]) {
                    let allel = match field.parse::<i8>().unwrap() {
                        -1 => None,
                        x => Some(x as u8),
                    };
                    geno1.push(allel);
                    cnt += 1;
                }
            });
            assert_eq!(cnt, n_valid_samples);
            line.clear();
        }

        let geno1 = geno1.finish();
        sites.finish(&genome);

        (geno1, sites)
    }

    fn read_data_file(
        data_file: impl AsRef<Path>,
        valid_samples: &Samples,
        min_snp_sep: u32,
        rec_rate: f64,
    ) -> (Genome, Matrix<u8>, Sites) {
        let mut sites = Sites::new();
        let mut gbuilder = GenomeBuilder::new();

        let mut line = String::with_capacity(100000);
        let mut f = std::fs::File::open(data_file.as_ref())
            .map(BufReader::new)
            .unwrap();

        // get all sample names
        f.read_line(&mut line).unwrap();
        let samples: Vec<_> = line
            .trim()
            .split('\t')
            .skip(2)
            .map(|x| x.to_owned())
            .collect();
        line.clear();

        let n_valid_samples = samples
            .iter()
            .filter(|s| valid_samples.m().contains_key(*s))
            .count();

        let mut geno1 = MatrixBuilder::<u8>::new(n_valid_samples);
        let mut last_chrname = String::new();
        let mut last_pos = 0;
        while f.read_line(&mut line).unwrap() != 0 {
            let mut fields = line.trim().split("\t");
            let chrname = fields.next().unwrap();
            let pos: u32 = fields.next().unwrap().parse().unwrap();

            // println!("pos: {pos}, last_pos: {last_pos}");
            if (chrname == last_chrname) && (last_pos + min_snp_sep > pos) {
                line.clear();
                continue;
            } else {
                last_chrname.clear();
                last_chrname.push_str(chrname);
                last_pos = pos;
            }

            let (_, gw_pos) = gbuilder.encode_pos(chrname, pos);

            sites.add(gw_pos, gw_pos as f64 * rec_rate * 100.0);

            let mut cnt = 0;
            fields.enumerate().for_each(|(i, field)| {
                if valid_samples.m().contains_key(&samples[i]) {
                    let allel = match field.parse::<i8>().unwrap() {
                        -1 => None,
                        x => Some(x as u8),
                    };
                    geno1.push(allel);
                    cnt += 1;
                }
            });
            assert_eq!(cnt, n_valid_samples);
            line.clear();
        }

        let genome = gbuilder.finish();
        let geno1 = geno1.finish();
        sites.finish(&genome);

        (genome, geno1, sites)
    }
    fn read_freq_file(
        freq_file: impl AsRef<Path>,
        args: &Arguments,
        genome: &Genome,
        sites: &Sites,
    ) -> Matrix<f64> {
        let mut f = std::fs::File::open(freq_file.as_ref())
            .map(BufReader::new)
            .unwrap();

        let mut freq = MatrixBuilder::<f64>::new(args.max_all as usize);
        let mut v = Vec::new();

        let mut line = String::with_capacity(100000);
        line.clear();

        let mut last_chrname = String::new();
        let mut last_pos = 0;
        while f.read_line(&mut line).unwrap() != 0 {
            let mut fields = line.trim().split("\t");
            let chrname = fields.next().unwrap();
            let pos: u32 = fields.next().unwrap().parse().unwrap();
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
            fields.enumerate().for_each(|(i, field)| {
                let af = field.parse::<f64>().unwrap();
                v[i] = af;
            });
            for af in v.iter() {
                freq.push(Some(*af));
            }
            line.clear();
        }

        let freq = freq.finish();
        freq
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
            for allele in geno.get_row_iter(row) {
                match allele {
                    Some(allele) => {
                        total += 1;
                        cnts[allele as usize] += 1;
                    }
                    None => {}
                };
            }
            for each in cnts.iter() {
                let af = Some(*each as f64 / total as f64);
                freq.push(af);
            }
        }

        let freq = freq.finish();
        // println!("geno.shape={}x{}", geno.get_nrows(), geno.get_ncols());
        // println!("freq.shape={}x{}", freq.get_nrows(), freq.get_ncols());
        freq
    }

    fn get_valid_pair_file(args: &Arguments, valid_samples: &Samples) -> Vec<(u32, u32)> {
        let mut v = vec![];
        if let Some(good_file) = args.good_file.as_ref() {
            let mut buf = String::new();
            std::fs::File::open(good_file)
                .map(BufReader::new)
                .unwrap()
                .read_to_string(&mut buf)
                .unwrap();
            let m = valid_samples.m();
            for line in buf.trim().split("\n") {
                let mut fields = line.split("\t");
                let s1 = fields.next().unwrap();
                let s2 = fields.next().unwrap();
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
        v
    }

    fn get_nall(freq1: &Matrix<f64>, freq2: Option<&Matrix<f64>>) -> Vec<u8> {
        let mut v = vec![];
        // println!("{:?}", freq1.get_row_raw_slice(0));
        // println!("{:?}", freq2.unwrap().get_row_raw_slice(0));
        for i in 0..freq1.get_nrows() {
            let it1 = freq1.get_row_iter(i).map(|x| x.unwrap() > 0.0);
            let n = match freq2 {
                Some(geno2) => {
                    let it2 = geno2.get_row_iter(i).map(|x| x.unwrap() > 0.0);
                    it1.zip(it2).filter(|(x, y)| *x || *y).count()
                }
                None => it1.filter(|x| *x).count(),
            };
            v.push(n as u8);
        }
        v
    }
    fn get_major_all(
        freq1: &Matrix<f64>,
        freq2: Option<&Matrix<f64>>,
        nsam1_valid: usize,
        nsam2_valid: Option<usize>,
    ) -> Vec<u8> {
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
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .unwrap()
                .0;
            v.push(major_all as u8);
            // if major_all > 1 {
            //     println!("{major_all}");
            // }
        }
        // println!("freq1.shape={},{}", freq1.get_nrows(), freq1.get_ncols());
        // println!("freq2.shape={},{}", freq2.get_nrows(), freq2.get_ncols());
        // println!("nsam1_valid={}, nsam2_valid={}", nsam1_valid, nsam2_valid);
        // println!("v.len={}", v.len());
        v
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
    ) -> Self {
        let prefix = match args.output.as_ref() {
            Some(output) => output,
            None => &args.data_file1,
        };
        let seg_fn = format!("{prefix}.hmm.txt");
        let frac_fn = format!("{prefix}.hmm_fract.txt");
        use std::io::Write;

        let mut seg_file = match buffer_size_segments {
            Some(bfsz) => std::fs::File::create(&seg_fn)
                .map(|seg_fn| BufWriter::with_capacity(bfsz, seg_fn))
                .unwrap(),
            None => std::fs::File::create(&seg_fn).map(BufWriter::new).unwrap(),
        };

        let mut frac_file = match buffer_size_frac {
            Some(bfsz) => std::fs::File::create(&frac_fn)
                .map(|frac_fn| BufWriter::with_capacity(bfsz, frac_fn))
                .unwrap(),
            None => std::fs::File::create(&frac_fn).map(BufWriter::new).unwrap(),
        };

        // write header:
        write!(
            &mut frac_file,
            "sample1\tsample2\tN_informative_sites\tdiscordance\tlog_p\tN_fit_iteration\tN_generation\tN_state_transition\tseq_shared_best_traj\tfract_sites_IBD\tfract_vit_sites_IBD\n"
        )
        .unwrap();

        write!(
            &mut seg_file,
            "sample1\tsample2\tchr\tstart\tend\tdifferent\tNsnp\n"
        )
        .unwrap();

        Self {
            seg_file: Arc::new(Mutex::new(seg_file)),
            frac_file: Arc::new(Mutex::new(frac_file)),
        }
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

    pub fn add_seg(&mut self, seg: SegRecord<'a>) {
        if self.segs.len() == self.segs.capacity() {
            self.flush_segs();
        }
        self.segs.push(seg);
    }

    pub fn add_frac(&mut self, frac: FracRecord<'a>) {
        if self.fracs.len() == self.fracs.len() {
            self.flush_frac();
        }
        self.fracs.push(frac);
    }

    pub fn flush_segs(&mut self) {
        use std::io::Write;
        let mut file = self.seg_file.lock().unwrap();
        for seg in self.segs.iter() {
            write!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
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
    }

    pub fn flush_frac(&mut self) {
        use std::io::Write;
        let mut file = self.frac_file.lock().unwrap();
        for frac in self.fracs.iter() {
            write!(
                file,
                "{}\t{}\t{}\t{:.4}\t{:0.5e}\t{}\t{:.3}\t{}\t{:.5}\t{:.5}\t{:.5}\n",
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
    }
}
