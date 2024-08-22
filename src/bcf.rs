use bcf_reader::*;
use itertools::{EitherOrBoth, Itertools};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    io::{BufReader, BufWriter},
};

use crate::{
    matrix::{Matrix, MatrixBuilder},
    samples::Samples,
    sites::SiteInfoRaw,
};

#[derive(Default, Serialize, Deserialize, Debug, Clone)]
pub struct BcfFilterArgs {
    pub min_depth: u32,
    pub min_ratio: f32,
    pub min_maf: f32,
    pub min_r1_r2: f32,
    pub filter_column_pass_only: bool,
    pub major_minor_alleles_must_be_snps: bool,
    pub min_site_nonmissing: f32,
    pub nonmissing_rates: Vec<(f32, f32)>,
    pub target_samples: Option<std::path::PathBuf>,
}

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("bcf reader error: {0:?}")]
    BcfReaderError(#[from] bcf_reader::Error),

    #[error("io error: {0:?}")]
    Io(#[from] std::io::Error),

    #[error("toml deserialize error: {0:?}")]
    TomlDeserializeError(#[from] toml::de::Error),

    #[error("genotype is empty: {file}:{line}")]
    GenotypeEmpty { file: &'static str, line: u32 },

    #[error("bcf genotype object is not ready")]
    BcfGenotypeNotReady,

    #[error("bcf genotype object is not empty")]
    BcfGenotypeNotEmpty,

    #[error("bincode error: {0:?}")]
    BincodeErr(#[from] bincode::Error),
}

impl BcfFilterArgs {
    pub fn new_from_toml_file(dom_gt_config_path: &str) -> Result<Self> {
        let toml_str = std::fs::read_to_string(dom_gt_config_path)?;
        let bcf_filter_args: BcfFilterArgs = toml::from_str(&toml_str)?;
        Ok(bcf_filter_args)
    }
    pub fn new_from_builtin() -> Result<Self> {
        let toml_str = r#"
min_depth = 5
min_ratio = 0.7
min_r1_r2 = 3.0
min_maf = 0.01
min_site_nonmissing = 0.3
filter_column_pass_only = true
major_minor_alleles_must_be_snps = true
nonmissing_rates= [
    [0.31, 0.155],
    [0.35, 0.175],
    [0.41, 0.205],
    [0.45, 0.225],
    [0.51, 0.255],
    [0.55, 0.275],
    [0.61, 0.305],
    [0.65, 0.325],
    [0.71, 0.355],
    [0.75, 0.375],
    [0.81, 0.405],
    [0.85, 0.45],
    [0.85, 0.50],
    [0.85, 0.55],
    [0.85, 0.60],
    [0.85, 0.65],
    [0.90, 0.70],
    [0.90, 0.80],
 ]
    "#;

        let dga: BcfFilterArgs = toml::from_str(toml_str)?;
        Ok(dga)
    }
}

#[derive(Default, Deserialize, Debug, Serialize)]
pub struct BcfGenotype {
    chr_vec: Vec<i32>,
    pos_vec: Vec<i32>,
    gt_vec: Vec<i8>,
    sample_vec: Vec<String>,
    chrname_map: HashMap<String, usize>,
    selected_samples: Vec<usize>,
    selected_sites: Vec<usize>,
    nsam_targeted: usize,
    args: BcfFilterArgs,
    ready_to_use: bool,
}

impl BcfGenotype {
    fn new(bcffilter_args: &BcfFilterArgs) -> Self {
        Self {
            chr_vec: vec![],
            pos_vec: vec![],
            gt_vec: vec![],
            sample_vec: vec![],
            chrname_map: HashMap::new(),
            selected_samples: vec![],
            selected_sites: vec![],
            nsam_targeted: 0,
            args: bcffilter_args.clone(),
            ready_to_use: false,
        }
    }
    fn read_dom(&mut self, bcf_fname: &str) -> Result<Header> {
        if self.ready_to_use {
            return Err(Error::BcfGenotypeNotEmpty);
        }
        let mut reader = get_bcf_gzip_reader(bcf_fname)?;
        let header = reader.read_header()?;
        // .map_err(|e| bcf_reader::Error::ParseHeaderError(e))?;
        let mut record = Record::default();
        let ad_key = header.get_idx_from_dictionary_str("FORMAT", "AD").unwrap();
        let nsam = header.get_samples().len();
        let mut chrname_map = HashMap::<String, usize>::new();
        for (id, dict) in header.dict_contigs().iter() {
            let chrname = dict["ID"].to_owned();
            chrname_map.insert(chrname, *id);
        }
        self.chrname_map = chrname_map;

        let is_sample_targeted =
            get_targeted_sample_indicator(&self.args, header.get_samples().as_slice())?;

        let nsam_in_targets: usize = is_sample_targeted.iter().map(|e| *e as usize).sum();
        self.nsam_targeted = nsam_in_targets;
        self.sample_vec.extend(
            header
                .get_samples()
                .iter()
                .zip(is_sample_targeted.iter())
                .filter(|(_sname, is_target)| **is_target)
                .map(|(sname, _is_target)| sname.to_owned()),
        );
        let mut nrec = 0;
        let mut nvalid = 0;

        let mut indv_ad = Vec::<(usize, u32)>::new(); // (allele_idx, ad)
        let mut allele_counts = Vec::<(usize, u32)>::new(); // (allele_idx, AC)
        while reader.read_record(&mut record).is_ok() {
            let nallele = record.n_allele() as usize;

            // not segrating or if not pass in filter columns
            if (nallele <= 1)
                || (self.args.filter_column_pass_only
                    && (!is_variant_pass_in_filter_columns(&header, &record)?))
            {
                if nrec % 1000 == 0 {
                    eprint!("\r{nrec}\t{nvalid}");
                }
                nrec += 1;
                continue;
            }

            let mut site_nonmiss_counter = 0;
            let chunks = record.fmt_field(ad_key).chunks(record.n_allele() as usize);
            for (mut nv_indv, _) in chunks
                .into_iter()
                .zip(is_sample_targeted.iter())
                .filter(|(_, is_target)| **is_target)
            {
                indv_ad.clear();
                nv_indv.try_for_each(|nv_res| -> Result<()> {
                    let val = nv_res?
                        .int_val()
                        .ok_or(bcf_reader::Error::NumericaValueEmptyInt)?;
                    let idx = indv_ad.len();
                    indv_ad.push((idx, val));
                    Ok(())
                })?;
                // indv_ad.extend(nv_indv.map(|nv| nv?.int_val().unwrap()).enumerate());
                indv_ad.sort_by_key(|x| u32::MAX - x.1);
                let mut dom_allele = -1i8;
                let total = indv_ad.iter().map(|x| x.1).sum::<u32>();
                let major = indv_ad[0].1;
                let minor = indv_ad[1].1;
                let r1 = major as f32 / total as f32;
                let r2_r1 = minor as f32 / major as f32;
                if (total > self.args.min_depth)
                    & (r1 >= self.args.min_ratio)
                    & (r2_r1 < 1.0 / self.args.min_r1_r2)
                {
                    dom_allele = indv_ad[0].0 as i8;
                    site_nonmiss_counter += 1;
                }
                self.gt_vec.push(dom_allele);
            }
            let non_missing_rate = site_nonmiss_counter as f32 / nsam as f32;

            let maf = get_maf_and_sorted_allele_count(
                &self.gt_vec[(self.gt_vec.len() - nsam_in_targets)..],
                nallele,
                &mut allele_counts,
            );
            let is_major_minor_snps = is_major_and_minor_allele_snps(&record, &allele_counts);

            // discard or keep current sites
            if (non_missing_rate < self.args.min_site_nonmissing)
                || (maf < self.args.min_maf)
                || (self.args.major_minor_alleles_must_be_snps && (!is_major_minor_snps))
            {
                // remove alleles for this site
                self.gt_vec.resize(self.gt_vec.len() - nsam_in_targets, 0);
            } else {
                self.chr_vec.push(record.chrom());
                self.pos_vec.push(record.pos());
                nvalid += 1;
            }

            nrec += 1;
            // if nrec % 1000 == 0 {
            //     eprintln!("\r{nrec}\t{nvalid}");
            // }
        }
        self.selected_samples.extend(0..self.sample_vec.len());
        self.selected_sites.extend(0..self.pos_vec.len());
        Ok(header)
    }

    fn read_first_ploidy(&mut self, bcf_fname: &str) -> Result<Header> {
        if self.ready_to_use {
            return Err(Error::BcfGenotypeNotEmpty);
        }
        let mut reader = get_bcf_gzip_reader(bcf_fname)?;
        let header = reader.read_header()?;
        // .map_err(|e| bcf_reader::Error::ParseHeaderError(e))?;
        let mut record = Record::default();
        let gt_key = header.get_idx_from_dictionary_str("FORMAT", "GT").unwrap();
        let nsam = header.get_samples().len();
        let mut chrname_map = HashMap::<String, usize>::new();
        for (id, dict) in header.dict_contigs().iter() {
            let chrname = dict["ID"].to_owned();
            chrname_map.insert(chrname, *id);
        }
        self.chrname_map = chrname_map;

        // -- targeted sample indicator
        let is_sample_targeted =
            get_targeted_sample_indicator(&self.args, header.get_samples().as_slice())?;
        let nsam_in_targets: usize = is_sample_targeted.iter().map(|e| *e as usize).sum();
        self.nsam_targeted = nsam_in_targets;

        // --
        self.sample_vec.extend(
            header
                .get_samples()
                .iter()
                .zip(is_sample_targeted.iter())
                .filter(|(_sname, is_target)| **is_target)
                .map(|(sname, _is_target)| sname.to_owned()),
        );
        let mut nrec = 0;
        let mut nvalid = 0;

        let mut allele_counts = Vec::<(usize, u32)>::new(); // (allele_idx, AC)
        while reader.read_record(&mut record).is_ok() {
            let nploidy = record.fmt_gt(&header).count() / nsam;
            let nallele = record.n_allele() as usize;

            // not segrating or if not pass in filter columns
            if (nallele <= 1)
                || (self.args.filter_column_pass_only
                    && (!is_variant_pass_in_filter_columns(&header, &record)?))
            {
                if nrec % 1000 == 0 {
                    eprint!("\r{nrec}\t{nvalid}");
                }
                nrec += 1;
                continue;
            }

            let mut site_nonmiss_counter = 0;
            let chunks = record.fmt_field(gt_key).chunks(nploidy);
            for (nv_indv, _) in chunks
                .into_iter()
                .zip(is_sample_targeted.iter())
                .filter(|(_, is_target)| **is_target)
            {
                let mut first_ploidy_allele = -1i8;
                nv_indv
                    .enumerate()
                    .try_for_each(|(iploidy, nv_res)| -> Result<()> {
                        let (_noploidy, dot, _phased, allele) = nv_res?.gt_val();
                        match iploidy {
                            0 => {
                                if !dot {
                                    first_ploidy_allele = allele as i8;
                                    site_nonmiss_counter += 1;
                                }
                            }
                            1 => {
                                if (!dot) && _noploidy {
                                    return Err(Error::BcfReaderError(bcf_reader::Error::Other(
                                        "genotype should be phased".to_owned(),
                                    )));
                                }
                            }
                            _ => {}
                        }
                        Ok(())
                    })?;
                self.gt_vec.push(first_ploidy_allele);
            }
            let non_missing_rate = site_nonmiss_counter as f32 / nsam as f32;

            let maf = get_maf_and_sorted_allele_count(
                &self.gt_vec[(self.gt_vec.len() - nsam_in_targets)..],
                nallele,
                &mut allele_counts,
            );
            let is_major_minor_snps = is_major_and_minor_allele_snps(&record, &allele_counts);

            // discard or keep current sites
            if (non_missing_rate < self.args.min_site_nonmissing)
                || (maf < self.args.min_maf)
                || (self.args.major_minor_alleles_must_be_snps && (!is_major_minor_snps))
            {
                // remove alleles for this site
                self.gt_vec.resize(self.gt_vec.len() - nsam_in_targets, 0);
            } else {
                self.chr_vec.push(record.chrom());
                self.pos_vec.push(record.pos());
                nvalid += 1;
            }

            nrec += 1;
            // if nrec % 1000 == 0 {
            //     eprintln!("\r{nrec}\t{nvalid}");
            // }
        }
        self.selected_samples.extend(0..self.sample_vec.len());
        self.selected_sites.extend(0..self.pos_vec.len());
        Ok(header)
    }
    fn read_each_ploidy(&mut self, bcf_fname: &str) -> Result<Header> {
        if self.ready_to_use {
            return Err(Error::BcfGenotypeNotEmpty);
        }
        let mut reader = get_bcf_gzip_reader(bcf_fname)?;
        let header = reader.read_header()?;
        // .map_err(|e| bcf_reader::Error::ParseHeaderError(e))?;
        let mut record = Record::default();
        let gt_key = header.get_idx_from_dictionary_str("FORMAT", "GT").unwrap();
        let nsam = header.get_samples().len();
        let mut chrname_map = HashMap::<String, usize>::new();
        for (id, dict) in header.dict_contigs().iter() {
            let chrname = dict["ID"].to_owned();
            chrname_map.insert(chrname, *id);
        }
        self.chrname_map = chrname_map;

        // -- targeted sample indicator
        let is_sample_targeted =
            get_targeted_sample_indicator(&self.args, header.get_samples().as_slice())?;
        let nsam_in_targets: usize = is_sample_targeted.iter().map(|e| *e as usize).sum();
        self.nsam_targeted = nsam_in_targets;

        // --
        let mut nrec = 0;
        let mut nvalid = 0;

        let mut allele_counts = Vec::<(usize, u32)>::new(); // (allele_idx, AC)
        let mut nploidy_all: Option<usize> = None;
        while reader.read_record(&mut record).is_ok() {
            // not segrating or if not pass in filter columns
            let nallele = record.n_allele() as usize;
            if (nallele <= 1)
                || (self.args.filter_column_pass_only
                    && (!is_variant_pass_in_filter_columns(&header, &record)?))
            {
                if nrec % 1000 == 0 {
                    eprint!("\r{nrec}\t{nvalid}");
                }
                nrec += 1;
                continue;
            }

            // check ploidy consistency
            let nploidy = record.fmt_gt(&header).count() / nsam;
            if nploidy_all.is_none() {
                nploidy_all = Some(nploidy);
            } else if nploidy_all != Some(nploidy) {
                return Err(Error::BcfReaderError(bcf_reader::Error::Other(
                    "number of ploidy inconsistent across sites".to_owned(),
                )));
            }

            let mut site_nonmiss_counter = 0;
            let chunks = record.fmt_field(gt_key).chunks(nploidy);
            for (nv_indv, _) in chunks
                .into_iter()
                .zip(is_sample_targeted.iter())
                .filter(|(_, is_target)| **is_target)
            {
                nv_indv
                    .enumerate()
                    .try_for_each(|(iploidy, nv_res)| -> Result<()> {
                        let mut first_ploidy_allele = -1i8;
                        let (noploidy, dot, _phased, allele) = nv_res?.gt_val();
                        if !dot {
                            first_ploidy_allele = allele as i8;
                            site_nonmiss_counter += 1;
                        }
                        self.gt_vec.push(first_ploidy_allele);
                        // checking phasing
                        if (iploidy != 0) && (!dot) && noploidy {
                            return Err(Error::BcfReaderError(bcf_reader::Error::Other(
                                "genotype should be phased".to_owned(),
                            )));
                        }
                        Ok(())
                    })?;
            }
            let non_missing_rate = site_nonmiss_counter as f32 / nsam as f32;

            let maf = get_maf_and_sorted_allele_count(
                &self.gt_vec[(self.gt_vec.len() - nsam_in_targets * nploidy)..],
                nallele,
                &mut allele_counts,
            );
            let is_major_minor_snps = is_major_and_minor_allele_snps(&record, &allele_counts);

            // discard or keep current sites
            if (non_missing_rate < self.args.min_site_nonmissing)
                || (maf < self.args.min_maf)
                || (self.args.major_minor_alleles_must_be_snps && (!is_major_minor_snps))
            {
                // remove alleles for this site
                self.gt_vec
                    .resize(self.gt_vec.len() - nsam_in_targets * nploidy, 0);
            } else {
                self.chr_vec.push(record.chrom());
                self.pos_vec.push(record.pos());
                nvalid += 1;
            }

            nrec += 1;
            // if nrec % 1000 == 0 {
            //     eprintln!("\r{nrec}\t{nvalid}");
            // }
        }

        if let Some(nploidy) = nploidy_all {
            header
                .get_samples()
                .iter()
                .zip(is_sample_targeted.iter())
                .filter(|(_sname, is_target)| **is_target)
                .for_each(|(sname, _is_target)| {
                    for iploidy in 0..nploidy {
                        self.sample_vec.push(format!("{sname}___{iploidy}"));
                    }
                });
        }
        self.selected_samples.extend(0..self.sample_vec.len());
        self.selected_sites.extend(0..self.pos_vec.len());
        Ok(header)
    }

    fn filter_by_missingness(&mut self, site_nonmissing_rate: f32, sample_nonmissing_rate: f32) {
        // swap out to avoid ownership issue
        let mut sites = vec![];
        let mut samples = vec![];
        std::mem::swap(&mut sites, &mut self.selected_sites);
        std::mem::swap(&mut samples, &mut self.selected_samples);

        let mut new_good_sites = vec![];
        for row in sites.iter() {
            let mut cnt = 0;
            for col in samples.iter() {
                if self.gt_vec[*row * self.nsam_targeted + col] != -1 {
                    cnt += 1;
                };
            }
            if cnt as f32 / samples.len() as f32 > site_nonmissing_rate {
                new_good_sites.push(*row);
            }
        }
        sites = new_good_sites;

        let mut new_good_samples = vec![];
        let mut count_ind = vec![0usize; samples.len()];
        for row in sites.iter() {
            for (icol, col) in samples.iter().enumerate() {
                if self.gt_vec[*row * self.nsam_targeted + col] != -1 {
                    count_ind[icol] += 1;
                };
            }
        }
        for (idx, cnt) in count_ind.iter().enumerate() {
            if *cnt as f32 / sites.len() as f32 > sample_nonmissing_rate {
                new_good_samples.push(samples[idx]);
            }
        }
        samples = new_good_samples;

        self.selected_samples = samples;
        self.selected_sites = sites;
    }

    fn calc_nonmiss(&self) -> f32 {
        let mut cnt = 0;
        let mut cnt2 = 0;
        for row in self.selected_sites.iter() {
            for col in self.selected_samples.iter() {
                if self.gt_vec[*row * self.nsam_targeted + col] != -1 {
                    cnt += 1;
                };
                cnt2 += 1;
            }
        }
        cnt as f32 / cnt2 as f32
    }

    /// only keep sites and samples that are selected
    fn consolidate(&self) -> Self {
        let mut new_dg = Self::new(&self.args);

        let nsites = self.selected_sites.len();
        new_dg.selected_sites.extend(0..nsites);
        new_dg
            .pos_vec
            .extend(self.selected_sites.iter().map(|i| self.pos_vec[*i]));
        new_dg
            .chr_vec
            .extend(self.selected_sites.iter().map(|i| self.chr_vec[*i]));

        let nsam = self.selected_samples.len();
        new_dg.selected_samples.extend(0..nsam);
        new_dg.sample_vec.extend(
            self.selected_samples
                .iter()
                .map(|i| self.sample_vec[*i].to_owned()),
        );
        new_dg.nsam_targeted = nsam;

        new_dg.gt_vec.reserve(nsites * nsam);
        for site in self.selected_sites.iter() {
            for sample in self.selected_samples.iter() {
                new_dg
                    .gt_vec
                    .push(self.gt_vec[*site * self.nsam_targeted + *sample]);
            }
        }

        new_dg.chrname_map = self.chrname_map.clone();
        new_dg.ready_to_use = true;

        new_dg
    }

    pub fn new_from_processing_bcf(
        bcf_read_mode: &crate::args::BcfReadMode,
        bcf_filter_args: &BcfFilterArgs,
        bcf_path: &str,
    ) -> Result<Self> {
        use crate::args::BcfReadMode::*;
        let mut bcf_gt = Self::new(bcf_filter_args);

        match bcf_read_mode {
            DominantAllele => {
                bcf_gt.read_dom(bcf_path)?;
            }
            FirstPloidy => {
                bcf_gt.read_first_ploidy(bcf_path)?;
            }
            EachPloidy => {
                bcf_gt.read_each_ploidy(bcf_path)?;
            }
        };

        println!(
            "before filtering: n good site = {}, n good samples = {}",
            bcf_gt.selected_sites.len(),
            bcf_gt.selected_samples.len()
        );

        let rates = bcf_gt.args.nonmissing_rates.clone();
        for (min_site_nonmiss_rate, min_sample_nonmiss_rate) in rates {
            bcf_gt.filter_by_missingness(min_site_nonmiss_rate, min_sample_nonmiss_rate);
            let overall_nonmiss = bcf_gt.calc_nonmiss();
            println!(
                "min_site_nonmiss={:.3}, min_sam_nonmiss={:.3}, nsite_left = {:>6}, nsam_left={:>6}, overall_call_target_sites_samples={:.3}",
                min_site_nonmiss_rate,
                min_sample_nonmiss_rate,
                bcf_gt.selected_sites.len(),
                bcf_gt.selected_samples.len(),
                overall_nonmiss,
            );
        }

        Ok(bcf_gt.consolidate())
    }

    pub fn get_samples(&self) -> &[String] {
        self.sample_vec.as_slice()
    }

    pub fn into_genotype_siteinfo(
        self,
        valid_samples: &Samples,
        min_snp_sep: u32,
    ) -> Result<(Matrix<u8>, SiteInfoRaw)> {
        let BcfGenotype {
            chr_vec,
            pos_vec,
            gt_vec: dom_vec,
            sample_vec,
            chrname_map,
            selected_samples: _,
            selected_sites: _,
            nsam_targeted: nsam,
            args: _,
            ready_to_use,
        } = self;

        if !ready_to_use {
            return Err(Error::BcfGenotypeNotReady);
        }

        // rebuild the chromosome map and update chr_vec
        let mut v: Vec<(String, usize)> = chrname_map.into_iter().collect();
        v.sort_by_key(|(_chrname, chrid)| *chrid);
        let id_map: HashMap<usize, usize> = v
            .iter()
            .enumerate()
            .map(|(newid, oldid)| (oldid.1, newid))
            .collect();

        // get sites info
        let nsites = pos_vec.len();
        let mut chr_idx_vec = Vec::<usize>::with_capacity(nsites);
        let mut chr_pos_vec = Vec::<u32>::with_capacity(nsites);
        let mut last_chrid = -1i32;
        let mut last_pos = 0i32;

        // valid columns
        let valid_col = (0..nsam)
            .filter(|i| valid_samples.m().contains_key(&sample_vec[*i]))
            .collect_vec();

        // build genotype
        let n_valid_samples = valid_samples.v().len();
        let mut geno1 = MatrixBuilder::<u8>::new(n_valid_samples);

        // dbg!(nsam, chr_vec.len(), nsam * chr_vec.len(), dom_vec.len());
        for (i, (chr_idx, pos)) in chr_vec.into_iter().zip(pos_vec.into_iter()).enumerate() {
            if (last_chrid == chr_idx) && (pos - last_pos < min_snp_sep as i32) {
                continue;
            }
            chr_idx_vec.push(id_map[&(chr_idx as usize)]);
            chr_pos_vec.push(pos as u32);
            let s = i * nsam;
            let e = (i + 1) * nsam;
            let gt = &dom_vec[s..e];

            for x in gt
                .iter()
                .enumerate()
                .merge_join_by(valid_col.iter(), |a, b| a.0.cmp(*b))
            {
                if let EitherOrBoth::Both((_, allele), _) = x {
                    let allele = match allele {
                        -1 => None,
                        x => Some(*x as u8),
                    };
                    geno1.push(allele);
                }
            }
            last_chrid = chr_idx;
            last_pos = pos;
        }

        let geno1 = geno1.finish();

        let chrname_vec: Vec<String> = v.into_iter().map(|(chrname, _)| chrname).collect();
        let chrname_map: HashMap<String, usize> = chrname_vec
            .iter()
            .enumerate()
            .map(|(id, s)| (s.to_owned(), id))
            .collect();
        let siteinfo = SiteInfoRaw::from_parts(chrname_vec, chrname_map, chr_pos_vec, chr_idx_vec);

        Ok((geno1, siteinfo))
    }

    pub fn save_to_file(&self, output: &str) -> Result<()> {
        let writer = std::fs::File::create(output).map(BufWriter::new)?;
        bincode::serialize_into(writer, &self)?;
        Ok(())
    }

    pub fn load_from_file(input: &str) -> Result<Self> {
        let reader = std::fs::File::open(input).map(BufReader::new)?;
        Ok(bincode::deserialize_from(reader)?)
    }

    fn restricted_to_a_chromsome(&self, chrname: &str) -> Self {
        let chr_id = self.chrname_map[chrname] as i32;

        let chr_nsites = self
            .chr_vec
            .iter()
            .filter(|chrid| **chrid == chr_id)
            .count();
        let mut pos_vec = Vec::with_capacity(chr_nsites);
        let mut gt_vec = Vec::with_capacity(self.nsam_targeted * chr_nsites);
        let mut chr_vec = Vec::with_capacity(chr_nsites);
        let selected_sites: Vec<_> = (0..chr_nsites).collect();

        self.chr_vec
            .iter()
            .zip(self.pos_vec.iter())
            .zip(self.gt_vec.chunks(self.nsam_targeted))
            .filter(|((chrid, _pos), _gt_row)| **chrid == chr_id)
            .for_each(|((chrid, pos), gt_row)| {
                pos_vec.push(*pos);
                gt_vec.extend_from_slice(gt_row);
                chr_vec.push(*chrid);
            });
        Self {
            chr_vec,
            pos_vec,
            gt_vec,
            sample_vec: self.sample_vec.clone(),
            chrname_map: self.chrname_map.clone(),
            selected_samples: self.selected_samples.clone(),
            selected_sites,
            nsam_targeted: self.nsam_targeted,
            args: self.args.clone(),
            ready_to_use: true,
        }
    }

    pub fn split_chromosomes_into_files(&self, output_prefix: &str) -> Result<()> {
        if let Some(parent) = std::path::Path::new(output_prefix).parent() {
            std::fs::create_dir_all(parent)?;
        }
        for (chrname, chrid) in self.chrname_map.iter() {
            if self.chr_vec.iter().any(|id| *id as usize == *chrid) {
                let bcf_gt_chr = self.restricted_to_a_chromsome(&chrname);
                let p = format!("{output_prefix}_{chrname}.bin");
                bcf_gt_chr.save_to_file(&p)?;
            }
        }
        Ok(())
    }
}

// ---- helper function

/// Define a type alias 'BcfGzipReader' which is a specific configuration of the generic 'BcfReader'.
/// 'BcfReader' is instantiated with 'ParMultiGzipReader' that uses a dynamic trait object for reading from a boxed value.
type BcfGzipReader = BcfReader<ParMultiGzipReader<Box<dyn std::io::Read>>>;

/// get bcf reader from stdin or bcf file
fn get_bcf_gzip_reader(bcf_fname: &str) -> Result<BcfGzipReader> {
    let reader: Box<dyn std::io::Read> = if bcf_fname == "-" {
        Box::new(std::io::stdin().lock())
    } else {
        Box::new(std::fs::File::open(bcf_fname).map(BufReader::new)?)
    };
    let reader = BcfReader::from_reader(ParMultiGzipReader::from_reader(reader, 15, None, None)?);
    Ok(reader)
}

/// get indicators for target samples
fn get_targeted_sample_indicator(
    args: &BcfFilterArgs,
    vcf_samples: &[String],
) -> Result<Vec<bool>> {
    let mut is_sample_targeted = vec![];
    let nsam = vcf_samples.len();
    is_sample_targeted.resize(nsam, true);
    if let Some(p) = args.target_samples.as_ref() {
        // read target samples
        let targets: std::collections::HashSet<_> = std::fs::read_to_string(p)?
            .trim()
            .split('\n')
            .map(String::from)
            .collect();
        // check if each sample in vcf  is in target set
        vcf_samples.iter().enumerate().for_each(|(idx, sname)| {
            if !targets.contains(sname) {
                is_sample_targeted[idx] = false;
            }
        });
    }
    Ok(is_sample_targeted)
}

/// get maf of current site
fn get_maf_and_sorted_allele_count(
    gt: &[i8],
    nallele: usize,
    allele_counts: &mut Vec<(usize, u32)>,
) -> f32 {
    allele_counts.clear();
    allele_counts.resize(nallele, (0, 0));
    let mut tot_allele_counts = 0u32;
    for a in gt.iter() {
        let a = *a;
        if a < 0 {
            continue;
        } else {
            allele_counts[a as usize].1 += 1;
            tot_allele_counts += 1;
        }
    }
    allele_counts.sort_by_key(|(_idx, cnt)| u32::MAX - cnt); // reverse sort
    let maf = allele_counts[1].1 as f32 / tot_allele_counts as f32;
    maf
}

fn is_major_and_minor_allele_snps(
    record: &bcf_reader::Record,
    allele_counts: &[(usize, u32)],
) -> bool {
    let rngs = record.alleles();
    let rng = &rngs[allele_counts[0].0];
    let major_allele = &record.buf_shared()[rng.start..rng.end];
    let rng = &rngs[allele_counts[1].0];
    let minor_allele = &record.buf_shared()[rng.start..rng.end];
    if (major_allele.len() > 1)
        || (major_allele[0] == b'*')
        || (minor_allele.len() > 1)
        || (minor_allele[0] == b'*')
    {
        false
    } else {
        true
    }
}

/// test if current site has no filter or pass filter
fn is_variant_pass_in_filter_columns(
    header: &bcf_reader::Header,
    record: &bcf_reader::Record,
) -> Result<bool> {
    let mut pass = true;
    record
        .filters()
        .into_iter()
        .try_for_each(|nv_res| -> Result<()> {
            let fmt_idx = nv_res?
                .int_val()
                .ok_or(bcf_reader::Error::NumericaValueEmptyInt)?
                as usize;
            let filter_name = header.dict_strings()[&fmt_idx]["ID"].as_str();
            if filter_name != "PASS" {
                pass = false;
            }
            Ok(())
        })?;

    Ok(pass)
}
