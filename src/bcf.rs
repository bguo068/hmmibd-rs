use bcf_reader::*;
use itertools::{EitherOrBoth, Itertools};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, io::BufReader};

use crate::{
    matrix::{Matrix, MatrixBuilder},
    samples::Samples,
    sites::SiteInfoRaw,
};

#[derive(Default, Serialize, Deserialize, Debug, Clone)]
pub struct DominantGenotypeArgs {
    pub min_depth: u32,
    pub min_ratio: f32,
    pub min_maf: f32,
    pub min_r1_r2: f32,
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

    #[error("Dominant genotype is empty: {file}:{line}")]
    DominantGenotypeEmpty { file: &'static str, line: u32 },
}

impl DominantGenotypeArgs {
    pub fn new_from_toml_file(dom_gt_config_path: &str) -> Result<Self> {
        let toml_str = std::fs::read_to_string(dom_gt_config_path)?;
        let dga: DominantGenotypeArgs = toml::from_str(&toml_str)?;
        Ok(dga)
    }
    pub fn new_from_builtin() -> Result<Self> {
        let toml_str = r#"
min_depth = 5
min_ratio = 0.7
min_r1_r2 = 3.0
min_maf = 0.01
min_site_nonmissing = 0.1
nonmissing_rates= [
    [0.05, 0.025],
    [0.11, 0.55 ],
    [0.15, 0.075],
    [0.21, 0.105],
    [0.25, 0.125],
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

        let dga: DominantGenotypeArgs = toml::from_str(toml_str)?;
        Ok(dga)
    }
}

#[derive(Default, Deserialize, Debug)]
pub struct DominantGenotype {
    chr_vec: Vec<i32>,
    pos_vec: Vec<i32>,
    dom_vec: Vec<i8>,
    sample_vec: Vec<String>,
    chrname_map: HashMap<String, usize>,
    selected_samples: Vec<usize>,
    selected_sites: Vec<usize>,
    nsam_targeted: usize,
    args: DominantGenotypeArgs,
}

impl DominantGenotype {
    fn new(bcffilter_args: &DominantGenotypeArgs) -> Self {
        Self {
            chr_vec: vec![],
            pos_vec: vec![],
            dom_vec: vec![],
            sample_vec: vec![],
            chrname_map: HashMap::new(),
            selected_samples: vec![],
            selected_sites: vec![],
            nsam_targeted: 0,
            args: bcffilter_args.clone(),
        }
    }
    fn read_dom(&mut self, bcf_fname: impl AsRef<std::path::Path>) -> Result<Header> {
        let reader = std::fs::File::open(bcf_fname.as_ref()).map(BufReader::new)?;
        let mut reader =
            BcfReader::from_reader(ParMultiGzipReader::from_reader(reader, 3, None, None)?);
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

        // get indicators for target samples
        let mut is_sample_targeted = vec![];
        is_sample_targeted.resize(nsam, true);
        if let Some(p) = self.args.target_samples.as_ref() {
            // read target samples
            let targets: std::collections::HashSet<_> = std::fs::read_to_string(p)?
                .trim()
                .split('\n')
                .map(String::from)
                .collect();
            // check if each sample in vcf  is in target set
            header
                .get_samples()
                .iter()
                .enumerate()
                .for_each(|(idx, sname)| {
                    if !targets.contains(sname) {
                        is_sample_targeted[idx] = false;
                    }
                });
        }
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

            // not segrating
            if nallele <= 1 {
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
                self.dom_vec.push(dom_allele);
            }
            let non_missing_rate = site_nonmiss_counter as f32 / nsam as f32;

            // calculate minor allele frequency based on dom_allele for this site
            allele_counts.clear();
            allele_counts.resize(nallele, (0, 0));
            let mut tot_allele_counts = 0u32;
            for a in self.dom_vec[(self.dom_vec.len() - nsam_in_targets)..].iter() {
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

            // at the end of the line, check if the site has too many missing dominant allele
            // if yes remove those
            if (non_missing_rate < self.args.min_site_nonmissing) || (maf < self.args.min_maf) {
                // remove alleles for this site
                self.dom_vec.resize(self.dom_vec.len() - nsam_in_targets, 0);
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
                if self.dom_vec[*row * self.nsam_targeted + col] != -1 {
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
                if self.dom_vec[*row * self.nsam_targeted + col] != -1 {
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
                if self.dom_vec[*row * self.nsam_targeted + col] != -1 {
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

        new_dg.dom_vec.reserve(nsites * nsam);
        for site in self.selected_sites.iter() {
            for sample in self.selected_samples.iter() {
                new_dg
                    .dom_vec
                    .push(self.dom_vec[*site * self.nsam_targeted + *sample]);
            }
        }

        new_dg.chrname_map = self.chrname_map.clone();

        new_dg
    }

    pub fn new_from_processing_bcf(
        dgt_args: &DominantGenotypeArgs,
        bcf_path: impl AsRef<std::path::Path>,
    ) -> Result<Self> {
        let mut dg = Self::new(dgt_args);
        dg.read_dom(bcf_path)?;

        println!(
            "before filtering: n good site = {}, n good samples = {}",
            dg.selected_sites.len(),
            dg.selected_samples.len()
        );

        let rates = dg.args.nonmissing_rates.clone();
        for (min_site_nonmiss_rate, min_sample_nonmiss_rate) in rates {
            dg.filter_by_missingness(min_site_nonmiss_rate, min_sample_nonmiss_rate);
            let overall_nonmiss = dg.calc_nonmiss();
            println!(
                "min_site_nonmiss={:.3}, min_sam_nonmiss={:.3}, nsite_left = {:>6}, nsam_left={:>6}, overall_call_target_sites_samples={:.3}",
                min_site_nonmiss_rate,
                min_sample_nonmiss_rate,
                dg.selected_sites.len(),
                dg.selected_samples.len(),
                overall_nonmiss,
            );
        }

        Ok(dg.consolidate())
    }

    pub fn get_samples(&self) -> &[String] {
        self.sample_vec.as_slice()
    }

    pub fn into_genotype_siteinfo(
        self,
        valid_samples: &Samples,
        min_snp_sep: u32,
    ) -> (Matrix<u8>, SiteInfoRaw) {
        let DominantGenotype {
            chr_vec,
            pos_vec,
            dom_vec,
            sample_vec,
            chrname_map,
            selected_samples: _,
            selected_sites: _,
            nsam_targeted: nsam,
            args: _,
        } = self;

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

        (geno1, siteinfo)
    }
}

// pub fn get_dominant_geotype(
//     dgt_args: &DominantGenotypeArgs,
//     bcf_path: impl AsRef<std::path::Path>,
// ) -> DominantGenotype {
//     let mut dg = DominantGenotype::new(dgt_args);
//     dg.read_dom(bcf_path);

//     println!(
//         "before filteing: n good site = {}, n good samples = {}",
//         dg.selected_sites.len(),
//         dg.selected_samples.len()
//     );

//     let rates = dg.args.nonmissing_rates.clone();
//     for (min_site_nonmiss_rate, min_sample_nonmiss_rate) in rates {
//         dg.filter_by_missingness(min_site_nonmiss_rate, min_sample_nonmiss_rate);
//         let overall_nonmiss = dg.calc_nonmiss();
//         println!(
//             "min_site_nonmiss={:.3}, min_sam_nonmiss={:.3}, nsite_left = {:>6}, nsam_left={:>6}, overall_call_target_sites_samples={:.3}",
//             min_site_nonmiss_rate,
//             min_sample_nonmiss_rate,
//             dg.selected_sites.len(),
//             dg.selected_samples.len(),
//             overall_nonmiss,
//         );
//     }

//     dg.consolidate()
// }
