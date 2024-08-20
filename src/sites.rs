use crate::args::RecombinationArg;
use crate::genome::GenomeInfo;
use crate::genome::{Genome, GenomeBuilder};
use crate::gmap::GeneticMap;
use std::collections::HashMap;
use std::ops::Range;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0:?}")]
    GenomeError(#[from] crate::genome::Error),
    #[error("{0:?}")]
    GmapError(#[from] crate::gmap::Error),
}

#[derive(Debug, Clone, Default)]
pub struct Sites {
    /// sorted genome-wide positions
    gw_pos_cm: Vec<f64>,
    gw_pos: Vec<u32>,
    // ranges of pos for chromosome
    chrom_ranges: Vec<Range<u32>>,
}

impl Sites {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn add(&mut self, gw_pos: u32, gw_pos_cm: f64) {
        self.gw_pos.push(gw_pos);
        self.gw_pos_cm.push(gw_pos_cm);
    }

    pub fn finish(&mut self, genome: &Genome) {
        // x>=boundary
        // x<bool

        // asert is sorted
        assert!(self
            .gw_pos
            .iter()
            .zip(self.gw_pos.iter().skip(1))
            .all(|(x, y)| x < y));
        for chrgwstart in genome.get_gwchrstarts() {
            let s = self.gw_pos.partition_point(|x| x < chrgwstart) as u32;
            if let Some(last) = self.chrom_ranges.last_mut() {
                last.end = s;
            }
            self.chrom_ranges.push(s..s);
        }
        if let Some(last) = self.chrom_ranges.last_mut() {
            last.end = self.gw_pos.len() as u32;
        }
    }

    pub fn has_gw_pos(&self, gw_pos: u32) -> bool {
        self.gw_pos.binary_search(&gw_pos).is_ok()
    }

    pub fn is_identical(&self, other: &Self) -> bool {
        (self.gw_pos == other.gw_pos) && (self.chrom_ranges == other.chrom_ranges)
    }

    pub fn get_pos_slice(&self) -> &[u32] {
        &self.gw_pos[..]
    }
    pub fn get_pos_cm_slice(&self) -> &[f64] {
        &self.gw_pos_cm[..]
    }
    pub fn get_chrom_pos_idx_ranges(&self, chrid: usize) -> (usize, usize) {
        let r = &self.chrom_ranges[chrid];
        (r.start as usize, r.end as usize)
    }
}

#[derive(PartialEq, Eq, Debug, Default)]
pub struct SiteInfoRaw {
    chrname_vec: Vec<String>,
    chrname_map: HashMap<String, usize>,
    chr_pos_vec: Vec<u32>,
    chr_idx_vec: Vec<usize>,
}
impl SiteInfoRaw {
    pub fn from_parts(
        chrname_vec: Vec<String>,
        chrname_map: HashMap<String, usize>,
        chr_pos_vec: Vec<u32>,
        chr_idx_vec: Vec<usize>,
    ) -> Self {
        Self {
            chrname_vec,
            chrname_map,
            chr_pos_vec,
            chr_idx_vec,
        }
    }
    pub fn new() -> Self {
        Self::default()
    }
    pub fn add_chr_name(&mut self, chrname: &str) {
        self.chrname_map
            .insert(chrname.to_owned(), self.chrname_vec.len());
        self.chrname_vec.push(chrname.to_owned());
    }
    pub fn add_chr_idx(&mut self, chr_name: &str) {
        let idx = self.chrname_map[chr_name];
        self.chr_idx_vec.push(idx);
    }
    pub fn add_chr_pos(&mut self, chr_pos: u32) {
        self.chr_pos_vec.push(chr_pos);
    }

    pub fn into_sites_and_genome(self, rec_args: &RecombinationArg) -> Result<(Sites, Genome)> {
        match rec_args.genome.as_ref() {
            Some(genome_toml_path) => {
                let ginfo =
                    GenomeInfo::from_toml_file(genome_toml_path).map_err(Error::GenomeError)?;
                let genome = Genome::from_genome_info(&ginfo);
                let gmap = GeneticMap::from_genome_info(&ginfo)?;
                let mut sites = Sites::new();

                // chrid map
                let idmap: HashMap<usize, usize> = self
                    .chrname_vec
                    .iter()
                    .enumerate()
                    .map(|(oldid, chrname)| (oldid, ginfo.idx[chrname]))
                    .collect();

                self.chr_idx_vec
                    .iter()
                    .zip(self.chr_pos_vec.iter())
                    .for_each(|(oldid, chr_pos)| {
                        let newid = idmap[oldid];
                        let gw_pos = ginfo.to_gw_pos(newid, *chr_pos);
                        let gw_pos_cm = gmap.get_cm(gw_pos);

                        sites.add(gw_pos, gw_pos_cm);
                    });
                sites.finish(&genome);
                Ok((sites, genome))
            }
            None => {
                let mut sites = Sites::new();
                let mut gbuilder = GenomeBuilder::new();
                let recom_rate = rec_args.rec_rate;
                self.chr_idx_vec
                    .iter()
                    .zip(self.chr_pos_vec.iter())
                    .try_for_each(|(oldid, chr_pos)| -> Result<()> {
                        let (_, gw_pos) =
                            gbuilder.encode_pos(&self.chrname_vec[*oldid], *chr_pos)?;
                        let gw_pos_cm = gw_pos as f64 * recom_rate * 100.0;
                        sites.add(gw_pos, gw_pos_cm);
                        Ok(())
                    })?;
                let genome = gbuilder.finish();
                sites.finish(&genome);

                Ok((sites, genome))
            }
        }
    }
}
