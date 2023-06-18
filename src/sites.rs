use crate::genome::Genome;
use std::ops::Range;
pub struct Sites {
    /// sorted genome-wide positions
    gw_pos: Vec<u32>,
    // ranges of pos for chromosome
    chrom_ranges: Vec<Range<u32>>,
}

impl Sites {
    pub fn new() -> Self {
        Self {
            gw_pos: vec![],
            chrom_ranges: vec![],
        }
    }
    pub fn add(&mut self, gw_pos: u32) {
        self.gw_pos.push(gw_pos);
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
            match self.chrom_ranges.last_mut() {
                Some(last) => last.end = s,
                None => {}
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
    pub fn get_chrom_pos_idx_ranges(&self, chrid: usize) -> (usize, usize) {
        let r = &self.chrom_ranges[chrid];
        (r.start as usize, r.end as usize)
    }
}
