// use crate::container::intervals::Intervals;
use serde::Deserialize;
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

#[derive(Deserialize, Debug)]
struct GenomeFile {
    name: String,
    chromsize: Vec<u32>,
    idx: HashMap<String, usize>,
    chromnames: Vec<String>,
    gmaps: Vec<String>,
}

impl GenomeFile {
    fn check(&self) {
        assert_eq!(self.chromsize.len(), self.chromnames.len());
        assert_eq!(self.chromsize.len(), self.gmaps.len());
        assert!(self.chromnames.iter().all(|x| self.idx.contains_key(x)));
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
        Self {
            name: String::new(),
            chromsize: vec![],
            chromnames: Vec::new(),
            idx: HashMap::new(),
            gwstarts: vec![],
            gmaps: vec![],
        }
    }

    pub fn from_toml_file<P>(path: P) -> Self
    where
        P: AsRef<Path>,
    {
        let mut s = String::new();
        let p: &Path = path.as_ref();
        std::fs::File::open(p)
            .unwrap()
            .read_to_string(&mut s)
            .unwrap();
        let mut gfile: GenomeFile = toml::from_str(&s).unwrap();
        gfile.check();
        let mut ginfo = Self::new();
        use std::mem::swap;
        swap(&mut gfile.name, &mut ginfo.name);
        swap(&mut gfile.chromsize, &mut ginfo.chromsize);
        swap(&mut gfile.chromnames, &mut ginfo.chromnames);
        swap(&mut gfile.idx, &mut ginfo.idx);
        swap(&mut gfile.gmaps, &mut ginfo.gmaps);
        ginfo.gwstarts.push(0u32);
        for chrlen in &ginfo.chromsize[0..(ginfo.chromsize.len() - 1)] {
            let last = ginfo.gwstarts.last().unwrap();
            ginfo.gwstarts.push(*chrlen + *last);
        }
        ginfo
    }
    pub fn to_gw_pos(&self, chrid: usize, pos: u32) -> u32 {
        self.gwstarts[chrid] as u32 + pos
    }
    pub fn to_chr_pos(&self, gw_pos: u32) -> (usize, &str, u32) {
        let chrid = self.gwstarts.partition_point(|x| *x <= (gw_pos)) - 1;
        let pos = gw_pos - self.gwstarts[chrid] as u32;
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
