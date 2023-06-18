use std::collections::HashMap;

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
}

pub struct GenomeBuilder {
    // chromsome and its index
    gwchrstarts: Vec<u32>,
    chromnames: Vec<String>,
    idx: HashMap<String, (u32, u32)>,
}

impl GenomeBuilder {
    pub fn new() -> Self {
        Self {
            gwchrstarts: vec![],
            chromnames: vec![],
            idx: HashMap::new(),
        }
    }

    /// convert string to idx and inform GENOME about the chromosme length
    pub fn encode_pos(&mut self, chrname: &str, pos: u32) -> (u32, u32) {
        if self.idx.contains_key(chrname) {
            let (chrid, max_pos) = self.idx[chrname];
            assert!(max_pos < pos);
            *self.idx.get_mut(chrname).unwrap() = (chrid, pos);
            let gw_chr_start = self.gwchrstarts.last().unwrap();
            (chrid, pos + gw_chr_start)
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

            (chrid, pos + this_gw_start)
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
