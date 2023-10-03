use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader, Read},
};

use crate::args::Arguments;

#[derive(Clone)]
pub struct Samples {
    v: Vec<String>,
    m: std::collections::HashMap<String, u32>,
    pop1_nsam: u32,
    pop2_nsam: u32,
}

impl Samples {
    pub fn from_slice(pop1_slice: &[String], pop2_slice: &[String]) -> Self {
        let pop1_nsam = pop1_slice.len() as u32;
        let pop2_nsam = pop2_slice.len() as u32;
        let total = pop1_nsam + pop2_nsam;
        let mut v = Vec::<String>::with_capacity(total as usize);
        let mut m = std::collections::HashMap::<String, u32>::with_capacity(total as usize);
        v.extend_from_slice(pop1_slice);
        v.extend_from_slice(pop2_slice);
        for (i, s) in pop1_slice.iter().enumerate() {
            m.insert(s.to_owned(), i as u32);
        }
        for (i, s) in pop2_slice.iter().enumerate() {
            m.insert(s.to_owned(), i as u32 + pop1_nsam);
        }

        Self {
            v,
            m,
            pop1_nsam,
            pop2_nsam,
        }
    }
    pub fn from_args(args: &Arguments) -> Self {
        let mut s = String::new();
        // get bad samples
        let mut bad_samples = HashSet::<String>::new();
        if let Some(bad_file) = args.bad_file.as_ref() {
            std::fs::File::open(bad_file)
                .map(BufReader::new)
                .unwrap()
                .read_to_string(&mut s)
                .unwrap();
            for x in s.trim().split("\t") {
                bad_samples.insert(x.to_owned());
            }
        }

        // get samples from data headers
        s.clear();
        let mut v = vec![];
        let mut m = HashMap::new();
        let mut pop1_nsam = 0u32;
        let mut pop2_nsam = 0u32;
        std::fs::File::open(&args.data_file1)
            .map(BufReader::new)
            .unwrap()
            .read_line(&mut s)
            .unwrap();
        for x in s.trim().split("\t").skip(2) {
            if !bad_samples.contains(x) {
                assert!(!m.contains_key(x), "duplicated sample names");
                m.insert(x.to_owned(), v.len() as u32);
                pop1_nsam += 1;

                v.push(x.to_owned());
            }
        }

        // data files
        s.clear();
        if let Some(data_file2) = args.data_file2.as_ref() {
            std::fs::File::open(data_file2)
                .map(BufReader::new)
                .unwrap()
                .read_line(&mut s)
                .unwrap();
            for x in s.trim().split("\t").skip(2) {
                if !bad_samples.contains(x) {
                    assert!(!m.contains_key(x), "duplicated sample names");
                    m.insert(x.to_owned(), v.len() as u32);
                    pop2_nsam += 1;
                    v.push(x.to_owned());
                }
            }
        }

        assert_eq!(pop1_nsam + pop2_nsam, v.len() as u32);

        Self {
            m,
            v,
            pop1_nsam,
            pop2_nsam,
        }
    }

    pub fn v(&self) -> &Vec<String> {
        &self.v
    }
    pub fn m(&self) -> &HashMap<String, u32> {
        &self.m
    }
    pub fn pop1_nsam(&self) -> u32 {
        self.pop1_nsam
    }
    pub fn pop2_nsam(&self) -> u32 {
        self.pop2_nsam
    }
}

#[test]
fn read_samples() {
    let args = Arguments::new_for_test();
    Samples::from_args(&args);
}
