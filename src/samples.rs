use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader, Read},
};

use crate::args::Arguments;
pub struct Samples {
    v: Vec<String>,
    m: std::collections::HashMap<String, u32>,
    s1: Vec<u32>,
    s2: Vec<u32>,
}

impl Samples {
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
        let mut s1 = vec![];
        let mut s2 = vec![];
        std::fs::File::open(&args.data_file1)
            .map(BufReader::new)
            .unwrap()
            .read_line(&mut s)
            .unwrap();
        for x in s.trim().split("\t").skip(2) {
            if !bad_samples.contains(x) {
                assert!(!m.contains_key(x), "duplicated sample names");
                m.insert(x.to_owned(), v.len() as u32);
                s1.push(v.len() as u32);

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
                    s2.push(v.len() as u32);
                    v.push(x.to_owned());
                }
            }
        }

        Self { m, v, s1, s2 }
    }

    pub fn v(&self) -> &Vec<String> {
        &self.v
    }
    pub fn m(&self) -> &HashMap<String, u32> {
        &self.m
    }
    pub fn s1(&self) -> &Vec<u32> {
        &self.s1
    }
    pub fn s2(&self) -> &Vec<u32> {
        &self.s2
    }
}

#[test]
fn read_samples() {
    let args = Arguments::new_for_test();
    Samples::from_args(&args);
}
