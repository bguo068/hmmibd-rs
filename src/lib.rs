// genotype table per chromosome
// - from file
// - remove bad samples
// - generate iterator of concordance or return a vector

// ProgArguments

// site table per chromosome
// - distance (cM or bp)
// - allele freuqency table
// - from file
// - from genotype table

// model
// - rate
// - k

// HmmRuner
// - pair
// - model
// - FBRunner
// - ViterbiRunner
// - BWRunner

// DataHandler
// - parameters
//     - esp
// - Genotype matrix
// - samples
// - pairs
// - allele frequency table
// - get_pair_data

// - PairData
//     - Discordance buffer
//     - allele frequency buffer
//     - nalt
//     - calculate_a (model, s1, s2)
//     - calculate_b (model, s1, s2, pos)
//     - gw_variable
//     - chr_variable
//        - finisht_fit
//        -

pub mod args;
pub mod bcf;
pub mod data;
pub mod genome;
pub mod genome2;
pub mod gmap;
pub mod hmm;
pub mod matrix;
pub mod model;
pub mod samples;
pub mod sites;
