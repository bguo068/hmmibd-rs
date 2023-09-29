use clap::{Args, Parser};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about, name = "hmmibd2")]
pub struct Arguments {
    /// File of genotype data.
    /// Format for genotype file: tab-delimited text file, with one single
    /// nucleotide polymorphism (SNP) per line. The first two columns are the
    /// chromosome and position, followed by one sample per column. A header
    /// line, giving the sample names, is required. Genotypes are coded by
    /// number: -1 for missing data, 0 for the first allele, 1 for the second,
    /// etc. SNPs and indels (if you trust them) can thus be treated on an equal
    /// footing. The variants must be in chromosome and position order, and can
    /// have between two and eight alleles (more, if you feel like changing
    /// max_allele).
    #[arg(short = 'i', long, required = true)]
    pub data_file1: PathBuf,

    /// Optional: file of genotype data from a second population
    #[arg(short = 'I', long)]
    pub data_file2: Option<PathBuf>,

    ///  Optional: File of allele frequencies for the sample population. Format:
    ///  tab-delimited, no header, one variant per row. Line format: <chromosome
    ///  (int)> <position (bp, int)> <allele 1 freq> <all 2 freq> [<all 3 freq>]
    ///  ... The genotype and frequency files must contain exactly the same
    ///  variants, in the same order. If no file is supplied, allele frequencies
    ///  are calculated from the input data file.
    #[arg(short = 'f', long)]
    pub freq_file1: Option<PathBuf>,

    ///  Optional: File of allele frequencies for the second population; same format as
    ///  for -f
    #[arg(short = 'F', long)]
    pub freq_file2: Option<PathBuf>,

    /// Optional: Maximum number of fit iterations
    #[arg(short = 'm', long, default_value = "5")]
    pub max_iter: u32,

    /// Optional: file of sample ids to exclude from all analysis. Format: no header, one
    /// id (string) per row. Note: b stands for "bad samples"
    #[arg(short = 'b', long)]
    pub bad_file: Option<PathBuf>,

    ///  Optional: File of sample pairs to analyze; all others are not processed by the HMM
    ///  (but are still used to calculate allele frequencies). Format: no header,
    ///  tab-delimited, two sample ids (strings) per row. Note: "g" stands for
    ///  "good pairs".
    #[arg(short = 'g', long)]
    pub good_file: Option<PathBuf>,

    ///  Optional: Cap on the number of generations (floating point). Sets the maximum
    ///  value for that parameter in the fit. This is useful if you are
    ///  interested in recent IBD and are working with a population with
    ///  substantial linkage disequilbrium. Specifying a small value will force
    ///  the program to assume little recombination and thus a low transition
    ///  rate; otherwise it will identify the small blocks of LD as ancient IBD,
    ///  and will force the number of generations to be large.
    #[arg(short = 'n', long, default_value = "Inf")]
    pub k_rec_max: f64,

    /// output put prefix, if not specified use the prefix as `-i` option.
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,

    /// error rate in genotype calls
    #[arg(long, default_value_t = 0.001)]
    pub eps: f64,

    /// minimum number of informative sites in a pairwise comparison (those
    /// with minor allele)
    #[arg(long, default_value_t = 10)]
    pub min_inform: u32,

    /// minimum discordance rate in comparison; set > 0 to skip identical pairs
    #[arg(long, default_value_t = 0.0)]
    pub min_discord: f64,

    /// maximum discordance rate in comparison; set < 1 to skip unrelated pairs
    #[arg(long, default_value_t = 1.0)]
    pub max_discord: f64,

    /// skip next snp(s) if too close to last one; in bp
    #[arg(long, default_value_t = 5)]
    pub min_snp_sep: u32,

    /// covergence criteria: min delta pi
    #[arg(long, default_value_t = 0.001)]
    pub fit_thresh_dpi: f64,

    /// covergence criteria: min delta k_rec
    #[arg(long, default_value_t = 0.01)]
    pub fit_thresh_dk: f64,

    /// covergence criteria: min (delta k_rec) / k_rec
    #[arg(long, default_value_t = 0.001)]
    pub fit_thresh_drelk: f64,

    /// max number of unique alleles per site
    #[arg(long, default_value_t = 8)]
    pub max_all: u32,

    #[command(flatten)]
    pub rec_args: RecombinationArg,

    /// filtering IBD segments: if set, none of IBD/DBD segments short than
    /// filt_min_seg_cm will not be written to hmm.txt files
    #[arg(long)]
    pub filt_min_seg_cm: Option<f64>,
    /// filtering IBD segments: if set, none of IBD/DBD segments from pairs with
    /// k_rec > filt_max_tmrca will not be written to hmm.txt files
    #[arg(long)]
    pub filt_max_tmrca: Option<f64>,
    /// filtering IBD segments: if set, no DBD (non-IBD) segments will not be written to hmm.txt files
    #[arg(long)]
    pub filt_ibd_only: bool,

    /// number of threads. 0 : use all cpus; non-zero: use the given numbers of threads
    #[arg(long, default_value_t = 0)]
    pub num_threads: usize,

    /// number of pairs of samples per chunk for parallelization
    #[arg(long, default_value_t = 120)]
    pub par_chunk_size: usize,
}

impl Arguments {
    pub fn new_for_test() -> Self {
        Self {
            data_file1: PathBuf::from("samp_data/pf3k_Cambodia_13.txt"),
            data_file2: Some(PathBuf::from("samp_data/pf3k_Ghana_13.txt")),
            freq_file1: Some(PathBuf::from("samp_data/freqs_pf3k_Cambodia_13.txt")),
            freq_file2: Some(PathBuf::from("samp_data/freqs_pf3k_Ghana_13.txt")),
            max_iter: 5,
            bad_file: None,
            good_file: None,
            k_rec_max: f64::MAX,
            output: None,
            eps: 0.001,
            min_inform: 10,
            min_discord: 0.0,
            max_discord: 1.0,
            min_snp_sep: 5,
            rec_args: RecombinationArg {
                rec_rate: 7.4e-7,
                genome: None,
            },
            fit_thresh_dpi: 0.001,
            fit_thresh_dk: 0.01,
            fit_thresh_drelk: 0.001,
            max_all: 8,
            filt_min_seg_cm: None,
            filt_max_tmrca: None,
            filt_ibd_only: false,
            num_threads: 1,
            par_chunk_size: 120,
        }
    }
}

#[derive(Args, Debug)]
#[group(required = false, multiple = false)]
pub struct RecombinationArg {
    /// recombination rate per generation per basepair. When used --genome should not be specified.
    #[arg(short = 'r', long, default_value_t = 7.4e-7)]
    pub rec_rate: f64,
    /// genome file which specifies chromosome size, names, and plink genetic map file paths.
    ///  When used --rec-rate should not be specified.
    #[arg(long)]
    pub genome: Option<PathBuf>,
}
