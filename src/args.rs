use clap::{Args, Parser, ValueEnum};

#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about, name = "hmmibd-rs", color=clap::ColorChoice::Always, styles=get_styles())]
pub struct Arguments {
    /// File of genotype data.
    /// (1) by default, format for genotype file: tab-delimited text file, with
    /// one single nucleotide polymorphism (SNP) per line. The first two columns
    /// are the chromosome and position, followed by one sample per column.
    /// A header line, giving the sample names, is required. Genotypes are
    /// coded by number: -1 for missing data, 0 for the first allele, 1 for the
    /// second, etc. SNPs and indels (if you trust them) can thus be treated on
    /// an equal footing. The variants must be in chromosome and position order,
    /// and can have between two and eight alleles (more, if you feel like
    /// changing max-allele).
    /// (2) When `--from-bcf` is specified, the input is expected in BCF
    /// format, either from a file or stdin. If the argument is a file path,
    /// the BCF genotype is read from that file. If the argument is "-", the BCF
    /// genotype is read from stdin.
    /// (3) When `--from-bin` is specified, a binary genotype input file
    /// is expected. See options `--bcf-to-bin-file-by-chromosome` or
    /// `bcf-to-bin-file` for generating binary genotype files.
    #[arg(short = 'i', long, required = true, help_heading = "input data")]
    pub data_file1: String,

    /// Optional: file of genotype data from a second population
    #[arg(
        short = 'I',
        long,
        group = "grp_data_file2",
        help_heading = "input data"
    )]
    pub data_file2: Option<String>,

    ///  Optional: File of allele frequencies for the sample population. Format:
    ///  tab-delimited, no header, one variant per row. Line format: <chromosome
    ///  (int)> <position (bp, int)> <allele 1 freq> <all 2 freq> <all 3 freq>
    ///  ... The genotype and frequency files must contain exactly the same
    ///  variants, in the same order. If no file is supplied, allele frequencies
    ///  are calculated from the input data file.
    #[arg(short = 'f', long, help_heading = "input data")]
    pub freq_file1: Option<String>,

    ///  Optional: File of allele frequencies for the second population; same format as
    ///  for -f
    #[arg(
        short = 'F',
        long,
        requires = "grp_data_file2",
        help_heading = "input data"
    )]
    pub freq_file2: Option<String>,

    /// Optional: file of sample ids to exclude from all analysis. Format: no header, one
    /// id (string) per row. Note: b stands for "bad samples"
    #[arg(short = 'b', long, help_heading = "input data")]
    pub bad_file: Option<String>,

    ///  Optional: File of sample pairs to analyze; all others are not processed by the HMM
    ///  (but are still used to calculate allele frequencies). Format: no header,
    ///  tab-delimited, two sample ids (strings) per row. Note: "g" stands for
    ///  "good pairs".
    #[arg(short = 'g', long, help_heading = "input data")]
    pub good_file: Option<String>,

    // ---- bcf file
    /// Optional: flag indicating whether the input file is of BCF format
    #[arg(
        long,
        default_value_t = false,
        group = "input_format",
        group = "grp_from_bcf",
        help_heading = "input data options"
    )]
    pub from_bcf: bool,

    /// Optional: flag indicating whether the input file is binary genotype data
    #[arg(
        long,
        default_value_t = false,
        group = "input_format",
        help_heading = "input data options"
    )]
    pub from_bin: bool,

    /// Optional: Read mode for the Bcf input. Currently "first-ploidy" and
    /// "each-ploidy" are experimental.
    #[arg(
        long,
        value_enum,
        default_value = "dominant-allele",
        requires("grp_from_bcf"),
        help_heading = "input data options"
    )]
    pub bcf_read_mode: BcfReadMode,

    /// Optional: When the --from_bcf flag is set and the --bcf-read-mode
    /// dominant-allele option (default setting) is used, this option
    /// allows specifying a custom BCF filtering configuration file. If not
    /// provided, a built-in filtering criteria will be used, and a temporary
    /// configuration file consistent with the built-in criteria named
    /// "tmp_bcf_filter_config.toml" will be generated in the current folder.
    #[arg(long, requires = "grp_from_bcf", help_heading = "input data options")]
    pub bcf_filter_config: Option<String>,

    /// Optional: output prefix. If not specified, the prefix from the `-i`
    /// option will be used. Note: When reading genotype from stdin, this option
    /// must be specified.
    #[arg(short = 'o', long, help_heading = "output data")]
    pub output: Option<String>,

    /// Optional: whether to suppress frac output. Toggle this on for minimal IO
    /// burden when only IBD segment output is needed
    #[arg(long, default_value_t = false, help_heading = "output options")]
    pub suppress_frac: bool,

    /// Optional, valid only with the `--from-bcf` option. When enabled,
    /// genotypes are written to a binary file for each chromosome. Chromosomes
    /// with no sites are ignored. When this is set, the HMM IBD inference step
    /// is skipped. This option is designed to separate the BCF processing step
    /// from the HMM inference step, especially when genotypes of different
    /// chromosomes are concatenated for sample filtering, and the user wants to
    /// split the filtering genotype into chromosomes and run each instance of
    /// hmmibd-rs on a chromosome in parallel.
    #[arg(long, requires = "grp_from_bcf", help_heading = "output options")]
    pub bcf_to_bin_file_by_chromosome: bool,

    /// Optional, valid only with the `--from-bcf` option. Similar to the
    /// `--bcf-to-bin-file-by-chromosome` option, but generates a single binary
    /// file containing all chromosomes.
    #[arg(long, requires = "grp_from_bcf", help_heading = "output options")]
    pub bcf_to_bin_file: bool,

    // ---- hmm options
    /// Optional: Maximum number of fit iterations
    #[arg(short = 'm', long, default_value = "5", help_heading = "hmm options")]
    pub max_iter: u32,

    ///  Optional: Cap on the number of generations (floating point). Sets the maximum
    ///  value for that parameter in the fit. This is useful if you are
    ///  interested in recent IBD and are working with a population with
    ///  substantial linkage disequilbrium. Specifying a small value will force
    ///  the program to assume little recombination and thus a low transition
    ///  rate; otherwise it will identify the small blocks of LD as ancient IBD,
    ///  and will force the number of generations to be large.
    #[arg(short = 'n', long, default_value = "Inf", help_heading = "hmm options")]
    pub k_rec_max: f64,

    /// Optional: error rate in genotype calls
    #[arg(long, default_value_t = 0.001, help_heading = "hmm options")]
    pub eps: f64,

    /// Optional: minimum number of informative sites in a pairwise comparison
    /// (those with minor allele)
    #[arg(long, default_value_t = 10, help_heading = "hmm options")]
    pub min_inform: u32,

    /// Optional: minimum discordance rate in comparison from 0.0 to 1.0;  set >
    /// 0 to skip identical pairs
    #[arg(long, default_value_t = 0.0, help_heading = "hmm options")]
    pub min_discord: f64,

    /// Optional: maximum discordance rate in comparison from 0.0 to 1.0;  set <
    /// 1 to skip unrelated pairs
    #[arg(long, default_value_t = 1.0, help_heading = "hmm options")]
    pub max_discord: f64,

    /// Optional: skip next snp(s) if too close to last one; in bp
    #[arg(long, default_value_t = 5, help_heading = "hmm options")]
    pub min_snp_sep: u32,

    /// Optional: covergence criteria: min delta pi
    #[arg(long, default_value_t = 0.001, help_heading = "hmm options")]
    pub fit_thresh_dpi: f64,

    /// Optional: covergence criteria: min delta k_rec
    #[arg(long, default_value_t = 0.01, help_heading = "hmm options")]
    pub fit_thresh_dk: f64,

    /// Optional: covergence criteria: min (delta k_rec) / k_rec
    #[arg(long, default_value_t = 0.001, help_heading = "hmm options")]
    pub fit_thresh_drelk: f64,

    /// Optional: Maximum number of unique alleles per site. For instance, if
    /// the input contains only biallelic variants, setting `--max-all 2` could
    /// reduce memory usage and enhance computation speed by utilizing a more
    /// compact memory layout.
    #[arg(long, default_value_t = 8, help_heading = "memory options")]
    pub max_all: u32,

    // Optional: recombination rate or rate map
    #[command(flatten)]
    pub rec_args: RecombinationArg,

    // ---- memory buffer
    /// Optional: output buffer size in bytes for IBD segments, by default it is
    /// 8Kb. When IO is slow or not working well for writing many small chunks
    /// of data, it can be beneficial to set this to a larger value to reduce of
    /// the number  of times to call system IO operation.
    #[arg(long, help_heading = "memory options")]
    pub buffer_size_segments: Option<usize>,

    /// Optional: output buffer size in bytes for IBD fraction records,  by
    /// default it is 8Kb. used similarly to --buffer-size-segments option
    #[arg(long, help_heading = "memory options")]
    pub buffer_size_frac: Option<usize>,

    // -- ibdseg filtering
    /// Optional: filtering IBD segments. If set, IBD/non-IBD segments short
    /// than FILT-MIN-SEG-CM cM will not be written to hmm.txt files
    #[arg(long, help_heading = "output segment filtering options")]
    pub filt_min_seg_cm: Option<f64>,

    /// Optional: filtering IBD segments. If set, IBD/non-IBD segments from
    /// pairs with k_rec > filt_max_tmrca will not be written to hmm.txt files
    #[arg(long, help_heading = "output segment filtering options")]
    pub filt_max_tmrca: Option<f64>,

    /// Optional: filtering IBD segments. If set, non-IBD segments will not be
    /// written to hmm.txt files
    #[arg(long, help_heading = "output segment filtering options")]
    pub filt_ibd_only: bool,

    // ---- parallelization
    /// Optional: number of threads. "0" : use all cpus; non-zero: use the given
    /// numbers of threads
    #[arg(long, default_value_t = 0, help_heading = "parallelization options")]
    pub num_threads: usize,

    /// Optional: parallelization mode. Mode "0", chunks of sample pairs: creates a
    /// vector of sample pairs, then splits it into equal lengths (specified in
    /// --par-chunk-size); different threads process these chunks in parallel.
    /// This model is intended for smaller sample sizes. Mode "1", pairs of sample
    /// chunks: creates a vector of samples, splits it into equal lengths, then
    /// creates a second vector of pairs of sample chunks; different threads
    /// process these pairs of sample chunks in parallel. As samples within each
    /// chunk are located close to each other but those from different chunks
    /// may be far apart, this mode first copies each pair of chunks to the
    /// thread heap so that all samples within these pairs of chunks are close
    /// to each other, aiming to better utilize memory cache locality. This mode
    /// is intended for working with larger sample sizes.
    #[arg(long, default_value_t = 0, help_heading = "parallelization options")]
    pub par_mode: u8,

    /// Optional: number of sample pairs per chunk for parallelization mode 0, or
    /// number of samples per chunk for parallelization mode 1
    #[arg(long, default_value_t = 120, help_heading = "parallelization options")]
    pub par_chunk_size: u32,
}

#[derive(Default, Debug, ValueEnum, Clone, Copy)]
pub enum BcfReadMode {
    #[default]
    DominantAllele,
    FirstPloidy,
    EachPloidy,
}

impl Arguments {
    pub fn new_for_test() -> Self {
        Self {
            data_file1: String::from("c/samp_data/pf3k_Cambodia_13.txt"),
            data_file2: Some(String::from("c/samp_data/pf3k_Ghana_13.txt")),
            from_bcf: false,
            from_bin: false,
            bcf_filter_config: None,
            freq_file1: Some(String::from("c/samp_data/freqs_pf3k_Cambodia_13.txt")),
            freq_file2: Some(String::from("c/samp_data/freqs_pf3k_Ghana_13.txt")),
            max_iter: 5,
            bad_file: None,
            good_file: None,
            k_rec_max: f64::MAX,
            output: Some("tmp_hmmibdrs".into()),
            suppress_frac: false,
            buffer_size_segments: Some(1000000),
            buffer_size_frac: Some(10000),
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
            bcf_read_mode: Default::default(),
            par_mode: 0,
            bcf_to_bin_file_by_chromosome: false,
            bcf_to_bin_file: false,
        }
    }
    pub fn new_for_test_bcf() -> Self {
        Self {
            data_file1: String::from("testdata/pf7_data/pf7_chr1_20samples_maf0.01.bcf"),
            data_file2: None,
            from_bcf: true,
            from_bin: false,
            bcf_filter_config: Some(String::from("testdata/pf7_data/dom_gt_config.toml")),
            freq_file1: None,
            freq_file2: None,
            max_iter: 5,
            bad_file: None,
            good_file: None,
            k_rec_max: f64::MAX,
            output: None,
            suppress_frac: false,
            buffer_size_segments: Some(1000000),
            buffer_size_frac: Some(10000),
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
            par_mode: 0,
            bcf_read_mode: Default::default(),
            bcf_to_bin_file_by_chromosome: false,
            bcf_to_bin_file: false,
        }
    }
}

#[derive(Args, Debug, Clone)]
#[group(required = false, multiple = false)]
pub struct RecombinationArg {
    /// Optional: constant recombination rate per generation per base pair.
    /// When used, --genome should not be specified.
    #[arg(
        short = 'r',
        long,
        default_value_t = 7.4e-7,
        group = "recombination_options",
        help_heading = "hmm options"
    )]
    pub rec_rate: f64,
    /// Optional: Genome file specifying chromosome sizes, names, and paths to
    /// PLINK genetic map files. When used, do not specify --rec-rate.
    #[arg(long, group = "recombination_options", help_heading = "hmm options")]
    pub genome: Option<String>,
}

pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}
