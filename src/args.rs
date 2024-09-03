use clap::{Args, Parser, ValueEnum};

#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about, name = "hmmibd-rs", color=clap::ColorChoice::Always, styles=get_styles())]
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
        help_heading = "input data option"
    )]
    pub from_bcf: bool,

    /// Optional: flag indicating whether the input file is binary genotype data
    #[arg(
        long,
        default_value_t = false,
        group = "input_format",
        help_heading = "input data option"
    )]
    pub from_bin: bool,

    /// Optional: Read mode for the Bcf input
    #[arg(
        long,
        value_enum,
        default_value = "dominant-allele",
        requires("grp_from_bcf"),
        help_heading = "input data option"
    )]
    pub bcf_read_mode: BcfReadMode,

    /// Optional: when --from_bcf flag and --bcf-read-mode dominant-allele option are used, use this to specify bcf filtering config file. If not provided, a builtin filtering criteria is used.
    #[arg(long, requires = "grp_from_bcf", help_heading = "input data option")]
    pub bcf_filter_config: Option<String>,

    /// output put prefix, if not specified use the prefix as `-i` option.
    #[arg(short = 'o', long, help_heading = "output data")]
    pub output: Option<String>,

    /// whether suppress frac output. Toggle this on for minimal IO burden when
    /// only IBD segment output is needed
    #[arg(long, default_value_t = false, help_heading = "output option")]
    pub suppress_frac: bool,

    /// Applied with the `--from-bcf` option. When enabled, genotypes are
    /// written to a binary file by chromosome. Chromosomes with no sites are
    /// ignored. When this is set, hmmibd inference is skipped.
    #[arg(long, requires = "grp_from_bcf", help_heading = "output option")]
    pub bcf_to_bin_file_by_chromosome: bool,

    /// Used with the --from-bcf option. When set, the genotype data for all
    /// chromosomes will be written to a single bin file. Chromosomes with zero
    /// sites will be ignored. When this is set, hmmibd will be skipped.
    #[arg(long, requires = "grp_from_bcf", help_heading = "output option")]
    pub bcf_to_bin_file: bool,

    // ---- hmm options
    /// Optional: Maximum number of fit iterations
    #[arg(short = 'm', long, default_value = "5", help_heading = "hmm option")]
    pub max_iter: u32,

    ///  Optional: Cap on the number of generations (floating point). Sets the maximum
    ///  value for that parameter in the fit. This is useful if you are
    ///  interested in recent IBD and are working with a population with
    ///  substantial linkage disequilbrium. Specifying a small value will force
    ///  the program to assume little recombination and thus a low transition
    ///  rate; otherwise it will identify the small blocks of LD as ancient IBD,
    ///  and will force the number of generations to be large.
    #[arg(short = 'n', long, default_value = "Inf", help_heading = "hmm option")]
    pub k_rec_max: f64,

    /// error rate in genotype calls
    #[arg(long, default_value_t = 0.001, help_heading = "hmm option")]
    pub eps: f64,

    /// minimum number of informative sites in a pairwise comparison (those
    /// with minor allele)
    #[arg(long, default_value_t = 10, help_heading = "hmm option")]
    pub min_inform: u32,

    /// minimum discordance rate in comparison; set > 0 to skip identical pairs
    #[arg(long, default_value_t = 0.0, help_heading = "hmm option")]
    pub min_discord: f64,

    /// maximum discordance rate in comparison; set < 1 to skip unrelated pairs
    #[arg(long, default_value_t = 1.0, help_heading = "hmm option")]
    pub max_discord: f64,

    /// skip next snp(s) if too close to last one; in bp
    #[arg(long, default_value_t = 5, help_heading = "hmm option")]
    pub min_snp_sep: u32,

    /// covergence criteria: min delta pi
    #[arg(long, default_value_t = 0.001, help_heading = "hmm option")]
    pub fit_thresh_dpi: f64,

    /// covergence criteria: min delta k_rec
    #[arg(long, default_value_t = 0.01, help_heading = "hmm option")]
    pub fit_thresh_dk: f64,

    /// covergence criteria: min (delta k_rec) / k_rec
    #[arg(long, default_value_t = 0.001, help_heading = "hmm option")]
    pub fit_thresh_drelk: f64,

    /// max number of unique alleles per site
    #[arg(long, default_value_t = 8, help_heading = "memory option")]
    pub max_all: u32,

    // ---- recombination rate and rate map
    #[command(flatten)]
    pub rec_args: RecombinationArg,

    // ---- memory buffer
    /// output buffer size for IBD segments, by default it is 8Kb.
    /// When IO is slow or not working well for writing many small chunks of data, it can
    /// be beneficial to set this to a larger value to reduce of the number  of
    /// times to call system IO operation
    #[arg(long, help_heading = "memory option")]
    pub buffer_size_segments: Option<usize>,

    /// output buffer size for IBD fraction records,  by default it is 8Kb.
    /// used similarly to --buffer-size-segments option
    #[arg(long, help_heading = "memory option")]
    pub buffer_size_frac: Option<usize>,

    // -- ibdseg filtering
    /// filtering IBD segments: if set, none of IBD/DBD segments short than
    /// filt_min_seg_cm will not be written to hmm.txt files
    #[arg(long, help_heading = "output segment filtering option")]
    pub filt_min_seg_cm: Option<f64>,
    /// filtering IBD segments: if set, none of IBD/DBD segments from pairs with
    /// k_rec > filt_max_tmrca will not be written to hmm.txt files
    #[arg(long, help_heading = "output segment filtering option")]
    pub filt_max_tmrca: Option<f64>,
    /// filtering IBD segments: if set, no DBD (non-IBD) segments will not be written to hmm.txt files
    #[arg(long, help_heading = "output segment filtering option")]
    pub filt_ibd_only: bool,

    // ---- parallelization
    /// number of threads. 0 : use all cpus; non-zero: use the given numbers of threads
    #[arg(long, default_value_t = 0, help_heading = "parallelization option")]
    pub num_threads: usize,

    /// par_mode: 0, chunks of sample pairs; 1, pairs of chunks.
    #[arg(long, default_value_t = 0, help_heading = "parallelization option")]
    pub par_mode: u8,

    /// number of pairs of samples per chunk for parallelization mode 0, or
    /// number of samples per chunk for parallelization mode 1
    #[arg(long, default_value_t = 120, help_heading = "parallelization option")]
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
    /// recombination rate per generation per basepair. When used --genome should not be specified.
    #[arg(
        short = 'r',
        long,
        default_value_t = 7.4e-7,
        group = "recombination_options",
        help_heading = "hmm option"
    )]
    pub rec_rate: f64,
    /// genome file which specifies chromosome size, names, and plink genetic map file paths.
    ///  When used --rec-rate should not be specified.
    #[arg(long, group = "recombination_options", help_heading = "hmm option")]
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
