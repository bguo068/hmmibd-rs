# hmmibd-rs

`hmmibd-rs` is a Rust rewrite of the original hmmIBD.c program, aimed at improving
readability and efficiency. The rewrite maintains the same hmm model and related
mathematical calculations as the C version. However, it introduces several
enhancements and features to enhance the functionality and usability of the
program.

## Features
The hmmibd-rs program offers the following features:

- Multithreading: The Rust version utilizes multithreading to accelerate the
inference process, allowing for faster computations and improved performance.
- Flexible Recombination Rate: Unlike the original C version, the Rust version
accepts recombination maps, enabling non-constant recombination rates. This
flexibility provides more accurate modeling and analysis of genetic data.
- Segment Filtering: The program includes the ability to filter IBD (Identity by
Descent) segments based on various criteria. This includes filtering by segment
type (IBD/non-IBD), length in centimorgans (cM), and estimated Time to Most
Recent Common Ancestor (TMRCA). This filtering capability allows users to focus
on specific segments of interest for further analysis.
- Command Line Configuration: In contrast to the C version, where some parameters are
hardcoded in the source code, the Rust version allows users to specify all
parameters from the command line. This includes the recombination rate and
fitting criteria, providing greater flexibility and customization

## Compile
```sh
cargo build --release --bin hmmibd2
```
The executable will be located at `target/release/hmmibd2`.

## Usage
```
Usage: hmmibd2 [OPTIONS] --data-file1 <DATA_FILE1>

Options:
  -i, --data-file1 <DATA_FILE1>
          File of genotype data. Format for genotype file:
          tab-delimited text file, with one single nucleotide
          polymorphism (SNP) per line. The first two columns are the
          chromosome and position, followed by one sample per column.
          A header line, giving the sample names, is required.
          Genotypes are coded by number: -1 for missing data, 0 for
          the first allele, 1 for the second, etc. SNPs and indels
          (if you trust them) can thus be treated on an equal
          footing. The variants must be in chromosome and position
          order, and can have between two and eight alleles (more, if
          you feel like changing max_allele)
  -I, --data-file2 <DATA_FILE2>
          Optional: file of genotype data from a second population
  -f, --freq-file1 <FREQ_FILE1>
          Optional: File of allele frequencies for the sample
          population. Format: tab-delimited, no header, one variant
          per row. Line format: <chromosome (int)> <position (bp,
          int)> <allele 1 freq> <all 2 freq> [<all 3 freq>] ... The
          genotype and frequency files must contain exactly the same
          variants, in the same order. If no file is supplied, allele
          frequencies are calculated from the input data file
  -F, --freq-file2 <FREQ_FILE2>
          Optional: File of allele frequencies for the second
          population; same format as for -f
  -m, --max-iter <MAX_ITER>
          Optional: Maximum number of fit iterations [default: 5]
  -b, --bad-file <BAD_FILE>
          Optional: file of sample ids to exclude from all analysis.
          Format: no header, one id (string) per row. Note: b stands
          for "bad samples"
  -g, --good-file <GOOD_FILE>
          Optional: File of sample pairs to analyze; all others are
          not processed by the HMM (but are still used to calculate
          allele frequencies). Format: no header, tab-delimited, two
          sample ids (strings) per row. Note: "g" stands for "good
          pairs"
  -n, --k-rec-max <K_REC_MAX>
          Optional: Cap on the number of generations (floating
          point). Sets the maximum value for that parameter in the
          fit. This is useful if you are interested in recent IBD and
          are working with a population with substantial linkage
          disequilbrium. Specifying a small value will force the
          program to assume little recombination and thus a low
          transition rate; otherwise it will identify the small
          blocks of LD as ancient IBD, and will force the number of
          generations to be large [default: Inf]
  -o, --output <OUTPUT>
          output put prefix, if not specified use the prefix as `-i`
          option
      --eps <EPS>
          error rate in genotype calls [default: 0.001]
      --min-inform <MIN_INFORM>
          minimum number of informative sites in a pairwise
          comparison (those with minor allele) [default: 10]
      --min-discord <MIN_DISCORD>
          minimum discordance rate in comparison; set > 0 to skip
          identical pairs [default: 0]
      --max-discord <MAX_DISCORD>
          maximum discordance rate in comparison; set < 1 to skip
          unrelated pairs [default: 1]
      --min-snp-sep <MIN_SNP_SEP>
          skip next snp(s) if too close to last one; in bp [default:
          5]
      --fit-thresh-dpi <FIT_THRESH_DPI>
          covergence criteria: min delta pi [default: 0.001]
      --fit-thresh-dk <FIT_THRESH_DK>
          covergence criteria: min delta k_rec [default: 0.01]
      --fit-thresh-drelk <FIT_THRESH_DRELK>
          covergence criteria: min (delta k_rec) / k_rec [default:
          0.001]
      --max-all <MAX_ALL>
          max number of unique alleles per site [default: 8]
  -r, --rec-rate <REC_RATE>
          recombination rate per generation per basepair. When used
          --genome should not be specified [default: 0.00000074]
      --genome <GENOME>
          genome file which specifies chromosome size, names, and
          plink genetic map file paths. When used --rec-rate should
          not be specified
      --filt-min-seg-cm <FILT_MIN_SEG_CM>
          filtering IBD segments: if set, none of IBD/DBD segments
          short than filt_min_seg_cm will not be written to hmm.txt
          files
      --filt-max-tmrca <FILT_MAX_TMRCA>
          filtering IBD segments: if set, none of IBD/DBD segments
          from pairs with k_rec > filt_max_tmrca will not be written
          to hmm.txt files
      --filt-ibd-only
          filtering IBD segments: if set, no DBD (non-IBD) segments
          will not be written to hmm.txt files
      --num-threads <NUM_THREADS>
          number of threads. 0 : use all cpus; non-zero: use the given 
          numbers of threads [default: 0]
      --par-chunk-size <PAR_CHUNK_SIZE>
          number of pairs of samples per chunk for parallelization 
          [default: 120]
      --par-mode <PAR_MODE>
          par_mode: 0, chunks of sample pairs; 1, pairs of chunks 
          [default: 0]
  -h, --help
          Print help
```

