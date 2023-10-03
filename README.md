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

See below for original documentation from the C version

# hmmIBD
A hidden Markov model for detecting segments of shared ancestry (identity by descent) in genetic sequence data.

## Overview

The program *hmmIBD* implements a hidden Markov model (HMM) for detecting genomic regions that are identical by descent (IBD) for pairs of haploid samples. It was written to find large IBD regions in sequenced haploid *P. falciparum* genomes, but it can be applied to other organisms (including phased diploids) and can find shorter IBD regions as well. The program takes as input a file of genotype calls for a set of samples, assumed to be from a single population. As of version 2.0.0, the program will also accept a second file of genotype calls, which are treated as coming from a different population with different allele frequencies. For a single population, all pairwise comparisons are made between the samples (unless otherwise specified with a -b or -g flag.) For two populations, all comparisons are made between samples from different populations. 

The details of the model and program are described in a manuscript in preparation. 

Under the HMM, each variant site is assumed to be in one of two hidden states, IBD or not-IBD.  To calculate the probability of each state, estimates of the allele frequencies for every variant are required.  By default, they are calculated from the input data, but a separate file of allele frequencies can be supplied by the user (preferable if analyzing a subset of the data).

The model has two free parameters, (1) the fraction of the genome that is IBD, and (2) the number of generations during which recombination has been breaking down IBD blocks. (Note: the former is generally estimated more accurately than the latter, and is relatively robust to the latter.) The program fits for optimal values of these parameters by an iterative estimation-maximization procedure. Iterations of the fit are capped at a user-settable maximum (default = 5). To accurately determine the IBD fraction for large shared chromosome segments, only a few iterations are needed, while for smaller, older blocks of IBD, the fit may continue to improve for 15 or more iterations.

The Viterbi algorithm calculates the best single set of state assignments given data under the HMM and outputs that set. The forward-backward algorithm sums the fraction of the genome that is IBD over all possible state assignments given data under the HMM, weighting each by the probability of that set of states. If you are interested in the IBD fraction, rather than precisely which parts of the genome are IBD, this is probably the output you want (see **fract_sites_IBD** in the Output section).

## Installation

The C source code must be compiled before use, and requires no special libraries. Using the current OS X built-in C compiler, I compile it with the command

```
cc -o hmmIBD -O3 -lm -Wall hmmIBD.c
```

which should yield the executable file *hmmIBD*, but other compilers should work as well.

## Execution

hmmIBD is run from the command line. It requires the user to supply two options when invoked; it can also take seven optional arguments:

```
hmmIBD -i <input filename (for pop1, if using 2 pops)> -o <output filename> 
     [-I <input file, pop2>] [-f <allele frequency file (pop1)>] [-F <allele freq file (pop2)>]
     [-b <file with samples to skip>] [-m <max fit iteration>] [-n <max N generation>] [-g <file with sample pairs to use>]
```
Required options:
- -i: File of genotype data. See below for format.
- -o: Output file name. Two output files will be produced, with ".hmm.txt" 
      and ".hmm_fract.txt" appended to the supplied name. See below for details

Optional options:
- -f: File of allele frequencies for the sample population. Format: tab-delimited, no header, one variant per row. Line format: `<chromosome (int)> <position (bp, int)> <allele 1 freq> <all 2 freq> [<all 3 freq>] ...` The genotype and frequency files must contain exactly the same variants, in the same order. If no file is supplied, allele frequencies are calculated from the input data file.
- -I: File of genotype data from a second population; same format as for -i. (added in 2.0.0)
- -F: File of allele frequencies for the second population; same format as for -f. (added in 2.0.0)
- -m: Maximum number of fit iterations (defaults to 5).
- -b: File of sample ids to exclude from all analysis. Format: no header, one id (string) per row. (Note: b stands for "bad samples".)
- -g: File of sample pairs to analyze; all others are not processed by the HMM 	(but are still used to calculate allele frequencies). Format: no header,	tab-delimited, two sample ids (strings) per row. (Note: "g" stands for 	"good pairs".)
- -n: Cap on the number of generations (floating point). Sets the maximum value for that parameter in the fit. This is useful if you are interested in recent IBD and are working with a population with substantial linkage disequilbrium. Specifying a small value will force the program to assume little recombination and thus a low transition rate; otherwise it will identify the small blocks of LD as ancient IBD, and will force the number of generations to be large.

## Input file formats

Format for genotype file: tab-delimited text file, with one single nucleotide polymorphism (SNP) per line. The first two columns are the chromosome and position, followed by one sample per column. A header line, giving the sample names, is required. Genotypes are coded by number: -1 for missing data, 0 for the first allele, 1 for the second, etc. SNPs and indels (if you trust them) can thus be treated on an equal footing. The variants must be in chromosome and position order, and can have between two and eight alleles (more, if you feel like changing *max_allele* in the code).

## Output files

Two output files are produced. The file *\<filename\>.hmm.txt* contains a list of all segments (where a segment is one or more contiguous variant sites in the same state) for each sample pair, the assigned state (IBD or not-IBD) for each, and the number of variants covered by the segment; note that an assigned state of 0 means IBD, while 1 means not-IBD. These segments represent the most probable state assignments.

The file *\<filename\>.hmm_fract.txt* summarizes results for each sample pair (including some information that may be of interest only to me). If you are only interested in the fraction of the genome that is IBD between a pair of samples, look at the last column (**fract_sites_IBD**). Columns:

- **N_informative_sites**: Number of sites with data for both samples, and with at least one copy of the minor allele.
- **discordance**: Fraction of informative sites that have different alleles in the two samples.
- **log_p**: natural logarithm of the probability of the final set of state assignments and the set of observations (calculated with the Viterbi algorithm).
- **N_fit_iteration**: Number of iterations carried out in EM fitting. 
- **N_generation**: Estimated number of generations of recombination (1 of 2 free parameters in fit).
- **N_state_transition**: Number of transitions between IBD and not-IBD states across entire genome.
- **seq_shared_best_traj**: Fraction of sequence IBD based on the best state assignment, calculated as (seq in IBD segments) / (seq in IBD segments + seq in not-IBD segments). Segments in which there is a state transition between IBD and not-IBD are ignored. 
- **fract_sites_IBD**: Fraction of variant sites called IBD calculated for all possible state assignments, weighted by their probability (equal to the marginal posterior probability of the IBD state calculated with the forward-backward algorithm).
- **fract_vit_sites_IBD**: Fraction of variant sites called IBD calculated for the best state assignment (i.e. the result of the Viterbi algorithm, as in seq_shared_best_traj).


## Some variables you might want to change

The following variables control program execution in various ways, and can be easily changed in the C code prior to compilation. 

- *eps* -- assumed genotype error rate (rate of calling allele i as allele j) (default = 0.1%).
- *min_inform* -- minimum number of informative sites (those with at least one copy of the minor allele in this pair) required for processing sample pair.
- *max_discord* -- maximum fraction of informative sites with discordant genotypes between the two samples; used to skip unrelated pairs .
- *min_discord* -- minimum discordance, useful for skipping identical pairs.
- *nchrom* -- number of chromosomes in genome (14 for P. falciparum and P. vivax)
- *min_snp_sep* -- number of bp separation required between SNPs, to avoid mutations spanning >1 bp (default = 10).

## Sample data

I have included 2 small sample data set of *P. falciparum* sequences, drawn from the Pf3K Cambodia and Ghana samples. Each data file has genotype data for chromosomes 1-3 from 10 samples. The corresponding frequency files contain allele frequency information for the same variants based on the entire Pf3K dataset for that country. Executing the commands 

```
hmmIBD -i samp_data/pf3k_Cambodia_13.txt -f samp_data/freqs_pf3k_Cambodia_13.txt -o mytest_1pop
```

```
hmmIBD -i samp_data/pf3k_Cambodia_13.txt -f samp_data/freqs_pf3k_Cambodia_13.txt -I samp_data/pf3k_Ghana_13.txt -F samp_data/freqs_pf3k_Ghana_13.txt -o mytest_2pop
```

should give similar results to that found in the output files in *samp_data/*.

AUTHOR

Steve Schaffner (sfs@broadinstitute.org)

Please feel free to get in touch with comments, suggestions and questions. 
