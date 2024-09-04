# hmmibd-rs

`hmmibd-rs` is an enhanced implementation of the original
[hmmIBD](https://github.com/glipsnort/hmmIBD) program, designed to improve
computational efficiency and incorporate new features. While retaining the same
HMM model and related mathematical calculations as the original C version,
`hmmibd-rs` integrates the recombination rate map and enables parallelization.
Furthermore, it introduces several features to enhance the program's
functionality and usability.

## Features

The hmmibd-rs program offers the following features:

- Multithreading: The Rust version utilizes multithreading to accelerate the
  inference process, allowing for faster computations and improved performance.
- Recombination Rate Map: Unlike the original C version, the Rust version allows
  users to specify recombination rate maps and incorporate non-uniform
  recombination rates along the chromosomes/genomes. This flexibility can improve
  the accuracy of IBD detection and IBD segment filtration by length, especially
  when the genome is known to recombine non-uniformly and a recombination rate map
  is available.
- Using BCF Files as Input: The Rust implementation allows the direct use of BCF
  format as input. When provided with unphased genotype data in BCF, `hmmibd-rs`
  will utilize 'FORMAT/AD' to infer the dominant allele for each site/sample and
  use these dominant alleles as phased genotypes. Additionally, `hmmibd-rs` allows
  further filtering by sample and site missingness.
- Segment Filtering: This program can filter Identity by Descent (IBD) segments
  based on various criteria, including segment type (IBD/non-IBD), length in
  centimorgans (cM), and estimated Time to Most Recent Common Ancestor (TMRCA)
  before reporting. This capability enables users to focus on specific segments of
  interest for further analysis and significantly reduce storage requirements when
  working with large data sets.
- Command Line Interface: Unlike the C version, where some parameters are
  hardcoded in the source code, the Rust version allows users to specify all
  parameters from the command line, including the recombination rate and fitting
  criteria. This offers greater flexibility and customization. Despite the
  changes, the Rust command line interface is designed to be compatible with
  the C implementation, allowing the same options from the C version to be used
  similarly in the Rust version.

## Compilation and installation

```sh
git clone https://github.com/bguo068/hmmibd-rs.git
cd hmmibd-rs
cargo build --release
```

The executable will be located at `target/release/hmmibd-rs` after the build.

This executable can be moved to any path specified in the `$PATH` environment
variable for convenience. Otherwise, you can specify the path to the executable
and run it.

## Usage

### Command line options

```
Usage: hmmibd-rs [OPTIONS] --data-file1 <DATA_FILE1>

Options:
  -h, --help     Print help
  -V, --version  Print version

input data:
  -i, --data-file1 <DATA_FILE1>
          File of genotype data. (1) by default, format for genotype file:
          tab-delimited text file, with one single nucleotide polymorphism (SNP)
          per line. The first two columns are the chromosome and position,
          followed by one sample per column. A header line, giving the sample
          names, is required. Genotypes are coded by number: -1 for missing
          data, 0 for the first allele, 1 for the second, etc. SNPs and indels
          (if you trust them) can thus be treated on an equal footing. The
          variants must be in chromosome and position order, and can have
          between two and eight alleles (more, if you feel like changing
          max-allele). (2) When `--from-bcf` is specified, the input is expected
          in BCF format, either from a file or stdin. If the argument is a file
          path, the BCF genotype is read from that file. If the argument is "-",
          the BCF genotype is read from stdin. (3) When `--from-bin` is
          specified, a binary genotype input file is expected. See options
          `--bcf-to-bin-file-by-chromosome` or `bcf-to-bin-file` for generating
          binary genotype files
  -I, --data-file2 <DATA_FILE2>
          Optional: file of genotype data from a second population
  -f, --freq-file1 <FREQ_FILE1>
          Optional: File of allele frequencies for the sample population.
          Format: tab-delimited, no header, one variant per row. Line format:
          <chromosome (int)> <position (bp, int)> <allele 1 freq> <all 2 freq>
          <all 3 freq> ... The genotype and frequency files must contain exactly
          the same variants, in the same order. If no file is supplied, allele
          frequencies are calculated from the input data file
  -F, --freq-file2 <FREQ_FILE2>
          Optional: File of allele frequencies for the second population; same
          format as for -f
  -b, --bad-file <BAD_FILE>
          Optional: file of sample ids to exclude from all analysis. Format: no
          header, one id (string) per row. Note: b stands for "bad samples"
  -g, --good-file <GOOD_FILE>
          Optional: File of sample pairs to analyze; all others are not
          processed by the HMM (but are still used to calculate allele
          frequencies). Format: no header, tab-delimited, two sample ids
          (strings) per row. Note: "g" stands for "good pairs"

input data options:
      --from-bcf
          Optional: flag indicating whether the input file is of BCF format
      --from-bin
          Optional: flag indicating whether the input file is binary genotype
          data
      --bcf-read-mode <BCF_READ_MODE>
          Optional: Read mode for the Bcf input. Currently "first-ploidy" and
          "each-ploidy" are experimental [default: dominant-allele] [possible
          values: dominant-allele, first-ploidy, each-ploidy]
      --bcf-filter-config <BCF_FILTER_CONFIG>
          Optional: When the --from_bcf flag is set and the --bcf-read-mode
          dominant-allele option (default setting) is used, this option allows
          specifying a custom BCF filtering configuration file. If not provided,
          a built-in filtering criteria will be used, and a temporary
          configuration file consistent with the built-in criteria named
          "tmp_bcf_filter_config.toml" will be generated in the current folder

output data:
  -o, --output <OUTPUT>  Optional: output prefix. If not specified, the prefix
                         from the `-i` option will be used. Note: When reading
                         genotype from stdin, this option must be specified

output options:
      --suppress-frac
          Optional: whether to suppress frac output. Toggle this on for minimal
          IO burden when only IBD segment output is needed
      --bcf-to-bin-file-by-chromosome
          Optional, valid only with the `--from-bcf` option. When enabled,
          genotypes are written to a binary file for each chromosome.
          Chromosomes with no sites are ignored. When this is set, the HMM IBD
          inference step is skipped. This option is designed to separate the BCF
          processing step from the HMM inference step, especially when genotypes
          of different chromosomes are concatenated for sample filtering, and
          the user wants to split the filtering genotype into chromosomes and
          run each instance of hmmibd-rs on a chromosome in parallel
      --bcf-to-bin-file
          Optional, valid only with the `--from-bcf` option. Similar to the
          `--bcf-to-bin-file-by-chromosome` option, but generates a single
          binary file containing all chromosomes

hmm options:
  -m, --max-iter <MAX_ITER>
          Optional: Maximum number of fit iterations [default: 5]
  -n, --k-rec-max <K_REC_MAX>
          Optional: Cap on the number of generations (floating point). Sets the
          maximum value for that parameter in the fit. This is useful if you are
          interested in recent IBD and are working with a population with
          substantial linkage disequilbrium. Specifying a small value will force
          the program to assume little recombination and thus a low transition
          rate; otherwise it will identify the small blocks of LD as ancient
          IBD, and will force the number of generations to be large [default:
          Inf]
      --eps <EPS>
          Optional: error rate in genotype calls [default: 0.001]
      --min-inform <MIN_INFORM>
          Optional: minimum number of informative sites in a pairwise comparison
          (those with minor allele) [default: 10]
      --min-discord <MIN_DISCORD>
          Optional: minimum discordance rate in comparison from 0.0 to 1.0;  set
          > 0 to skip identical pairs [default: 0]
      --max-discord <MAX_DISCORD>
          Optional: maximum discordance rate in comparison from 0.0 to 1.0;  set
          < 1 to skip unrelated pairs [default: 1]
      --min-snp-sep <MIN_SNP_SEP>
          Optional: skip next snp(s) if too close to last one; in bp [default:
          5]
      --fit-thresh-dpi <FIT_THRESH_DPI>
          Optional: covergence criteria: min delta pi [default: 0.001]
      --fit-thresh-dk <FIT_THRESH_DK>
          Optional: covergence criteria: min delta k_rec [default: 0.01]
      --fit-thresh-drelk <FIT_THRESH_DRELK>
          Optional: covergence criteria: min (delta k_rec) / k_rec [default:
          0.001]
  -r, --rec-rate <REC_RATE>
          Optional: constant recombination rate per generation per base pair.
          When used, --genome should not be specified [default: 0.00000074]
      --genome <GENOME>
          Optional: Genome file specifying chromosome sizes, names, and paths to
          PLINK genetic map files. When used, do not specify --rec-rate

memory options:
      --max-all <MAX_ALL>
          Optional: Maximum number of unique alleles per site. For instance, if
          the input contains only biallelic variants, setting `--max-all 2`
          could reduce memory usage and enhance computation speed by utilizing a
          more compact memory layout [default: 8]
      --buffer-size-segments <BUFFER_SIZE_SEGMENTS>
          Optional: output buffer size in bytes for IBD segments, by default it
          is 8Kb. When IO is slow or not working well for writing many small
          chunks of data, it can be beneficial to set this to a larger value to
          reduce of the number  of times to call system IO operation
      --buffer-size-frac <BUFFER_SIZE_FRAC>
          Optional: output buffer size in bytes for IBD fraction records,  by
          default it is 8Kb. used similarly to --buffer-size-segments option

output segment filtering options:
      --filt-min-seg-cm <FILT_MIN_SEG_CM>
          Optional: filtering IBD segments. If set, IBD/non-IBD segments short
          than FILT-MIN-SEG-CM cM will not be written to hmm.txt files
      --filt-max-tmrca <FILT_MAX_TMRCA>
          Optional: filtering IBD segments. If set, IBD/non-IBD segments from
          pairs with k_rec > filt_max_tmrca will not be written to hmm.txt files
      --filt-ibd-only
          Optional: filtering IBD segments. If set, non-IBD segments will not be
          written to hmm.txt files

parallelization options:
      --num-threads <NUM_THREADS>
          Optional: number of threads. "0" : use all cpus; non-zero: use the
          given numbers of threads [default: 0]
      --par-mode <PAR_MODE>
          Optional: parallelization mode. Mode "0" processes chunks of sample
          pairs by creating a vector of sample pairs  and splitting it into
          equal lengths, as specified by --par-chunk-size. Different threads
          handle these chunks in parallel without coping data. This mode is
          suitable for smaller sample sizes. Mode "1" processes pairs of sample
          chunks by creating a vector of samples, splitting it into chunks, and
          then creating a second vector of pairs of these chunks. Different
          threads process these pairs of sample chunks in parallel. This mode is
          ideal for larger sample sizes, where genotype data for sample pairs
          may be dispersed, potentially slowing down the data reading step. To
          address this, Mode "1" first copies each pair of chunks to the thread
          heap, ensuring that samples within these pairs are close together,
          thus benefiting from better memory cache locality [default: 0]
      --par-chunk-size <PAR_CHUNK_SIZE>
          Optional: number of sample pairs per chunk for parallelization mode 0,
          or number of samples per chunk for parallelization mode 1 [default:
          120]
```

### Examples

1. Use `hmmibd-rs` similarly to `hmmIBD`

```sh
target/release/hmmibd-rs \
    -i c/samp_data/pf3k_Cambodia_13.txt  \
    -f c/samp_data/freqs_pf3k_Cambodia_13.txt  \
    -o tmp  -r 6.66667e-7
target/release/hmmibd-rs \
    -i c/samp_data/pf3k_Cambodia_13.txt \
    -I c/samp_data/pf3k_Ghana_13.txt \
    -f c/samp_data/freqs_pf3k_Cambodia_13.txt \
    -F c/samp_data/freqs_pf3k_Ghana_13.txt \
    -o tmp  -r 6.66667e-7
```

2. Specifiy multi-threading options

```sh
target/release/hmmibd-rs \
    -i c/samp_data/pf3k_Cambodia_13.txt  \
    -f c/samp_data/freqs_pf3k_Cambodia_13.txt \
    --par-mode 1 \
    --par-chunk-size 10 \
    -o tmp  -r 6.66667e-7
```
Please note there are two options for `--par-mode`. Mode 0 is generally suitable
for small sample sizes, while mode 1 is better for large sample sizes. See the
corresponding help message by running `hmmibd-rs -h`.

3. Use recombination rate map

```sh
# for data with a single chromsome
target/release/hmmibd-rs \
    -i testdata/sim_data/gt_chr1.txt \
    --genome testdata/sim_data/genome.toml \
    -o tmp

# for data with multiple chromsomes
cat testdata/sim_data/gt_chr{1,2,3}.txt >  testdata/sim_data/gt_all_chrs.txt
target/release/hmmibd-rs \
    -i testdata/sim_data/gt_all_chrs.txt \
    --genome testdata/sim_data/genome.toml \
    -o tmp
```

Note: The `genome.toml` file contains information about chromosome names, sizes,
and Plink genetic map files. Examples:
[genome file](testdata/sim_data/genome.toml),
and plink files for
[constant rate](testdata/sim_data/map/) and
[nonuniform rate](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip)

4. Filter IBD results before reporting

```sh
target/release/hmmibd-rs \
    -i testdata/sim_data/gt_all_chrs.txt \
    --suppress-frac \
    --filt-ibd-only \
    --filt-min-seg-cm 2.0 \
    -o tmp
```

5. Use BCF as input from files

   - Create a bcf_filt_config file. Details of this config files, see section
     `BCF Processing Details and BCF Filtering Configuration` section below.

   ```sh
   # create bcf_filt_config file
   echo '
   min_depth = 5
   min_ratio = 0.7
   min_maf = 0.01
   filter_column_pass_only = true
   major_minor_alleles_must_be_snps = true
   min_r1_r2 = 3.0
   min_site_nonmissing = 0.1
   nonmissing_rates= [
       [0.41, 0.205],
       [0.51, 0.255],
       [0.61, 0.305],
       [0.71, 0.355],
       [0.81, 0.405],
       [0.85, 0.45],
       [0.85, 0.50],
   ]
   ' > tmp_bcf_filt_config.toml
   ```

   - Use bcf file as input for IBD detection

   ```sh
   target/release/hmmibd-rs \
       -i testdata/pf7_data/pf7_chr1_20samples_maf0.01.bcf \
       --from-bcf \
       --bcf-read-mode dominant-allele \
       --bcf-filter-config  tmp_bcf_filt_config.toml \
       -o tmp
   ```

6. Use BCF as input from stdin. This can be useful in certain cases:
   - When genotype data exists in multiple files (such as each chromosome
     having a BCF file), you might want to use `bcf concat` to combine all these
     files and pipe the data to `hmmibd-rs` without writing the concatenated
     files to disk first.
   - When the input BCF requires additional manipulation, you can first modify
     the input BCF as needed and then pipe it into `hmmibd-rs`, thereby avoiding
     the creation of intermediate BCF files.

```sh
bcftools view testdata/pf7_data/pf7_chr1_20samples_maf0.01.bcf -Ou | \
target/release/hmmibd-rs \
    -i - \
    --from-bcf \
    --bcf-read-mode dominant-allele \
    --bcf-filter-config  tmp_bcf_filt_config.toml \
    -o tmp
```

7. Use binary genotype as input

   - create bcf binary files

   ```sh
   target/release/hmmibd-rs \
       -i testdata/pf7_data/pf7_chr1_20samples_maf0.01.bcf \
       --from-bcf \
       --bcf-read-mode dominant-allele \
       --bcf-filter-config  tmp_bcf_filt_config.toml \
       --bcf-to-bin-file-by-chromsome
       -o tmp_bcf_bin_file
   ```

   - use bcf binary files for IBD calling

   ```sh
   target/release/hmmibd-rs \
       -i tmp_bcf_bin_file_Pf3D7_01_v3.bin  \
       --from-bin \
       -o tmp
   ```

8. More complete examples

Please check the
[hmmibd-rs-bench-empirical](https://github.com/bguo068/hmmibd-rs-bench-empirical)
repository.

Detailed documentation of all related command line options can be obtained via
`target/release/hmmibd-rs --help`.

## BCF Processing Details and BCF Filtering Configuration

### Requirment on the input BCF file

Currently, `hmmibd-rs` assumes the following about the input BCF files:

- Each position (CHROM/POS pair) occurs only once throughout the file. This
  implies that each site corresponds to a single BCF record. If multiallelic
  records in BCF files have been converted to multiple biallelic records, you will
  need to use `bcftools norm -m+` to merge these records back into multiallelic
  records.
- Sites are sorted by position within each chromosome.
- When using `--bcf-read-mode dominant-allele` (set by default), the `FORMAT/AD`
  fields must be present in the BCF file.

### Part 1. Read through the bcf files and filter sites

1. Read genotype data from the (FORMAT/AD) columns for targeted samples:

   - In the config file, add a line
     `target_samples = "path/to/sample_list_file"`
     to specify targeted samples. In the sample list
     file, each sample name should be on a separate line. For example, targeted
     samples could be the monoclonal samples with `Fws > 0.95`.
   - If `target_samples` is not specified, all samples in the BCF file will be
     included.

2. If a site has empty ALT columns, it is skipped.

3. At the current site, for each sample, use FORMAT/AD to determine the dominant allele:

   - The major and minor alleles (for a sample) are defined as the alleles with
     the largest and second largest AD values, respectively.
   - The following thresholds are used to determine if the major allele is
     dominant:
     - `min_depth`: Minimum total AD for the sample. To be dominant, the
       major allele must have `total-AD` ≥ `min_depth`.
     - `min_ratio`: Minimum value
       for the ratio of `major-allele-AD` to `total-AD`. To be dominant,
       `major-allele-AD / total-AD` must be ≥ `min_ratio`.
     - `min_r1_r2`: Minimum
       ratio of major allele AD to minor allele AD. To be dominant,
       `minor-allele-AD / major-allele-AD` must be <= 1.0 / `min_r1_r2`.
     - If any of these conditions are false, the sample is treated as having a
       missing genotype for this site. Otherwise, the major allele is designated
       as the dominant allele representing the genotype for this sample at the
       current site.

4. Once the genotype represented by dominant alleles or missing status is
   determined for all samples at the current site, the following tests are used to
   filter the site:
   - Calculate the `non_missing_rate` as the fraction of samples with dominant
     alleles. If the site has `non_missing_rate` < `min_site_nonmissing`, the
     site will be filtered out.
   - Calculate the minor allele frequency `maf` for the site using dominant
     allele genotypes. If `maf` < `min_maf`, the site will be filtered out.
   - Calculate the allele counts across samples for the site and determine the
     major and minor alleles (for a site) as the alleles with the largest and
     second largest sample counts. If `major_minor_alleles_must_be_snps` is set
     to `true` in the configuration, the site will be filtered out if either the
     major or minor allele is not a SNP.

### Part 2. Iteratively Filter Sites and Samples by Genotype Missingness

This part involves an iterative process of removing sites and samples with the
highest genotype missing rates (or lowest non-missing rates). The thresholds can
be configured via the variable `nonmissing_rates`, which is a list of
two-element arrays. The first element of each array is the
`min_site_nonmiss_rate`, and the second is the `min_sample_nonmiss_rate`.

For example:

```toml
nonmissing_rates= [
    [0.31, 0.155],
    [0.35, 0.175],
    ...
    [0.90, 0.80],
]
```

We can gradually increase either or both of the two thresholds over iterations
to avoid aggressively removing sites or samples by an abrupt stringent
threshold. The iterative filtering process should be adjusted according to the
data to ensure that a desired number of samples and sites are retained for IBD
detection.

For a full example of the bcf_filter_config.toml file, please check [this
file](https://github.com/bguo068/hmmibd-rs-bench-empirical/blob/main/input/bcf_filter_config.toml)

## Output

The output files follow the same format as `hmmIBD`:

> Two output files are produced. The file _\<filename\>.hmm.txt_ contains a list
> of all segments (where a segment is one or more contiguous variant sites in
> the same state) for each sample pair, the assigned state (IBD or not-IBD) for
> each, and the number of variants covered by the segment; note that an assigned
> state of 0 means IBD, while 1 means not-IBD. These segments represent the most
> probable state assignments.
>
> The file _\<filename\>.hmm_fract.txt_ summarizes results for each sample pair
> (including some information that may be of interest only to me). If you are only
> interested in the fraction of the genome that is IBD between a pair of samples,
> look at the last column (**fract_sites_IBD**). Columns:
>
> - **N_informative_sites**: Number of sites with data for both samples, and with
>   at least one copy of the minor allele.
> - **discordance**: Fraction of informative sites that have different alleles in
>   the two samples.
> - **log_p**: natural logarithm of the probability of the final set of state
>   assignments and the set of observations (calculated with the Viterbi algorithm).
> - **N_fit_iteration**: Number of iterations carried out in EM fitting.
> - **N_generation**: Estimated number of generations of recombination (1 of 2
>   free parameters in fit).
> - **N_state_transition**: Number of transitions between IBD and not-IBD states
>   across entire genome.
> - **seq_shared_best_traj**: Fraction of sequence IBD based on the best state
>   assignment, calculated as (seq in IBD segments) / (seq in IBD segments + seq in
>   not-IBD segments). Segments in which there is a state transition between IBD and
>   not-IBD are ignored.
> - **fract_sites_IBD**: Fraction of variant sites called IBD calculated for all
>   possible state assignments, weighted by their probability (equal to the marginal
>   posterior probability of the IBD state calculated with the forward-backward
>   algorithm).
> - **fract_vit_sites_IBD**: Fraction of variant sites called IBD calculated
>   for the best state assignment (i.e. the result of the Viterbi algorithm, as in
>   seq_shared_best_traj).
