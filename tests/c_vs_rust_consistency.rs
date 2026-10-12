//! Integration tests: hmmibd-rs must reproduce the original C hmmIBD (c/hmmIBD.c) exactly,
//! for both output files, on the bundled sample data.
//!
//! The C program is compiled with `cc` (or `$CC`) into the cargo target directory on first
//! use. Each case runs both programs on the same input and compares
//!   * `.hmm.txt`       line by line (order-insensitive), and
//!   * `.hmm_fract.txt` field by field with numeric parsing, so that the only accepted
//!     difference is float formatting (C prints `-2.54409e+03`, Rust prints `-2.54409e3`).
//! Cases are chosen so that each previously divergent code path is exercised:
//! seq_shared_best_traj / N_fit_iteration / N_generation (every case), the pair filter on
//! informative sites (`few_variants`), and the newline-delimited `-b` file (`bad_samples`).
//!
//! The C program has fixed-size filename buffers (out prefix + extension must fit in 64
//! bytes, input paths in 128), so both programs are run from the repository root with short
//! relative paths; temporary files live in `target/cvr/`.
//!
//! The compiled C program is the oracle; the files in c/samp_data/output_* were produced by an
//! older hmmIBD release and differ slightly from what c/hmmIBD.c produces, so they are not used.
//!
//! Run with `cargo test --release --test c_vs_rust_consistency`.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Once;

const C_SRC: &str = "c/hmmIBD.c";
const GT1: &str = "c/samp_data/pf3k_Cambodia_13.txt";
const FQ1: &str = "c/samp_data/freqs_pf3k_Cambodia_13.txt";
const GT2: &str = "c/samp_data/pf3k_Ghana_13.txt";
const FQ2: &str = "c/samp_data/freqs_pf3k_Ghana_13.txt";

static BUILD_C: Once = Once::new();

/// Repository root; every command runs with this as its working directory.
fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Short relative directory for the compiled C binary and all temporary files: target/cvr/
const WORK: &str = "target/cvr";

fn work_dir() -> PathBuf {
    let dir = root().join(WORK);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// Path of a work file as the programs see it (relative to the repository root).
fn rel(name: &str) -> String {
    format!("{WORK}/{name}")
}

fn c_binary() -> PathBuf {
    let bin = work_dir().join("hmmIBD_c");
    BUILD_C.call_once(|| {
        let cc = std::env::var("CC").unwrap_or_else(|_| "cc".to_string());
        let status = Command::new(&cc)
            .current_dir(root())
            .args(["-O2", "-o"])
            .arg(&bin)
            .args([C_SRC, "-lm"])
            .status()
            .unwrap_or_else(|e| panic!("cannot run C compiler `{cc}` ({e}); set $CC"));
        assert!(status.success(), "compiling {C_SRC} failed");
    });
    bin
}

fn run(bin: &Path, args: &[&str], out_prefix: &str) {
    let out = Command::new(bin)
        .current_dir(root())
        .args(args)
        .arg("-o")
        .arg(out_prefix)
        .output()
        .unwrap_or_else(|e| panic!("cannot run {}: {e}", bin.display()));
    assert!(
        out.status.success(),
        "{} {:?} failed:\n{}",
        bin.display(),
        args,
        String::from_utf8_lossy(&out.stderr)
    );
}

fn body_lines(path: &Path) -> Vec<String> {
    let path = root().join(path);
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines().skip(1).map(str::to_owned).collect()
}

/// A .hmm_fract.txt row with numeric columns parsed, so formatting differences vanish.
#[derive(Debug, PartialEq, PartialOrd)]
struct FracRow {
    sample1: String,
    sample2: String,
    n_informative: u64,
    discordance: f64,
    log_p: f64,
    n_fit_iteration: u64,
    n_generation: f64,
    n_state_transition: u64,
    seq_shared_best_traj: f64,
    fract_sites_ibd: f64,
    fract_vit_sites_ibd: f64,
}

fn parse_frac(path: &Path) -> Vec<FracRow> {
    let mut rows: Vec<FracRow> = body_lines(path)
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 11, "{}: bad frac line {l:?}", path.display());
            let d = |i: usize| f[i].parse::<f64>().unwrap_or_else(|_| panic!("{}: field {i} {:?}", path.display(), f[i]));
            let u = |i: usize| f[i].parse::<u64>().unwrap_or_else(|_| panic!("{}: field {i} {:?}", path.display(), f[i]));
            FracRow {
                sample1: f[0].to_owned(),
                sample2: f[1].to_owned(),
                n_informative: u(2),
                discordance: d(3),
                log_p: d(4),
                n_fit_iteration: u(5),
                n_generation: d(6),
                n_state_transition: u(7),
                seq_shared_best_traj: d(8),
                fract_sites_ibd: d(9),
                fract_vit_sites_ibd: d(10),
            }
        })
        .collect();
    rows.sort_by(|a, b| a.partial_cmp(b).unwrap());
    rows
}

/// Run both programs with identical hmmIBD-style arguments and assert both outputs agree.
/// `expect_rows`: if Some(n), additionally assert the number of pairs written.
fn assert_same_outputs(case: &str, args: &[&str], expect_rows: Option<usize>) {
    work_dir();
    let c_prefix = rel(&format!("{case}_c"));
    let rs_prefix = rel(&format!("{case}_rs"));
    run(&c_binary(), args, &c_prefix);
    let mut rs_args = args.to_vec();
    rs_args.extend(["--num-threads", "1"]);
    run(Path::new(env!("CARGO_BIN_EXE_hmmibd-rs")), &rs_args, &rs_prefix);

    // segments: same multiset of lines
    let c_seg: BTreeSet<String> = body_lines(Path::new(&format!("{c_prefix}.hmm.txt"))).into_iter().collect();
    let rs_seg: BTreeSet<String> = body_lines(Path::new(&format!("{rs_prefix}.hmm.txt"))).into_iter().collect();
    assert_eq!(c_seg, rs_seg, "[{case}] .hmm.txt differs");

    // fraction file: same rows after numeric parsing
    let c_frac = parse_frac(Path::new(&format!("{c_prefix}.hmm_fract.txt")));
    let rs_frac = parse_frac(Path::new(&format!("{rs_prefix}.hmm_fract.txt")));
    assert_eq!(c_frac.len(), rs_frac.len(), "[{case}] .hmm_fract.txt row counts differ");
    for (c, r) in c_frac.iter().zip(rs_frac.iter()) {
        assert_eq!(c, r, "[{case}] .hmm_fract.txt row differs");
    }
    if let Some(n) = expect_rows {
        assert_eq!(c_frac.len(), n, "[{case}] unexpected number of pairs");
    }
}

#[test]
fn single_pop_with_freq_file() {
    assert_same_outputs("single_freq", &["-i", GT1, "-f", FQ1], Some(45));
}

#[test]
fn single_pop_freqs_inferred() {
    assert_same_outputs("single_nofreq", &["-i", GT1], Some(45));
}

#[test]
fn two_pops_with_freq_files() {
    assert_same_outputs("two_pop", &["-i", GT1, "-I", GT2, "-f", FQ1, "-F", FQ2], Some(100));
}

#[test]
fn max_iter_1() {
    // every pair stops after one iteration: exercises N_fit_iteration reporting when unconverged
    assert_same_outputs("max_iter_1", &["-i", GT1, "-f", FQ1, "-m", "1"], Some(45));
}

#[test]
fn bad_samples_newline_delimited() {
    // one id per line, as documented for hmmIBD; both programs must drop the same 2 samples (45 -> 28 pairs)
    fs::write(work_dir().join("bad.txt"), "PH0047-Cx\nPH0051-Cx\n").unwrap();
    assert_same_outputs("bad_samples", &["-i", GT1, "-f", FQ1, "-b", &rel("bad.txt")], Some(28));
}

#[test]
fn few_variants_pair_filter() {
    // 40 variants: no pair has >= 10 informative sites, so the pair filter must remove all of them in both
    let dir = work_dir();
    let take = |src: &str, n: usize| {
        fs::read_to_string(root().join(src)).unwrap().lines().take(n).map(|l| format!("{l}\n")).collect::<String>()
    };
    fs::write(dir.join("gt40.txt"), take(GT1, 41)).unwrap();
    fs::write(dir.join("fq40.txt"), take(FQ1, 40)).unwrap();
    assert_same_outputs("few_variants", &["-i", &rel("gt40.txt"), "-f", &rel("fq40.txt")], Some(0));
}

#[test]
fn good_pairs_file() {
    // pairs listed in header order. (A pair listed in reverse order is run by both programs but
    // printed with the samples swapped by hmmibd-rs; that formatting difference is not covered here.)
    fs::write(work_dir().join("good.txt"), "PH0047-Cx\tPH0051-Cx\nPH0055-C\tPH0057-C\n").unwrap();
    assert_same_outputs("good_pairs", &["-i", GT1, "-f", FQ1, "-g", &rel("good.txt")], Some(2));
}
