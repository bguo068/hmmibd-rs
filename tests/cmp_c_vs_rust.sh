#! /usr/bin/env bash
# Compare hmmibd-rs against the original C hmmIBD (c/hmmIBD.c) on the bundled data.
# Both output files are compared for every case:
#   .hmm.txt        exact lines, order-insensitive
#   .hmm_fract.txt  field-wise after numeric normalisation (C prints -2.54409e+03, Rust -2.54409e3)
# Usage (any directory):  tests/cmp_c_vs_rust.sh          (fast cases, ~1 min)
#                         FULL=1 tests/cmp_c_vs_rust.sh   (also the full sim_data set; the C program takes ~30 min)
# Exit status is 0 only if every comparison passes.
#
# The C program has fixed-size filename buffers (output prefix + extension must fit in 64 bytes,
# input paths in 128), so everything runs from the repository root with short relative paths and
# temporary files are kept in target/cvr_sh/.
set -u
cd "$(dirname "$0")/.."
W=target/cvr_sh; rm -rf "$W"; mkdir -p "$W"; trap 'rm -rf "$W"' EXIT
nfail=0

# Compile / build
${CC:-cc} c/hmmIBD.c -lm -O2 -o "$W/hmmIBD" 2>/dev/null > /dev/null || { echo "cannot compile c/hmmIBD.c"; exit 1; }
cargo build -q --release --bin hmmibd-rs > /dev/null || { echo "cargo build failed"; exit 1; }
RS=target/release/hmmibd-rs; C=$W/hmmIBD

# normalise a frac file: drop header, re-format every numeric column through awk so formatting differences vanish, sort
norm_frac() {
  sed 1d "$1" | awk -F'\t' 'BEGIN{OFS="\t"} {printf "%s\t%s\t%d\t%.4f\t%.5e\t%d\t%.3f\t%d\t%.5f\t%.5f\t%.5f\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | sort
}
norm_seg() { sed 1d "$1" | sort; }

# compare <name> <expected_rows|-> <hmmIBD args...>   (runs both programs with the same args)
compare() {
  local name=$1 want=$2; shift 2
  "$C"  "$@" -o "$W/${name}_c"  > /dev/null 2>&1 || { echo "FAIL $name: hmmIBD run failed";    nfail=$((nfail+1)); return; }
  "$RS" "$@" -o "$W/${name}_rs" --num-threads 1 > /dev/null 2>&1 || { echo "FAIL $name: hmmibd-rs run failed"; nfail=$((nfail+1)); return; }
  local dseg dfrac rows
  dseg=$(diff <(norm_seg "$W/${name}_c.hmm.txt") <(norm_seg "$W/${name}_rs.hmm.txt") | wc -l)
  dfrac=$(diff <(norm_frac "$W/${name}_c.hmm_fract.txt") <(norm_frac "$W/${name}_rs.hmm_fract.txt") | wc -l)
  rows=$(sed 1d "$W/${name}_rs.hmm_fract.txt" | wc -l | tr -d ' ')
  if [ "$dseg" -eq 0 ] && [ "$dfrac" -eq 0 ] && { [ "$want" = "-" ] || [ "$rows" -eq "$want" ]; }; then
    echo "pass $name ($rows pairs)"
  else
    echo "FAIL $name: seg diff lines=$dseg frac diff lines=$dfrac pairs=$rows (expected $want)"; nfail=$((nfail+1))
  fi
}

D=c/samp_data
# Original data provided in c/samp_data
compare samp_data_single_pop        45  -i $D/pf3k_Cambodia_13.txt -f $D/freqs_pf3k_Cambodia_13.txt -r 6.66667e-7
compare samp_data_single_pop_nofreq 45  -i $D/pf3k_Cambodia_13.txt
compare samp_data_two_pops          100 -i $D/pf3k_Cambodia_13.txt -I $D/pf3k_Ghana_13.txt -f $D/freqs_pf3k_Cambodia_13.txt -F $D/freqs_pf3k_Ghana_13.txt -r 6.66667e-7
compare samp_data_max_iter_1        45  -i $D/pf3k_Cambodia_13.txt -f $D/freqs_pf3k_Cambodia_13.txt -m 1
# bad-sample file, one id per line (documented format): both must drop the same two samples
printf 'PH0047-Cx\nPH0051-Cx\n' > "$W/bad.txt"
compare samp_data_bad_samples       28  -i $D/pf3k_Cambodia_13.txt -f $D/freqs_pf3k_Cambodia_13.txt -b "$W/bad.txt"
# good-pair file (pairs in header order; a reversed pair is run by both but printed swapped by hmmibd-rs)
printf 'PH0047-Cx\tPH0051-Cx\nPH0055-C\tPH0057-C\n' > "$W/good.txt"
compare samp_data_good_pairs        2   -i $D/pf3k_Cambodia_13.txt -f $D/freqs_pf3k_Cambodia_13.txt -g "$W/good.txt"
# 40 variants: no pair reaches 10 informative sites, the pair filter must remove all of them in both programs
head -41 $D/pf3k_Cambodia_13.txt > "$W/gt40.txt"; head -40 $D/freqs_pf3k_Cambodia_13.txt > "$W/fq40.txt"
compare samp_data_few_variants      0   -i "$W/gt40.txt" -f "$W/fq40.txt"

# --par-mode 1 (chunked samples) must give the same result as the default mode
"$RS" -i $D/pf3k_Cambodia_13.txt -f $D/freqs_pf3k_Cambodia_13.txt -r 6.66667e-7 --par-mode 1 -o "$W/parmode1" > /dev/null 2>&1
if [ "$(diff <(norm_seg "$W/samp_data_single_pop_c.hmm.txt") <(norm_seg "$W/parmode1.hmm.txt") | wc -l)" -eq 0 ] && \
   [ "$(diff <(norm_frac "$W/samp_data_single_pop_c.hmm_fract.txt") <(norm_frac "$W/parmode1.hmm_fract.txt") | wc -l)" -eq 0 ]; then
  echo "pass samp_data_single_pop_parmode1"; else echo "FAIL samp_data_single_pop_parmode1"; nfail=$((nfail+1)); fi

# Simulated data, subset (first 98 samples, three chromosomes)
S=testdata/sim_data
cut -f 1-100 $S/gt_chr1.txt          >  "$W/gt.txt"
cut -f 1-100 $S/gt_chr2.txt | sed 1d >> "$W/gt.txt"
cut -f 1-100 $S/gt_chr3.txt | sed 1d >> "$W/gt.txt"
cat $S/frq_chr{1,2,3}.txt           >  "$W/frq.txt"
compare sim_data_subset             -   -i "$W/gt.txt" -f "$W/frq.txt" -r 6.66667e-7

if [ "${FULL:-0}" = 1 ]; then
  cat $S/gt_chr1.txt          >  "$W/gt_full.txt"
  cat $S/gt_chr2.txt | sed 1d >> "$W/gt_full.txt"
  cat $S/gt_chr3.txt | sed 1d >> "$W/gt_full.txt"
  compare sim_data_full             -   -i "$W/gt_full.txt" -f "$W/frq.txt" -r 6.66667e-7   # ~30 min for the C program
fi

echo; if [ "$nfail" -eq 0 ]; then echo "all comparisons passed"; else echo "$nfail comparison(s) FAILED"; fi
[ "$nfail" -eq 0 ]
