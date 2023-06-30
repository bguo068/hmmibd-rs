#! /usr/bin/env bash
set -e

# Compile

gcc ../hmmIBD.c -lm -O2 -o hmmIBD 2>/dev/null > /dev/null
cargo build -q --release --bin hmmibd2 > /dev/null && cp ../target/release/hmmibd2 ./



# Compare c vs rust versions using original data provided in samp_data folder

./hmmIBD -i ../samp_data/pf3k_Cambodia_13.txt -I ../samp_data/pf3k_Ghana_13.txt \
    -f ../samp_data/freqs_pf3k_Cambodia_13.txt  -F ../samp_data/freqs_pf3k_Ghana_13.txt \
    -o orig -r 6.66667e-7 > /dev/null
./hmmibd2 -i ../samp_data/pf3k_Cambodia_13.txt -I ../samp_data/pf3k_Ghana_13.txt \
    -f ../samp_data/freqs_pf3k_Cambodia_13.txt  -F ../samp_data/freqs_pf3k_Ghana_13.txt \
    -o new  -r 6.66667e-7 > /dev/null
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass samp_data; 
else echo failed samp_data ; fi



# Compare c vs rust versions using simulated data

# use a subset of data
cut -f 1-100 ../sim_data/gt_chr1.txt          >  gt.txt
cut -f 1-100 ../sim_data/gt_chr2.txt | sed 1d >> gt.txt
cut -f 1-100 ../sim_data/gt_chr3.txt | sed 1d >> gt.txt
cat ../sim_data/frq_chr{1,2,3}.txt           > frq.txt
./hmmIBD  -i gt.txt -f frq.txt -r 6.66667e-7 -o orig > /dev/null
./hmmibd2 -i gt.txt -f frq.txt -r 6.66667e-7 -o new  > /dev/null
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass sim_data subset; 
else echo failed sim_data subset ; fi

# use full data sets
cat ../sim_data/gt_chr1.txt          >  gt.txt
cat ../sim_data/gt_chr2.txt | sed 1d >> gt.txt
cat ../sim_data/gt_chr3.txt | sed 1d >> gt.txt
cat ../sim_data/frq_chr{1,2,3}.txt           > frq.txt
time ./hmmIBD  -i gt.txt -f frq.txt -r 6.66667e-7 -o orig > /dev/null # 30.0min
time ./hmmibd2 -i gt.txt -f frq.txt -r 6.66667e-7 -o new --max-all 2 > /dev/null # 4.5 min
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass sim_data full; 
else echo failed sim_data full ; fi