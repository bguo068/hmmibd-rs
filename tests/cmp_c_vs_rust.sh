#! /usr/bin/env bash
set -e

# Compile

gcc ../c/hmmIBD.c -lm -O2 -o hmmIBD 2>/dev/null > /dev/null
cargo build -q --release --bin hmmibd-rs > /dev/null && cp ../target/release/hmmibd-rs ./

# Compare c vs rust versions using original data provided in c/samp_data folder

./hmmIBD -i ../c/samp_data/pf3k_Cambodia_13.txt  \
    -f ../c/samp_data/freqs_pf3k_Cambodia_13.txt  \
    -o orig -r 6.66667e-7 > /dev/null
./hmmibd-rs -i ../c/samp_data/pf3k_Cambodia_13.txt  \
    -f ../c/samp_data/freqs_pf3k_Cambodia_13.txt  \
    -o new  -r 6.66667e-7 > /dev/null
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass c/samp_data_single_pop; 
else echo failed c/samp_data_c/samp_data_single_pop ; fi


./hmmIBD -i ../c/samp_data/pf3k_Cambodia_13.txt  \
    -f ../c/samp_data/freqs_pf3k_Cambodia_13.txt  \
    -o orig -r 6.66667e-7 > /dev/null
./hmmibd-rs -i ../c/samp_data/pf3k_Cambodia_13.txt  \
    -f ../c/samp_data/freqs_pf3k_Cambodia_13.txt --par-mode 1 \
    -o new  -r 6.66667e-7 > /dev/null
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass c/samp_data_single_pop_parmode1; 
else echo failed c/samp_data_c/samp_data_single_pop_parmode1 ; fi


# Compare c vs rust versions using original data provided in c/samp_data folder

./hmmIBD -i ../c/samp_data/pf3k_Cambodia_13.txt -I ../c/samp_data/pf3k_Ghana_13.txt \
    -f ../c/samp_data/freqs_pf3k_Cambodia_13.txt  -F ../c/samp_data/freqs_pf3k_Ghana_13.txt \
    -o orig -r 6.66667e-7 > /dev/null
./hmmibd-rs -i ../c/samp_data/pf3k_Cambodia_13.txt -I ../c/samp_data/pf3k_Ghana_13.txt \
    -f ../c/samp_data/freqs_pf3k_Cambodia_13.txt  -F ../c/samp_data/freqs_pf3k_Ghana_13.txt \
    -o new  -r 6.66667e-7 > /dev/null
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass c/samp_data; 
else echo failed c/samp_data ; fi



# Compare c vs rust versions using simulated data

# use a subset of data
cut -f 1-100 ../testdata/sim_data/gt_chr1.txt          >  gt.txt
cut -f 1-100 ../testdata/sim_data/gt_chr2.txt | sed 1d >> gt.txt
cut -f 1-100 ../testdata/sim_data/gt_chr3.txt | sed 1d >> gt.txt
cat ../testdata/sim_data/frq_chr{1,2,3}.txt           > frq.txt
./hmmIBD  -i gt.txt -f frq.txt -r 6.66667e-7 -o orig > /dev/null
./hmmibd-rs -i gt.txt -f frq.txt -r 6.66667e-7 -o new  > /dev/null
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass testdata/sim_data subset; 
else echo failed testdata/sim_data subset ; fi

# use full data sets
cat ../testdata/sim_data/gt_chr1.txt          >  gt.txt
cat ../testdata/sim_data/gt_chr2.txt | sed 1d >> gt.txt
cat ../testdata/sim_data/gt_chr3.txt | sed 1d >> gt.txt
cat ../testdata/sim_data/frq_chr{1,2,3}.txt           > frq.txt
time ./hmmIBD  -i gt.txt -f frq.txt -r 6.66667e-7 -o orig > /dev/null # 30.0min
time ./hmmibd-rs -i gt.txt -f frq.txt -r 6.66667e-7 -o new --max-all 2 > /dev/null # 4.5 min
diff <(sed 1d new.hmm.txt | sort ) <(sed 1d orig.hmm.txt | sort ) > diff.txt
if [ `cat diff.txt | wc -l`  -eq 0 ] ; then echo pass testdata/sim_data full; 
else echo failed testdata/sim_data full ; fi

time ./hmmibd-rs -i ../testdata/sim_data/gt_chr1.txt -r 6.66667e-7 -o chr1 --max-all 2 > /dev/null # 4.5 min
time ./hmmibd-rs -i ../testdata/sim_data/gt_chr2.txt -r 6.66667e-7 -o chr2 --max-all 2 > /dev/null # 4.5 min
time ./hmmibd-rs -i ../testdata/sim_data/gt_chr3.txt -r 6.66667e-7 -o chr3 --max-all 2 > /dev/null # 4.5 min
