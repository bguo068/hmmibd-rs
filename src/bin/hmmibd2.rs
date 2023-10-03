use clap::Parser;
use hmmibd_rs::args::Arguments;
use hmmibd_rs::data::{InputData, OutputBuffer, OutputFiles};
use hmmibd_rs::hmm::HmmRunner;
use rayon::prelude::*;
fn main() {
    let cli = Arguments::parse();
    let num_threads = cli.num_threads;
    let par_chunk_size = cli.par_chunk_size;

    let out = OutputFiles::new_from_args(&cli);
    let input = InputData::from_args(&cli);
    // println!("{:?}", &input.sites);
    eprintln!("{:?}", &input.args);

    // use local theadpool instead of global threadpool
    // so that the number of threads used in this instance
    // is not affecting other instances of the same program run
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    if cli.par_mode == 0 {
        let runner = HmmRunner::new(&input);

        pool.install(|| {
            input
                .pairs
                .par_chunks(par_chunk_size as usize)
                .for_each(|chunk| {
                    for pair in chunk.iter() {
                        let pair = (pair.0 as usize, pair.1 as usize);
                        let mut out = OutputBuffer::new(&out, 1000000, 10000);
                        runner.run_hmm_on_pair(pair, &mut out);
                        out.flush_frac();
                        out.flush_segs();
                        // println!("finish pair: {:4}, {:4}", pair.0, pair.1);
                    }
                });
        })
    } else if cli.par_mode == 1 {
        eprintln!("par-mode == 1. WARN: still need to verify accuracy of IBD segments");
        use std::sync::{Arc, RwLock};
        let chunks_done = Arc::new(RwLock::new((0usize, 0usize)));
        let chunk_pairs = input.get_chunk_pairs();
        let chunks_all = chunk_pairs.len();

        pool.install(|| {
            chunk_pairs.par_iter().for_each(|chunkpair| {
                let mut out = OutputBuffer::new(&out, 1000000, 10000);
                let chunkpair_input = input.clone_inputdata_for_chunkpair(*chunkpair);
                let runner = HmmRunner::new(&chunkpair_input);
                for pair in chunkpair_input.pairs.iter() {
                    let pair = (pair.0 as usize, pair.1 as usize);
                    runner.run_hmm_on_pair(pair, &mut out);
                }

                out.flush_frac();
                out.flush_segs();
                {
                    let (last, n) = { *chunks_done.read().unwrap() };

                    if n - last > chunks_all / 100 {
                        println!(
                            "finished {n} / {chunks_all}: {:.3}%       ",
                            (n * 100) as f64 / chunks_all as f64
                        );
                        let chunks_done = &mut *chunks_done.write().unwrap();
                        chunks_done.0 = chunks_done.1;
                    }
                }
                {
                    chunks_done.write().unwrap().1 += 1;
                }
            })
        })
    } else {
        eprintln!(
            "--par-mode {} is not supported. Please check `hmmibd2 --help` for help info.",
            cli.par_mode
        );
        std::process::exit(-1);
    }
}
