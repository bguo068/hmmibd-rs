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
    let input = InputData::from_args(cli);
    // println!("{:?}", &input.sites);
    println!("{:?}", &input.args);

    let runner = HmmRunner::new(&input);

    // use local theadpool instead of global threadpool
    // so that the number of threads used in this instance
    // is not affecting other instances of the same program run
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    pool.install(|| {
        input.pairs.par_chunks(par_chunk_size).for_each(|chunk| {
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
}
