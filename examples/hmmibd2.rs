use clap::Parser;
use hmmibd_rs::args::Arguments;
use hmmibd_rs::data::{InputData, OutputFiles, OutputBuffer};
use hmmibd_rs::hmm::HmmRunner;
use rayon::prelude::*;
fn main() {
    let cli = Arguments::parse();

    let out = OutputFiles::new_from_args(&cli);
    let input = InputData::from_args(cli);

    let runner = HmmRunner::new(&input);

    input.pairs.par_chunks(120).for_each(|chunk|{
        for pair in chunk.iter() {
        let pair = (pair.0 as usize, pair.1 as usize);
        let mut out = OutputBuffer::new(&out, 1000000, 10000);
        runner.run_hmm_on_pair(pair, &mut out);
        out.flush_frac();
        out.flush_segs();
        // println!("finish pair: {:4}, {:4}", pair.0, pair.1);
    }

    });
    

}
