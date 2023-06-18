use clap::Parser;
use hmmibd_rs::args::Arguments;
use hmmibd_rs::data::{InputData, OutputFiles};
use hmmibd_rs::hmm::HmmRunner;
fn main() {
    let cli = Arguments::parse();

    let mut out = OutputFiles::new_from_args(&cli);
    let input = InputData::from_args(cli);

    let runner = HmmRunner::new(&input);

    for pair in input.pairs.iter() {
        let pair = (pair.0 as usize, pair.1 as usize);
        runner.run_hmm_on_pair(pair, &mut out);
        // println!("finish pair: {:4}, {:4}", pair.0, pair.1);
    }
}
