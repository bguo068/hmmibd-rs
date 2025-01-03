use anyhow::{anyhow, Result};
use clap::Parser;
use hmmibd_rs::args::Arguments;
use hmmibd_rs::data::{self, InputData, OutputBuffer, OutputFiles};
use hmmibd_rs::hmm::HmmRunner;
use rayon::prelude::*;

fn main() -> Result<()> {
    let cli = Arguments::parse();
    let num_threads = cli.num_threads;
    let par_chunk_size = cli.par_chunk_size;
    let suppress_frac = cli.suppress_frac;
    let buffer_size_segments = cli.buffer_size_segments;
    let buffer_size_frac = cli.buffer_size_frac;

    if cli.data_file1 == "-" {
        if !cli.from_bcf {
            return Err(data::Error::CliArgError(
                "Reading from stdin is currently only supported when --from-bcf is used",
            )
            .into());
        }

        if cli.output.is_none() {
            return Err(data::Error::CliArgError(
                "When reading from stdin is specified, --output should be specified as it cannot be inferred",
            ).into());
        }
    }

    let input = InputData::from_args(&cli)?;
    let outfiles = OutputFiles::new_from_args(&cli, buffer_size_segments, buffer_size_frac)?;

    eprintln!("{:#?}", &input.args);

    // use local theadpool instead of global threadpool
    // so that the number of threads used in this instance
    // is not affecting other instances of the same program run
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()?;

    if cli.par_mode == 0 {
        let runner = HmmRunner::new(&input);

        pool.install(|| -> Result<()> {
            input
                .pairs
                .par_chunks(par_chunk_size as usize)
                .try_for_each(|chunk| -> Result<()> {
                    // a and b are reused for all pair iterations.  a and b are
                    // in a larger scope to avoid repeated calculations which
                    // reduce run time by > 20%

                    // transition prob matrix, nsites x 2 x2
                    let mut a = vec![];
                    // emission prob matrix,  nsites x  x 2 (O is known)
                    let mut b = vec![];

                    for pair in chunk.iter() {
                        let pair = (pair.0 as usize, pair.1 as usize);
                        let mut out = OutputBuffer::new(&outfiles, 5, 1);
                        runner.run_hmm_on_pair(pair, &mut a, &mut b, &mut out, suppress_frac)?;
                        out.flush_frac()?;
                        out.flush_segs()?;
                        // println!("finish pair: {:4}, {:4}", pair.0, pair.1);
                    }
                    Ok(())
                })?;
            Ok(())
        })?
    } else if cli.par_mode == 1 {
        use std::sync::{Arc, RwLock};
        let chunks_done = Arc::new(RwLock::new((0usize, 0usize)));
        let chunk_pairs = input.get_chunk_pairs();
        let chunks_all = chunk_pairs.len();

        pool.install(|| -> anyhow::Result<()> {
            chunk_pairs
                .par_iter()
                .try_for_each(|chunkpair| -> anyhow::Result<()> {
                    // a and b are reused for all pair iterations.  a and b are
                    // in a larger scope to avoid repeated calculations which
                    // reduce run time by > 20%

                    // transition prob matrix, nsites x 2 x2
                    let mut a = vec![];
                    // emission prob matrix,  nsites x  x 2 (O is known)
                    let mut b = vec![];

                    let chunkpair_input = input.clone_inputdata_for_chunkpair(*chunkpair)?;
                    let max_npairs = par_chunk_size as usize * par_chunk_size as usize - 1;
                    let mut outbuffer = OutputBuffer::new(&outfiles, 10 * max_npairs, max_npairs);
                    let pairs = &chunkpair_input.pairs;
                    let runner = HmmRunner::new(&chunkpair_input);

                    for pair in pairs.iter() {
                        let pair = (pair.0 as usize, pair.1 as usize);
                        runner.run_hmm_on_pair(
                            pair,
                            &mut a,
                            &mut b,
                            &mut outbuffer,
                            suppress_frac,
                        )?;
                    }

                    outbuffer.flush_frac()?;
                    outbuffer.flush_segs()?;
                    {
                        let (last, n) = {
                            *chunks_done
                                .read()
                                .map_err(|e| anyhow!("cannot read rwlock with error : {e:#?}"))?
                        };

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
                        chunks_done
                            .write()
                            .map_err(|e| anyhow!("cannot write due to RwLock with error: {e:#?} "))?
                            .1 += 1;
                    };
                    Ok(())
                })?;
            Ok(())
        })?;
    } else {
        eprintln!(
            "--par-mode {} is not supported. Please check `hmmibd-rs --help` for help info.",
            cli.par_mode
        );
        std::process::exit(-1);
    }
    Ok(())
}
