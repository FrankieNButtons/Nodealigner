use clap::ArgMatches;
use gfa_reader::Gfa;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Arc, Mutex};

/// Extract paths and node coordinates from GFA
/// Writes a TSV with columns: node, start, end, seq, length, path (supports P & W by converting W to paths).
pub fn extract_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let gfa_file = matches.get_one::<String>("gfa").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let threads = matches
        .get_one::<String>("threads")
        .unwrap_or(&String::new())
        .clone();
    let num_threads: usize = threads.parse().unwrap_or(1);

    println!("[info] [extract] Running with arguments:");
    println!("    --gfa     : {gfa_file}");
    println!("    --output  : {output_file}");
    println!("    --threads : {threads}");

    if num_threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .ok();
        println!("[info] Rayon thread pool set to {num_threads} threads");
    }

    println!("[info] Reading GFA and extracting path information...");
    let mut graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file_multi(gfa_file, num_threads);

    // Convert GFA2 Walks (W) into regular Paths so we handle both P & W uniformly.
    // The separator synthesizes a PanSN-style path name from (sample,hap,seq,coords).
    // Overlaps are set to "*" per the walk specification.
    let converted_before = graph.walk.len();
    graph.walk_to_path("#");
    if converted_before > 0 {
        println!("[info] Converted {converted_before} walks (W) into path entries");
    } else {
        println!("[info] No walks (W) found to convert; proceeding with native P paths only");
    }

    let out_file = File::create(output_file)?;
    let out = Arc::new(Mutex::new(BufWriter::new(out_file)));
    {
        let mut guard = out.lock().unwrap();
        writeln!(guard, "node\tstart\tend\tseq\tlength\tpath")?;
    }
    // Parallel, streamed write: per-path local buffer -> append under lock; order not guaranteed
    graph.paths.par_iter().for_each(|path| {
        let mut start = 0usize;
        let mut local_buf = String::with_capacity(path.nodes.len().saturating_mul(32));
        for node in &path.nodes {
            let seq = graph.get_sequence_by_id(node);
            let len = seq.len();
            let end = start + len;
            use std::fmt::Write as _;
            let _ = writeln!(
                &mut local_buf,
                "{}\t{}\t{}\t{}\t{}\t{}",
                node, start, end, seq, len, path.name
            );
            start = end;
        }
        if !local_buf.is_empty() {
            let mut guard = out.lock().unwrap();
            let _ = guard.write_all(local_buf.as_bytes());
        }
    });
    println!("[info] Extraction complete. Output written to {output_file}.");
    Ok(())
}
