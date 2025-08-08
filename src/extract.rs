use gfa_reader::Gfa;
use std::fs::File;
use std::io::{BufWriter, Write};
use clap::ArgMatches;
use rayon::prelude::*;

/// Extract paths and node coordinates from GFA
/// Writes a TSV with columns: node, start, end, path
pub fn extract_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let gfa_file = matches.get_one::<String>("gfa").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let threads = matches.get_one::<String>("threads").unwrap_or(&String::new()).clone();
    let num_threads: usize = threads.parse().unwrap_or(1);

    println!("[info] [extract] Running with arguments:");
    println!("    --gfa     : {gfa_file}");
    println!("    --output  : {output_file}");
    println!("    --threads : {threads}");

    if num_threads > 1 {
        rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().ok();
        println!("[info] Rayon thread pool set to {num_threads} threads");
    }

    println!("[info] Reading GFA and extracting path information...");
    let graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file_multi(gfa_file, num_threads);

    let mut out = BufWriter::new(File::create(output_file)?);
    writeln!(out, "node\tstart\tend\tpath")?;

    let mut all_lines: Vec<String> = graph.paths.par_iter()
        .flat_map(|path| {
            let mut start = 0usize;
            let mut lines = Vec::with_capacity(path.nodes.len());
            for node in path.nodes.iter() {
                let len = graph.get_sequence_by_id(node).len();
                let end = start + len;
                lines.push(format!("{}\t{}\t{}\t{}", node, start, end, path.name));
                start = end;
            }
            lines
        })
        .collect();

    // Optionally sort output for consistency
    all_lines.sort();
    for line in all_lines {
        writeln!(out, "{line}")?;
    }
    println!("[info] Extraction complete. Output written to {output_file}.");
    Ok(())
}