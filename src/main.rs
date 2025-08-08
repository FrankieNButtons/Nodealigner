mod extract;
mod filter;
mod io_utils;

use clap::{Arg, Command};
use std::collections::HashSet;
use rayon::prelude::*;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let app = Command::new("gfa2bin-aligner")
        .version("0.2.0")
        .about("Graph VCF toolkit: combine and extract modes")
        .subcommand(
            Command::new("combine")
                .about("Combine VCF with alignment TSV, replacing #CHROM by path, with filter/sort/threads. Optionally use reference.tsv as fallback.")
                .arg(Arg::new("vcf").short('v').long("vcf").help("Input VCF file").required(true))
                .arg(Arg::new("tsv").short('t').long("tsv").help("TSV file (alignment)").required(true))
                .arg(Arg::new("reference").short('r').long("reference").help("Optional reference.tsv for extra #CHROM mapping (produced by extract)").num_args(1))
                .arg(Arg::new("skip").short('s').long("skip").help("Comma-separated path keywords to skip.").num_args(1))
                .arg(Arg::new("output").short('o').long("output").help("Output VCF file path (default: <input>.replaced.vcf)"))
                .arg(Arg::new("sort").long("sort").help("Sort VCF records by POS (ascending)").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("threads").short('T').long("threads").help("Number of threads for Rayon").num_args(1))
        )
        .subcommand(
            Command::new("extract")
                .about("Extract all paths from a GFA into a four-column TSV: node, start, end, path")
                .arg(Arg::new("gfa").short('g').long("gfa").help("Input GFA file").required(true))
                .arg(Arg::new("output").short('o').long("output").help("Output TSV file (4 columns: node, start, end, path)").default_value("./reference.tsv"))
                .arg(Arg::new("threads").short('T').long("threads").help("Number of threads for parsing GFA").num_args(1))
        );
    let matches = app.get_matches();

    match matches.subcommand() {
        Some(("combine", sub_m)) => combine_main(sub_m),
        Some(("extract", sub_m)) => extract::extract_main(sub_m),
        _ => {
            println!("Please choose a subcommand: combine or extract. Use --help for details.");
            Ok(())
        }
    }
}

fn combine_main(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let vcf_path = matches.get_one::<String>("vcf").unwrap();
    let tsv_path = matches.get_one::<String>("tsv").unwrap();

    let default_output = {
        let p = Path::new(vcf_path);
        if let Some(stem) = p.file_stem().and_then(|s| s.to_str()) {
            if let Some(parent) = p.parent() {
                parent.join(format!("{stem}.replaced.vcf")).to_str().unwrap().to_owned()
            } else {
                format!("{stem}.replaced.vcf")
            }
        } else {
            "output.vcf".to_owned()
        }
    };
    let output_path = matches.get_one::<String>("output").map(|s| s.to_owned()).unwrap_or(default_output);

    let skip_keywords = matches.get_one::<String>("skip").unwrap_or(&String::new()).clone();
    let sort_enabled = matches.get_flag("sort");
    let threads = matches.get_one::<String>("threads").unwrap_or(&String::new()).clone();

    let reference_path = matches.get_one::<String>("reference").map(|s| s.as_str());

    println!("[info] [combine] Running with arguments:");
    println!("    --vcf      : {vcf_path}");
    println!("    --tsv      : {tsv_path}");
    println!("    --reference: {:?}", reference_path);
    println!("    --output   : {output_path}");
    println!("    --skip     : {skip_keywords}");
    println!("    --sort     : {sort_enabled}");
    println!("    --threads  : {threads}");
    println!("    #CHROM field will always be replaced by the corresponding path from TSV, or reference.tsv as fallback if provided.");

    let num_threads: Option<usize> = matches.get_one::<String>("threads")
        .and_then(|s| s.parse().ok());
    if let Some(n) = num_threads {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().unwrap();
        println!("[info] Rayon thread pool set to {n} threads");
    }

    let skip_keywords_set: HashSet<String> = matches.get_one::<String>("skip")
        .map(|s| s.split(',').map(|x| x.trim().to_string()).filter(|s| !s.is_empty()).collect())
        .unwrap_or_default();

    println!("[info] Reading TSV and building node-to-path map...");
    let mut node2path = io_utils::read_tsv_node_path(tsv_path)?;
    // 这里会自动打印debug: node2path contains key 1? true/false

    if let Some(ref_path) = reference_path {
        println!("[info] Reading reference.tsv: {ref_path}");
        let ref_map = io_utils::read_reference_tsv(ref_path)?;
        let mut count = 0;
        for (node, path) in ref_map {
            if !node2path.contains_key(&node) {
                node2path.insert(node, path);
                count += 1;
            }
        }
        println!("[info] reference.tsv merged, {} new node-paths added.", count);
    }

    println!("[info] Reading VCF file...");
    let (headers, records) = io_utils::read_vcf_records(vcf_path)?;
    println!("[info] VCF read complete; loaded {} records.", records.len());

    println!("[info] Filtering VCF records...");
    let mut records = filter::filter_records_by_path(records, &node2path, &skip_keywords_set);
    println!("[info] Filtering complete; {} records remaining.", records.len());

    if sort_enabled {
        println!("[info] Sorting VCF records by POS...");
        records.par_sort_by_key(|rec| rec.pos);
        println!("[info] Sorting complete.");
    }

    println!("[info] Replacing #CHROM with path from TSV/reference...");
    for rec in &mut records {
        if rec.pos == 1 {
            // println!("[debug] current record POS == 1, node2path contains key 1? {}", node2path.contains_key(&1));
        }
        if let Some(path) = node2path.get(&rec.pos) {
            if rec.pos == 1 {
                // println!("[debug] FOUND mapping for node=1, path={}", path);
            }
            let mut fields: Vec<&str> = rec.raw.split('\t').collect();
            if fields.len() > 1 {
                fields[0] = path;
                rec.raw = fields.join("\t");
            }
        }
    }
    println!("[info] Replacement of #CHROM complete.");

    println!("[info] Writing output VCF...");
    io_utils::write_vcf_records(&output_path, &headers, &records)?;
    println!("[info] All operations complete. Output written to {output_path}.");

    Ok(())
}