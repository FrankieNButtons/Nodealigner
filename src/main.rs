mod extract;
mod io_stream;
mod sort_main;

use clap::{Arg, Command};
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use std::io::{BufRead, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let app = Command::new("gfa2bin-aligner")
        .version("0.2.0")
        .about("Graph VCF toolkit: combine and extract modes")
        .subcommand(
            Command::new("combine")
                .about("Combine VCF with alignment TSV, replacing #CHROM by path, with filter/sort/threads. Optionally use reference.tsv as fallback.")
                .arg(Arg::new("vcf").short('v').long("vcf").help("Input VCF file").required(true))
                .arg(Arg::new("alignment").short('a').long("alignment").help("TSV file (alignment)").required(true))
                .arg(Arg::new("reference").short('r').long("reference").help("Optional reference.tsv for extra #CHROM mapping (produced by extract)").num_args(1))
                .arg(Arg::new("skip").short('s').long("skip").help("Comma-separated substrings. A record is dropped if its raw #CHROM contains any of them.").num_args(1))
                .arg(Arg::new("ignore").long("ignore").help("Ignore/normalize CHROM level [0-5] (applied after --skip): 0=keep, 1=has 'chr', 2=token [0-9XYM], 3=no suffix, 4=only chr{1..22,X,Y,M}, 5=only {1..22,X,Y,M}").num_args(1).default_value("4"))
                .arg(Arg::new("output").short('o').long("output").help("Output VCF file path (default: <input>.replaced.vcf)"))
                .arg(Arg::new("sort").long("sort").help("Sort VCF records (default by POS ascending)").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("prefix").short('p').long("prefix").help("Column to sort by: keyword (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) or 0-based index").default_value("POS"))
                .arg(Arg::new("reverse").long("reverse").help("Sort descending (big to small)").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("threads").short('T').long("threads").help("Number of threads for Rayon").num_args(1))
        )
        .subcommand(
            Command::new("extract")
                .about("Extract all paths from a GFA into a four-column TSV: node, start, end, path")
                .arg(Arg::new("gfa").short('g').long("gfa").help("Input GFA file").required(true))
                .arg(Arg::new("output").short('o').long("output").help("Output TSV file (4 columns: node, start, end, path)").default_value("./reference.tsv"))
                .arg(Arg::new("threads").short('T').long("threads").help("Number of threads for parsing GFA").num_args(1))
        )
        .subcommand(
            Command::new("sort")
                .about("Sort a VCF file by the specified column index")
                .arg(Arg::new("vcf").short('v').long("vcf").help("Input VCF file").required(true))
                .arg(Arg::new("prefix").short('p').long("prefix").help("Column to sort by: keyword (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) or 0-based index").default_value("POS"))
                .arg(Arg::new("reverse").long("reverse").help("Sort descending").action(clap::ArgAction::SetTrue))
        );
    let matches = app.get_matches();

    match matches.subcommand() {
        Some(("combine", sub_m)) => combine_main(sub_m),
        Some(("extract", sub_m)) => extract::extract_main(sub_m),
        Some(("sort", sub_m)) => sort_main::sort_main(sub_m),
        _ => {
            println!("Please choose a subcommand: combine or extract. Use --help for details.");
            Ok(())
        }
    }
}

fn combine_main(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let vcf_path = matches.get_one::<String>("vcf").unwrap();
    let tsv_path = matches.get_one::<String>("alignment").expect("--alignment/-a is required");

    // helper: insert ".sorted" before trailing ".vcf"; if no .vcf, append ".sorted.vcf"
    fn with_sorted_suffix(p: &str) -> String {
        if let Some(stripped) = p.strip_suffix(".vcf") {
            format!("{}.sorted.vcf", stripped)
        } else {
            format!("{}.sorted.vcf", p)
        }
    }

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
    // final output name rule: if we sorted, insert ".sorted" before ".vcf"
    // (combine without --sort keeps the original base name)
    // NOTE: we compute this early only for logging; actual file writes below use this too
    let mut final_output_path = output_path.clone();

    let skip_keywords = matches.get_one::<String>("skip").unwrap_or(&String::new()).clone();
    let ignore_level: u8 = matches
        .get_one::<String>("ignore")
        .and_then(|s| s.parse::<u8>().ok())
        .unwrap_or(0);
    let sort_enabled = matches.get_flag("sort");
    if sort_enabled { final_output_path = with_sorted_suffix(&output_path); }
    let threads = matches.get_one::<String>("threads").unwrap_or(&String::new()).clone();
    let prefix_key = matches.get_one::<String>("prefix").map(|s| s.to_string()).unwrap_or_else(|| "POS".to_string());
    let reverse = matches.get_flag("reverse");

    let reference_path = matches.get_one::<String>("reference").map(|s| s.as_str());

    println!("[info] [combine] Running with arguments:");
    println!("    --vcf      : {vcf_path}");
    println!("    --alignment: {tsv_path}");
    println!("    --reference: {:?}", reference_path);
    println!("    --output   : {output_path}{}", if sort_enabled { "  -> final: ".to_string() + &final_output_path } else { String::new() });
    println!("    --skip     : {skip_keywords}");
    println!("    --ignore   : {ignore_level}");
    println!("    --sort     : {sort_enabled}");
    println!("    --threads  : {threads}");
    println!("    --prefix   : {prefix_key}");
    println!("    --reverse  : {reverse}");
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
    let mut node2path = io_stream::read_tsv_node_path(tsv_path)?;
    // 这里会自动打印debug: node2path contains key 1? true/false

    if let Some(ref_path) = reference_path {
        println!("[info] Reading reference.tsv: {ref_path}");
        let ref_map = io_stream::read_reference_tsv(ref_path)?;
        let mut count = 0;
        for (node, path) in ref_map {
            if !node2path.contains_key(&node) {
                node2path.insert(node, path);
                count += 1;
            }
        }
        println!("[info] reference.tsv merged, {} new node-paths added.", count);
    }

    // --- Streaming pass to temp file ---
    let tmp_out = format!("{output_path}.tmp");
    println!("[info] Streaming CHROM replacement & CHROM-skip to temp: {}", tmp_out);
    let stats = io_stream::stream_replace_chrom_to_tmp(
        vcf_path,
        &tmp_out,
        &node2path,
        &skip_keywords_set,
        ignore_level,
    )?;
    println!(
        "[info] Streaming complete: total={}, replaced={}, skipped={}, unmapped={}",
        stats.total, stats.replaced, stats.skipped, stats.unmapped
    );

    // --- Sort or finalize ---
    if sort_enabled {
        // Determine column index from prefix keyword or numeric index
        let col_index: usize = match prefix_key.as_str() {
            "CHROM" | "chrom" | "#CHROM" | "#chrom" => 0,
            "POS" | "pos" => 1,
            "ID" | "id" => 2,
            "REF" | "ref" => 3,
            "ALT" | "alt" => 4,
            "QUAL" | "qual" => 5,
            "FILTER" | "filter" => 6,
            "INFO" | "info" => 7,
            "FORMAT" | "format" => 8,
            other => other.parse::<usize>().unwrap_or(1), // fallback to POS
        };
        println!("[info] Sorting temp VCF by column '{}' -> index {}{}", prefix_key, col_index, if reverse { " (reverse)" } else { "" });

        // Read temp file, split header vs data
        let file = std::fs::File::open(&tmp_out)?;
        let reader = std::io::BufReader::new(file);
        let mut header_lines: Vec<String> = Vec::new();
        let mut data_lines: Vec<String> = Vec::new();
        for line in reader.lines() {
            let l = line?;
            if l.starts_with('#') { header_lines.push(l); } else { data_lines.push(l); }
        }

        // Sort data lines by the chosen column
        data_lines.sort_by(|a, b| {
            let a_fields: Vec<&str> = a.split('\t').collect();
            let b_fields: Vec<&str> = b.split('\t').collect();
            let a_key = a_fields.get(col_index).copied().unwrap_or("");
            let b_key = b_fields.get(col_index).copied().unwrap_or("");
            let ord = match (a_key.parse::<f64>(), b_key.parse::<f64>()) {
                (Ok(na), Ok(nb)) => na.partial_cmp(&nb).unwrap_or(std::cmp::Ordering::Equal),
                _ => a_key.cmp(b_key),
            };
            if reverse { ord.reverse() } else { ord }
        });

        // Write headers + sorted records to output
        let mut out = std::io::BufWriter::new(std::fs::File::create(&final_output_path)?);
        for h in header_lines { writeln!(out, "{}", h)?; }
        for d in data_lines { writeln!(out, "{}", d)?; }
        println!("[info] Sorting done: wrote {} records to {}", stats.total - stats.skipped, final_output_path);

        // Remove temp
        let _ = fs::remove_file(&tmp_out);
    } else {
        fs::rename(&tmp_out, &output_path)?;
    }

    println!("[info] All operations complete. Output written to {output_path}.");

    Ok(())
}