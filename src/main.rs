mod extract;
mod header;
mod io_stream;
mod sort_main;

use clap::{Arg, Command};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, Write};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let app = Command::new("gfa2bin-aligner")
        .version("0.0.4")
        .about("Graph VCF toolkit extends from `gfa2bin`: align and extract modes. Additionally, when reference.tsv is provided, CHROM is set to path, POS may be replaced by the node’s start coordinate from reference.tsv, and ID may be set to the original POS (implementation depends on io_stream).")
        .subcommand(
            Command::new("align")
                .about("Align VCF with alignment TSV, replacing #CHROM by path, with filter/sort/threads. Optionally use reference.tsv as fallback.")
                .arg(Arg::new("vcf").short('v').long("vcf").help("Input VCF file").required(true))
                .arg(Arg::new("alignment").short('a').long("alignment").help("TSV file (alignment)").required(true))
                .arg(Arg::new("reference").short('r').long("reference").help("reference.tsv for CHROM mapping/header synthesis; required unless --no-header. Supports 4- or 6-column TSV (last column always path); optionally uses third column start to set VCF POS during alignment.").num_args(1).required_unless_present("no-header"))
                .arg(Arg::new("skip").short('s').long("skip").help("Comma-separated substrings. A record is dropped if its raw #CHROM contains any of them.").num_args(1))
                .arg(Arg::new("ignore").long("ignore").help("Ignore/normalize CHROM level [0-5] (applied after --skip): 0=keep, 1=has 'chr', 2=token [0-9XYM], 3=no suffix, 4=only chr{1..22,X,Y,M}, 5=only {1..22,X,Y,M}").num_args(1).default_value("4"))
                .arg(Arg::new("output").short('o').long("output").help("Output VCF file path (default: <input>.replaced.vcf)"))
                .arg(Arg::new("sort").long("sort").help("Sort VCF records (default by POS ascending)").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("prefix").short('p').long("prefix").help("Column to sort by: keyword (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) or 0-based index (default: POS)").default_value("POS"))
                .arg(Arg::new("reverse").long("reverse").help("Sort descending (big to small)").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("threads").short('T').long("threads").help("Number of threads for Rayon").num_args(1))
                .arg(Arg::new("no-header").long("no-header").help("Do not synthesize a header on the combined VCF (by default, header is added using reference.tsv)").action(clap::ArgAction::SetTrue))
        )
        .subcommand(
            Command::new("extract")
                .about("Extract all paths from a GFA into TSV (4 or 6 columns): node, start, end, [seq, length,] path")
                .arg(Arg::new("gfa").short('g').long("gfa").help("Input GFA file").required(true))
                .arg(Arg::new("output").short('o').long("output").help("Output TSV file (4 or 6 columns). Last column is path.").default_value("./reference.tsv"))
                .arg(Arg::new("threads").short('T').long("threads").help("Number of threads for parsing GFA").num_args(1))
        )
        .subcommand(
            Command::new("header")
                .about("Synthesize a VCF header using reference.tsv (columns: node, start, end[, seq, length], path) and keys inferred from the VCF body; writes a new VCF with merged header")
                .arg(Arg::new("vcf").short('v').long("vcf").help("Input VCF/VCF-like file (plain text)").required(true))
                .arg(Arg::new("reference").short('r').long("reference").help("reference.tsv produced by 'extract' (4 or 6 cols: node, start, end[, seq, length], path; last column is path)").required(true))
                .arg(Arg::new("output").short('o').long("output").help("Output VCF path (default: <input>.withheader.vcf)"))
                .arg(Arg::new("threads").short('T').long("threads").help("Number of threads for Rayon").num_args(1))
                .arg(
                    Arg::new("ignore")
                        .long("ignore")
                        .help("Ignore/normalize CHROM level [0-5] (applied after --skip): 0=keep, 1=has 'chr', 2=token [0-9XYM], 3=no suffix, 4=only chr{1..22,X,Y,M}, 5=only {1..22,X,Y,M}")
                        .num_args(1)
                        .default_value("4")
                )
        )
        .subcommand(
            Command::new("sort")
                .about("Sort a VCF file; default is by POS (ascending). Use --prefix to choose column and --reverse for descending.")
                .arg(Arg::new("vcf").short('v').long("vcf").help("Input VCF file").required(true))
                .arg(Arg::new("prefix").short('p').long("prefix").help("Column to sort by: keyword (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) or 0-based index (default: POS)").default_value("POS"))
                .arg(Arg::new("reverse").long("reverse").help("Sort descending").action(clap::ArgAction::SetTrue))
                .arg(Arg::new("output").short('o').long("output").help("Output VCF path (default: <input>.sorted.vcf)"))
        );
    let matches = app.get_matches();

    match matches.subcommand() {
        Some(("align", sub_m)) => align_main(sub_m),
        Some(("extract", sub_m)) => extract::extract_main(sub_m),
        Some(("header", sub_m)) => header::header_main(sub_m),
        Some(("sort", sub_m)) => sort_main::sort_main(sub_m),
        _ => {
            println!(
                "Please choose a subcommand: align, extract, header, or sort. Use --help for details."
            );
            Ok(())
        }
    }
}

fn align_main(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let vcf_path = matches.get_one::<String>("vcf").unwrap();
    let tsv_path = matches
        .get_one::<String>("alignment")
        .expect("--alignment/-a is required");

    // helper: insert ".sorted" before trailing ".vcf"; if no .vcf, append ".sorted.vcf"
    fn with_sorted_suffix(p: &str) -> String {
        if let Some(stripped) = p.strip_suffix(".vcf") {
            format!("{}.sorted.vcf", stripped)
        } else {
            format!("{}.sorted.vcf", p)
        }
    }

    // helper: build a headed output path in the *original VCF's directory*,
    // inserting ".headed" before ".vcf" and preserving the candidate's base name.
    fn headed_in_original_dir(original_vcf: &str, candidate_path: &str) -> String {
        use std::path::Path;
        let orig_dir = Path::new(original_vcf)
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| Path::new(".").to_path_buf());
        let candidate_name = Path::new(candidate_path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("output.vcf");
        let name = if let Some(stripped) = candidate_name.strip_suffix(".vcf") {
            format!("{}.headed.vcf", stripped)
        } else {
            format!("{}.headed.vcf", candidate_name)
        };
        orig_dir.join(name).to_string_lossy().into_owned()
    }

    let default_output = {
        let p = Path::new(vcf_path);
        if let Some(stem) = p.file_stem().and_then(|s| s.to_str()) {
            if let Some(parent) = p.parent() {
                parent
                    .join(format!("{stem}.replaced.vcf"))
                    .to_str()
                    .unwrap()
                    .to_owned()
            } else {
                format!("{stem}.replaced.vcf")
            }
        } else {
            "output.vcf".to_owned()
        }
    };
    let output_path = matches
        .get_one::<String>("output")
        .map(|s| s.to_owned())
        .unwrap_or(default_output);
    // final output name rule: if we sorted, insert ".sorted" before ".vcf"
    // (align without --sort keeps the original base name)
    // NOTE: we compute this early only for logging; actual file writes below use this too
    let mut final_output_path = output_path.clone();

    let skip_keywords = matches
        .get_one::<String>("skip")
        .unwrap_or(&String::new())
        .clone();
    let ignore_level: u8 = matches
        .get_one::<String>("ignore")
        .and_then(|s| s.parse::<u8>().ok())
        .unwrap_or(4);
    let sort_enabled = matches.get_flag("sort");
    if sort_enabled {
        final_output_path = with_sorted_suffix(&output_path);
    }
    let threads = matches
        .get_one::<String>("threads")
        .unwrap_or(&String::new())
        .clone();
    let prefix_key = matches
        .get_one::<String>("prefix")
        .map(|s| s.to_string())
        .unwrap_or_else(|| "POS".to_string());
    let reverse = matches.get_flag("reverse");
    let no_header = matches.get_flag("no-header");

    // Note: we will later pull `start` for POS replacement and swap POS→ID in io_stream
    let reference_path = matches.get_one::<String>("reference").map(|s| s.as_str());

    println!("[info] [align] Running with arguments:");
    println!("    --vcf      : {vcf_path}");
    println!("    --alignment: {tsv_path}");
    println!("    --reference: {:?}", reference_path);
    println!(
        "    --output   : {output_path}{}",
        if sort_enabled {
            "  -> final: ".to_string() + &final_output_path
        } else {
            String::new()
        }
    );
    println!("    --skip     : {skip_keywords}");
    println!("    --ignore   : {ignore_level}");
    println!("    --sort     : {sort_enabled}");
    println!("    --threads  : {threads}");
    println!("    --prefix   : {prefix_key}");
    println!("    --reverse  : {reverse}");
    println!("    --no-header: {}", no_header);
    println!(
        "    #CHROM field will always be replaced by the corresponding path from TSV, or reference.tsv as fallback if provided."
    );

    let num_threads: Option<usize> = matches
        .get_one::<String>("threads")
        .and_then(|s| s.parse().ok());
    if let Some(n) = num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap();
        println!("[info] Rayon thread pool set to {n} threads");
    }

    let skip_keywords_set: HashSet<String> = matches
        .get_one::<String>("skip")
        .map(|s| {
            s.split(',')
                .map(|x| x.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect()
        })
        .unwrap_or_default();

    let mut node2path: HashMap<u64, String> = HashMap::new();

    if let Some(ref_path) = reference_path {
        println!("[info] Reading reference.tsv first: {ref_path}");
        let ref_map = io_stream::read_reference_tsv(ref_path)?; // also seeds REF_NODE2* caches
        let ref_count = ref_map.len();
        node2path = ref_map;
        println!(
            "[info] reference.tsv loaded: {} node-paths",
            ref_count
        );
    } else {
        println!("[info] No reference.tsv provided; will rely on alignment TSV for path mapping");
    }

    println!("[info] Reading alignment TSV and adding missing node-paths: {tsv_path}");
    let aln_map = io_stream::read_tsv_node_path(tsv_path)?;
    let mut add_from_aln = 0usize;
    for (node, path) in aln_map {
        if !node2path.contains_key(&node) {
            node2path.insert(node, path);
            add_from_aln += 1;
        }
    }
    println!(
        "[info] alignment.tsv merged: {} new node-paths added (only for nodes missing in reference)",
        add_from_aln
    );

    // --- Streaming pass to temp file ---
    let tmp_out = format!("{output_path}.tmp");
    println!(
        "[info] Streaming CHROM replacement & CHROM-skip to temp: {}",
        tmp_out
    );
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
        println!(
            "[info] Sorting temp VCF by column '{}' -> index {}{}",
            prefix_key,
            col_index,
            if reverse { " (reverse)" } else { "" }
        );

        // Read temp file, split header vs data
        let file = std::fs::File::open(&tmp_out)?;
        let reader = std::io::BufReader::new(file);
        let mut header_lines: Vec<String> = Vec::new();
        let mut data_lines: Vec<String> = Vec::new();
        for line in reader.lines() {
            let l = line?;
            if l.starts_with('#') {
                header_lines.push(l);
            } else {
                data_lines.push(l);
            }
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
        for h in header_lines {
            writeln!(out, "{}", h)?;
        }
        for d in data_lines {
            writeln!(out, "{}", d)?;
        }
        println!(
            "[info] Sorting done: wrote {} records to {}",
            stats.total - stats.skipped,
            final_output_path
        );

        // Remove temp
        let _ = fs::remove_file(&tmp_out);
    } else {
        fs::rename(&tmp_out, &output_path)?;
    }

    // Decide which file we produced
    let align_out = if sort_enabled {
        final_output_path.clone()
    } else {
        output_path.clone()
    };

    let output_path_log = if !no_header {
        // Reference is required unless --no-header per CLI; unwrap is safe here
        let ref_path = matches
            .get_one::<String>("reference")
            .expect("reference.tsv required unless --no-header");
        let threads_opt: Option<usize> = matches
            .get_one::<String>("threads")
            .and_then(|s| s.parse().ok());
        let headed_output = headed_in_original_dir(vcf_path, &align_out);
        println!(
            "[info] Auto-running 'header' on aligned output: {}",
            headed_output
        );
        let _hdr_out = header::header_run(
            &align_out,
            ref_path,
            threads_opt,
            Some(&headed_output),
            ignore_level,
        )?;
        // Remove the intermediate replaced/sorted file after headering
        if let Err(e) = fs::remove_file(&align_out) {
            eprintln!(
                "[warn] Failed to remove intermediate file {}: {}",
                align_out, e
            );
        }
        headed_output
    } else {
        println!("[info] --no-header set: skipping automatic header synthesis");
        align_out
    };

    println!("[info] All operations complete. Output written to {output_path_log}.");

    Ok(())
}
