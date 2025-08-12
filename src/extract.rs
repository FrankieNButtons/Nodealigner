use clap::ArgMatches;
use gfa_reader::Gfa;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Extracts the token that follows a case-insensitive "chr" occurrence.
/// Returns (token_opt, has_suffix_after_token, found_chr_anywhere)
/// - token is uppercased (e.g., "x" -> "X")
/// - token is either a run of digits (e.g., "12") or a single letter in {X,Y,M}
fn extract_chr_token(raw: &str) -> (Option<String>, bool, bool) {
    let lower = raw.to_ascii_lowercase();
    if let Some(idx) = lower.find("chr") {
        let after = &raw[idx + 3..];
        let mut chars = after.chars();
        if let Some(first) = chars.next() {
            // Accept a run of digits OR a single X/Y/M letter
            if first.is_ascii_digit() {
                let mut token = String::new();
                token.push(first);
                for c in chars.clone() {
                    if c.is_ascii_digit() {
                        token.push(c);
                    } else {
                        break;
                    }
                }
                let token_len = token.len();
                let has_suffix = after[token_len + 0..].chars().next().is_some();
                return (Some(token), has_suffix, true);
            } else {
                let up = first.to_ascii_uppercase();
                if up == 'X' || up == 'Y' || up == 'M' {
                    let token = up.to_string();
                    let has_suffix = after[1..].chars().next().is_some();
                    return (Some(token), has_suffix, true);
                }
            }
        }
        return (None, false, true);
    }
    (None, false, false)
}

/// Standard human chromosome token check: {1..22, X, Y, M}
fn is_std_human_chr_token(tok: &str) -> bool {
    let t = tok.to_ascii_uppercase();
    if t == "X" || t == "Y" || t == "M" {
        return true;
    }
    if t.chars().all(|c| c.is_ascii_digit()) {
        if let Ok(n) = t.parse::<u32>() {
            return (1..=22).contains(&n);
        }
    }
    false
}

/// Apply ignore rules 0..=5 (see doc for details)
fn apply_ignore_rules(raw: &str, level: u8) -> Option<String> {
    match level {
        0 => Some(raw.to_string()),
        1 => {
            let contains_chr = raw.to_ascii_lowercase().contains("chr");
            if contains_chr {
                Some(raw.to_string())
            } else {
                None
            }
        }
        2 => {
            let (tok_opt, _has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr {
                return None;
            }
            match tok_opt {
                Some(t) => {
                    let ok = matches!(t.as_str(), "X" | "Y" | "M")
                        || t.chars().all(|c| c.is_ascii_digit());
                    if ok { Some(raw.to_string()) } else { None }
                }
                none => none,
            }
        }
        3 => {
            let (tok_opt, has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr {
                return None;
            }
            if tok_opt.is_none() {
                return None;
            }
            if has_suffix {
                return None;
            }
            Some(raw.to_string())
        }
        4 => {
            let (tok_opt, _has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr {
                return None;
            }
            let t = tok_opt?;
            if !is_std_human_chr_token(&t) {
                return None;
            }
            Some(format!("chr{}", t.to_ascii_uppercase()))
        }
        5 => {
            let (tok_opt, _has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr {
                return None;
            }
            let t = tok_opt?;
            if !is_std_human_chr_token(&t) {
                return None;
            }
            Some(t.to_ascii_uppercase())
        }
        _ => Some(raw.to_string()),
    }
}

/// Extract paths and node coordinates from GFA
/// Writes a TSV with columns: node, start, end, seq, length, path (supports P & W by converting W to paths).
pub fn extract_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let gfa_file = matches
        .get_one::<String>("gfa")
        .expect("--gfa is required")
        .to_string();

    // Default output: same directory as the input GFA, file name "reference.tsv"
    let output_file: String = matches
        .get_one::<String>("output")
        .cloned()
        .unwrap_or_else(|| {
            let p = Path::new(&gfa_file);
            let dir = p.parent().unwrap_or_else(|| Path::new("."));
            dir.join("reference.tsv").to_string_lossy().into_owned()
        });

    let output_file_clone = output_file.clone();

    // Optional threads
    let threads = matches
        .get_one::<String>("threads")
        .cloned()
        .unwrap_or_default();
    let num_threads: usize = threads.parse().unwrap_or(1);

    // Optional ignore level (we'll register the CLI in main.rs later). Accept either "ignore" or "ignore-level".
    let ignore_level: u8 = matches
        .get_one::<String>("ignore")
        .or_else(|| matches.get_one::<String>("ignore-level"))
        .and_then(|s| s.parse::<u8>().ok())
        .unwrap_or(0);

    println!("[info] [extract] Running with arguments:");
    println!("    --gfa     : {}", gfa_file);
    println!("    --output  : {}", output_file);
    println!(
        "    --threads : {}",
        if threads.is_empty() {
            "1 (default)"
        } else {
            threads.as_str()
        }
    );
    println!("    --ignore  : {}", ignore_level);

    if num_threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .ok();
        println!("[info] Rayon thread pool set to {num_threads} threads");
    }

    println!("[info] Reading GFA and extracting path information...");
    let mut graph: Gfa<u32, (), ()> = Gfa::parse_gfa_file_multi(&gfa_file, num_threads);

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
        // Apply ignore rules to path name; skip entire path if it does not pass.
        let maybe_name = apply_ignore_rules(&path.name, ignore_level);
        let out_name: String = match maybe_name {
            Some(s) => s,
            _none => return,
        };
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
                node, start, end, seq, len, out_name
            );
            start = end;
        }
        if !local_buf.is_empty() {
            let mut guard = out.lock().unwrap();
            let _ = guard.write_all(local_buf.as_bytes());
        }
    });
    println!("[info] Extraction complete. Output written to {output_file_clone}.");
    Ok(())
}
