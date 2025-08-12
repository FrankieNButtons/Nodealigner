// src/sort_main.rs
use clap::ArgMatches;
use flate2::read::MultiGzDecoder;
use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Build <orig_stem>.sorted.vcf in the *original VCF's directory*.
/// If `candidate_name` is provided (from --output), use its base name but
/// still place it in the original directory and force the suffix to `.sorted.vcf`.
fn sorted_in_original_dir(original_vcf: &str, candidate_name: Option<&str>) -> String {
    use std::path::Path;
    let p = Path::new(original_vcf);
    let dir = p.parent().unwrap_or_else(|| Path::new("."));

    let base = if let Some(name) = candidate_name {
        let fname = Path::new(name)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("output.vcf");
        if fname.ends_with(".vcf.gz") {
            fname.trim_end_matches(".vcf.gz").to_string() + ".sorted.vcf"
        } else if fname.ends_with(".vcf") {
            fname.trim_end_matches(".vcf").to_string() + ".sorted.vcf"
        } else {
            fname.to_string() + ".sorted.vcf"
        }
    } else {
        // Derive from the original file's base name
        let fname = p
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("output.vcf");
        if fname.ends_with(".vcf.gz") {
            fname.trim_end_matches(".vcf.gz").to_string() + ".sorted.vcf"
        } else if fname.ends_with(".vcf") {
            fname.trim_end_matches(".vcf").to_string() + ".sorted.vcf"
        } else {
            let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("output");
            format!("{stem}.sorted.vcf")
        }
    };

    dir.join(base).to_string_lossy().into_owned()
}

/// Extract a token immediately following any "chr" substring (case-insensitive).
/// Returns (Some(token), has_suffix, found_chr)
#[inline]
fn extract_chr_token(raw: &str) -> (Option<String>, bool, bool) {
    let raw_lower = raw.to_ascii_lowercase();
    if let Some(idx) = raw_lower.find("chr") {
        let after = &raw[idx + 3..];
        let mut token = String::new();
        let mut has_suffix = false;
        for c in after.chars() {
            if c.is_ascii_alphanumeric() {
                token.push(c);
            } else {
                has_suffix = true;
                break;
            }
        }
        if token.is_empty() {
            return (None, has_suffix, true);
        }
        let tok_norm = if token.chars().all(|c| c.is_ascii_digit()) {
            token
        } else {
            token.to_ascii_uppercase()
        };
        (Some(tok_norm), has_suffix, true)
    } else {
        (None, false, false)
    }
}

/// Map any CHROM string to a rank and a normalized display like chr1..chr22, chrX, chrY, chrM
/// Rank: 1..=22 -> 1..=22, X->23, Y->24, M/MT->25.
/// Returns None if the chromosome cannot be interpreted (those will sort after known ones).
#[inline]
fn chrom_rank_and_display(raw: &str) -> Option<(u8, String)> {
    // Prefer token after "chr"
    if let (Some(tok), _suf, _found_chr) = extract_chr_token(raw) {
        let t = tok.to_ascii_uppercase();
        if let Ok(n) = t.parse::<u8>() {
            if (1..=22).contains(&n) {
                return Some((n, format!("chr{}", n)));
            }
        }
        return match t.as_str() {
            "X" => Some((23, "chrX".to_string())),
            "Y" => Some((24, "chrY".to_string())),
            "M" | "MT" => Some((25, "chrM".to_string())),
            _ => None,
        };
    }

    // Accept bare tokens: "1", "2", ..., "22", "X", "Y", "M", "MT"
    let t = raw.trim();
    if t.chars().all(|c| c.is_ascii_digit()) {
        if let Ok(n) = t.parse::<u8>() {
            if (1..=22).contains(&n) {
                return Some((n, format!("chr{}", n)));
            }
        }
    }
    match t.to_ascii_uppercase().as_str() {
        "X" => Some((23, "chrX".to_string())),
        "Y" => Some((24, "chrY".to_string())),
        "M" | "MT" => Some((25, "chrM".to_string())),
        _ => None,
    }
}

/// Universal comparator for VCF lines (body only).
/// Sorts by CHROM rank (1..22, X, Y, M), then by POS (numeric), then by ID (string).
#[inline]
fn cmp_vcf_records(a: &str, b: &str) -> Ordering {
    let fa: Vec<&str> = a.split('\t').collect();
    let fb: Vec<&str> = b.split('\t').collect();

    let ca = *fa.get(0).unwrap_or(&"");
    let cb = *fb.get(0).unwrap_or(&"");

    let ra = chrom_rank_and_display(ca);
    let rb = chrom_rank_and_display(cb);

    match (ra, rb) {
        (Some((r1, _)), Some((r2, _))) => {
            if r1 != r2 {
                return r1.cmp(&r2);
            }
            let pa = fa.get(1).and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            let pb = fb.get(1).and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            if pa != pb {
                return pa.cmp(&pb);
            }
            let ia = *fa.get(2).unwrap_or(&"");
            let ib = *fb.get(2).unwrap_or(&"");
            ia.cmp(ib)
        }
        (Some(_), _) => std::cmp::Ordering::Less,
        (_, Some(_)) => std::cmp::Ordering::Greater,
        (_, _) => {
            // Fallback: lexicographic CHROM, then numeric POS
            let ord = ca.cmp(cb);
            if ord != std::cmp::Ordering::Equal {
                return ord;
            }
            let pa = fa.get(1).and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            let pb = fb.get(1).and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            pa.cmp(&pb)
        }
    }
}

pub fn sort_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // Adjust these argument names if your CLI uses different flags
    let input = matches
        .get_one::<String>("vcf")
        .or_else(|| matches.get_one::<String>("input"))
        .ok_or("Missing --vcf/--input")?;
    let candidate_name = matches.get_one::<String>("output").map(|s| s.as_str());
    let output = sorted_in_original_dir(input, candidate_name);

    println!("[info] [sort] --vcf {input}");
    println!("[info] [sort] --output {output}");

    // Open input (supports .vcf and .vcf.gz)
    let infile = File::open(input)?;
    let ext = Path::new(input)
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    let reader: Box<dyn BufRead> = if ext.eq_ignore_ascii_case("gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(infile)))
    } else {
        Box::new(BufReader::new(infile))
    };

    // Collect header and body
    let mut pre_header: Vec<String> = Vec::new(); // lines starting with "##"
    let mut col_header: Option<String> = None; // line starting with "#CHROM"
    let mut body: Vec<String> = Vec::new(); // variant lines

    for line in reader.lines() {
        let l = line?;
        if l.starts_with("##") {
            pre_header.push(l);
        } else if l.starts_with("#CHROM\t") || l.starts_with("#CHROM ") {
            col_header = Some(l);
        } else {
            // everything after #CHROM is body; if a malformed file had data before, we still treat as body
            body.push(l);
        }
    }

    // Sort body with the universal comparator
    // For very large files, consider external merge sort; this keeps it simple and fast for typical sizes.
    body.sort_by(|a, b| cmp_vcf_records(a, b));

    // Write output
    let mut w = BufWriter::new(File::create(&output)?);
    for h in &pre_header {
        writeln!(w, "{}", h)?;
    }
    if let Some(h) = &col_header {
        writeln!(w, "{}", h)?;
    } else {
        // Fallback if there's no column header (rare, non-compliant VCF)
        writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    }
    for l in &body {
        writeln!(w, "{}", l)?;
    }
    w.flush()?;

    println!("[info] [sort] Done → {output}");
    Ok(())
}
