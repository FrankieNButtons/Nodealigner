use clap::ArgMatches;
use std::cmp::Ordering;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

fn keyword_to_index(prefix: &str) -> Option<usize> {
    match prefix {
        "CHROM" | "chrom" | "#CHROM" | "#chrom" => Some(0),
        "POS" | "pos" => Some(1),
        "ID" | "id" => Some(2),
        "REF" | "ref" => Some(3),
        "ALT" | "alt" => Some(4),
        "QUAL" | "qual" => Some(5),
        "FILTER" | "filter" => Some(6),
        "INFO" | "info" => Some(7),
        "FORMAT" | "format" => Some(8),
        _ => None,
    }
}

pub fn sort_main(matches: &ArgMatches) -> Result<(), Box<dyn Error>> {
    // clap v4 accessors
    let vcf_path = matches
        .get_one::<String>("vcf")
        .map(|s| s.as_str())
        .ok_or("--vcf is required")?;
    let prefix_key = matches
        .get_one::<String>("prefix")
        .map(|s| s.to_string())
        .unwrap_or_else(|| "POS".to_string());
    let reverse = matches.get_flag("reverse");

    // Determine column index from keyword or numeric string (default POS)
    let col_index: usize = keyword_to_index(&prefix_key)
        .or_else(|| prefix_key.parse::<usize>().ok())
        .unwrap_or(1);

    eprintln!(
        "[info] Sorting VCF: {} by '{}' (index {}){}",
        vcf_path,
        prefix_key,
        col_index,
        if reverse { " (reverse)" } else { "" }
    );

    let file = File::open(vcf_path)?;
    let reader = BufReader::new(file);

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

    data_lines.sort_by(|a, b| {
        let a_fields: Vec<&str> = a.split('\t').collect();
        let b_fields: Vec<&str> = b.split('\t').collect();
        let a_key = a_fields.get(col_index).copied().unwrap_or("");
        let b_key = b_fields.get(col_index).copied().unwrap_or("");
        let ord = match (a_key.parse::<f64>(), b_key.parse::<f64>()) {
            (Ok(na), Ok(nb)) => na.partial_cmp(&nb).unwrap_or(Ordering::Equal),
            _ => a_key.cmp(b_key),
        };
        if reverse { ord.reverse() } else { ord }
    });

    // Write to stdout
    let mut out = std::io::BufWriter::new(std::io::stdout());
    for h in header_lines { writeln!(out, "{}", h)?; }
    for d in &data_lines { writeln!(out, "{}", d)?; }

    eprintln!("[info] Sorting complete: {} records", data_lines.len());
    Ok(())
}