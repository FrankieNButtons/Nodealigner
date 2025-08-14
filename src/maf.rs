use log::info;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Filter VCF by non-0/0 genotype proportion ("MAF-like") and print to stdout.
///
/// Keeps records where (#samples with GT != 0/0 and != 0|0) / (#samples with non-missing GT) > thresh.
/// Missing GTs (./.) are excluded from the denominator.
pub fn maf_main(matches: &clap::ArgMatches, output_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let vcf = matches.get_one::<String>("vcf").unwrap().as_str();
    let thresh = matches
        .get_one::<String>("thresh")
        .map(|s| s.parse::<f64>().unwrap_or(0.05))
        .unwrap_or(0.05);
    let _threads = matches.get_one::<String>("threads").map(|s| s.as_str()).unwrap_or("1");

    info!("Running 'gfa2bin-aligner maf'");
    // [INFO] style logging for arguments, aligned to match align_main example
    println!("[INFO]     --vcf     = {}", vcf);
    println!("[INFO]     --thresh  = {}", thresh);
    println!("[INFO]     --threads = {}", _threads);
    println!("[INFO]     --output  = {}", output_path);

    if vcf.ends_with(".gz") {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Compressed VCF (.gz) is not supported in this subcommand yet. Please decompress first (e.g., `bgzip -d file.vcf.gz`).",
        )));
    }

    let file = File::open(vcf)?;
    let reader = BufReader::new(file);

    let mut header_lines = Vec::new();
    let mut kept_chroms = HashSet::new();
    let mut variant_lines = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('#') {
            // Collect header lines
            header_lines.push(line);
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 {
            // No sample columns; nothing to filter on -> skip
            continue;
        }

        // FORMAT column and GT index
        let format = cols[8];
        let mut gt_idx: Option<usize> = None;
        for (i, key) in format.split(':').enumerate() {
            if key == "GT" {
                gt_idx = Some(i);
                break;
            }
        }
        let gt_idx = match gt_idx {
            Some(i) => i,
            _none => {
                // Without GT we can't compute; skip
                continue;
            }
        };

        let mut denom = 0usize;   // non-missing
        let mut count_11 = 0usize; // GT == 1/1 or 1|1
        let mut count_00 = 0usize; // GT == 0/0 or 0|0
        for sample in cols.iter().skip(9) {
            if sample.is_empty() {
                continue;
            }
            let parts: Vec<&str> = sample.split(':').collect();
            if gt_idx >= parts.len() {
                continue; // malformed
            }
            let gt = parts[gt_idx];
            // Missing?
            if gt == "./." || gt == ".|." || gt == "." {
                continue;
            }
            denom += 1;
            if gt == "1/1" || gt == "1|1" {
                count_11 += 1;
            } else if gt == "0/0" || gt == "0|0" {
                count_00 += 1;
            }
        }

        if denom == 0 {
            continue;
        }

        // If all samples are 1/1 or 0/0 (homozygous for either allele), remove the line
        if count_11 + count_00 == denom {
            continue;
        }

        // Keep the line if the proportion of 1/1 among non-missing samples is greater than thresh
        let prop_11 = (count_11 as f64) / (denom as f64);
        if prop_11 > thresh {
            kept_chroms.insert(cols[0].to_string());
            variant_lines.push(line);
        }
    }

    let file_out = File::create(output_path)?;
    let mut out = BufWriter::new(file_out);

    for header in header_lines {
        if header.starts_with("##contig=") {
            // Extract ID from contig line
            if let Some(id_start) = header.find("ID=") {
                let id_sub = &header[id_start + 3..];
                let id_end = id_sub.find(|c: char| c == ',' || c == '>').unwrap_or(id_sub.len());
                let chrom_id = &id_sub[..id_end];
                if !kept_chroms.contains(chrom_id) {
                    continue;
                }
            }
        }
        writeln!(out, "{}", header)?;
    }

    for variant_line in variant_lines {
        writeln!(out, "{}", variant_line)?;
    }

    out.flush()?;
    Ok(())
}
