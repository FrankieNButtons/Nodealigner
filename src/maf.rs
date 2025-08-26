use log::info;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use rayon::prelude::*;

/// Filter VCF by per-GT proportions and print to stdout.
///
/// For each variant line, consider only non-missing diploid genotypes whose alleles are in {0,1}.
/// Normalize phased-symbol to unphased-symbol ("|" -> "/") **but keep allele order** so that 0/1 and 1/0 are
/// treated as distinct categories. Let denom be the count of such valid calls. For each GT category in
/// {"0/0","0/1","1/0","1/1"} that actually appears (count > 0), compute its proportion p = count / denom.
/// Keep the line iff every appearing category's p lies within the closed interval [thresh, 1 - thresh].
/// Missing (./.) and non-{0,1} allele genotypes are excluded from denom.
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
    println!("[INFO]     rule     = for GT in {{0/0,0/1,1/0,1/1}} that appear: each proportion in [thresh, 1-thresh]");

    if vcf.ends_with(".gz") {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Compressed VCF (.gz) is not supported in this subcommand yet. Please decompress first (e.g., `bgzip -d file.vcf.gz`).",
        )));
    }

    let file = File::open(vcf)?;
    let reader = BufReader::new(file);
    let all_lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

    let mut header_lines: Vec<String> = Vec::new();
    let mut variants: Vec<(usize, String)> = Vec::new();
    for (idx, line) in all_lines.into_iter().enumerate() {
        if line.starts_with('#') {
            header_lines.push(line);
        } else {
            variants.push((idx, line));
        }
    }

    // Process variant lines in parallel, but keep (idx) to restore original order
    let kept: Vec<(usize, String, String)> = variants
        .par_iter()
        .filter_map(|(idx, line)| {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 10 {
                return None; // no samples
            }

            // FORMAT column and GT index
            let format = cols[8];
            let mut gt_idx_opt: Option<usize> = None;
            for (i, key) in format.split(':').enumerate() {
                if key == "GT" {
                    gt_idx_opt = Some(i);
                    break;
                }
            }
            let gt_idx = if let Some(i) = gt_idx_opt {
                i
            } else {
                return None;
            };

            let mut denom = 0usize; // # valid non-missing GT with alleles in {0,1}
            let mut c_00 = 0usize;
            let mut c_01 = 0usize; // strictly 0/1
            let mut c_10 = 0usize; // strictly 1/0
            let mut c_11 = 0usize;

            for sample in cols.iter().skip(9) {
                if sample.is_empty() { continue; }
                let parts: Vec<&str> = sample.split(':').collect();
                if gt_idx >= parts.len() { continue; }
                let gt = parts[gt_idx];

                // Skip missing GT
                if gt == "./." || gt == ".|." || gt == "." { continue; }

                // Normalize phased to unphased (keep allele order) without borrowing a temporary
                let gt_norm: String = if gt.contains('|') { gt.replace('|', "/") } else { gt.to_string() };

                let ab: Vec<&str> = gt_norm.split('/').collect();
                if ab.len() != 2 { continue; }
                let a = ab[0];
                let b = ab[1];

                // Only keep biallelic {0,1} diploid genotypes
                if (a != "0" && a != "1") || (b != "0" && b != "1") { continue; }

                denom += 1;
                if a == "0" && b == "0" {
                    c_00 += 1;
                } else if a == "0" && b == "1" {
                    c_01 += 1;
                } else if a == "1" && b == "0" {
                    c_10 += 1;
                } else if a == "1" && b == "1" {
                    c_11 += 1;
                }
            }

            if denom == 0 { return None; }

            let check = |cnt: usize, denom: usize, thresh: f64| -> bool {
                let p = (cnt as f64) / (denom as f64);
                p >= thresh && p <= (1.0 - thresh)
            };

            // For each category that appears, require its proportion to be within [thresh, 1-thresh]
            if (c_00 > 0 && !check(c_00, denom, thresh))
                || (c_01 > 0 && !check(c_01, denom, thresh))
                || (c_10 > 0 && !check(c_10, denom, thresh))
                || (c_11 > 0 && !check(c_11, denom, thresh))
            {
                return None;
            }

            // Keep: return (idx, chrom, line)
            Some((*idx, cols[0].to_string(), line.clone()))
        })
        .collect::<Vec<(usize, String, String)>>();

    // Restore original VCF order and gather kept chroms
    let mut kept = kept; // make mutable
    kept.sort_by_key(|(idx, _, _)| *idx);

    let mut kept_chroms: HashSet<String> = HashSet::new();
    for (_, chrom, _) in kept.iter() {
        kept_chroms.insert(chrom.clone());
    }

    let variant_lines: Vec<String> = kept.into_iter().map(|(_, _, line)| line).collect();

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
