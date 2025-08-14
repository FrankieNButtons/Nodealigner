use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

/// Stream the file twice:
/// 1) detect the earliest sample column where, for the first `same` variant lines,
///    col[i] == col[i-1] holds (i.e., duplicated content start).
/// 2) re-emit the VCF trimming to that column (exclusive).
///
/// Behavior:
/// - If cut_idx == 0 (no samples at all) or cut_idx <= 10 → keep only fixed cols 1..9
/// - If no duplication detected → keep all columns
/// - Header/meta lines (#...) are always preserved.
/// - Works for .vcf and .vcf.gz by sniffing extension.
pub fn run_cleanning(vcf_file: &str, same: &usize) -> Result<()> {
    if *same == 0 {
        // Nothing to detect; just stream through unchanged
        let mut inp = open_maybe_gz(vcf_file)?;
        let mut out = std::io::stdout().lock();
        std::io::copy(&mut inp, &mut out).context("streaming VCF")?;
        return Ok(());
    }

    // First pass: find cut index
    let cut_idx = detect_cut_idx(vcf_file, *same)
        .with_context(|| format!("detecting duplicated content start in {vcf_file}"))?;

    // Second pass: emit trimmed VCF to stdout
    let reader = BufReader::new(open_maybe_gz(vcf_file)?);
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());

    let mut chrom_seen = false;
    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                chrom_seen = true;
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 9 {
                    bail!("#CHROM header has fewer than 9 columns");
                }
                if let Some(cut) = cut_idx {
                    // 0-based: keep [0..9) fixed cols and [9..cut) samples
                    if cut <= 10 {
                        writeln!(out, "{}", fields[..9].join("\t"))?;
                    } else {
                        let mut kept = Vec::with_capacity(9 + (cut - 9));
                        kept.extend_from_slice(&fields[..9]);
                        kept.extend_from_slice(&fields[9..cut]);
                        writeln!(out, "{}", kept.join("\t"))?;
                    }
                } else {
                    writeln!(out, "{line}")?;
                }
            } else {
                writeln!(out, "{line}")?;
            }
            continue;
        }

        // Body line
        if !chrom_seen {
            bail!("Encountered variant line before #CHROM header");
        }
        if let Some(cut) = cut_idx {
            if cut <= 10 {
                // Only fixed columns
                let mut it = line.split('\t');
                // Collect first 9 fields
                let mut first9 = Vec::with_capacity(9);
                for _ in 0..9 {
                    match it.next() {
                        Some(v) => first9.push(v),
                        None => bail!("Variant line has fewer than 9 fields"),
                    }
                }
                writeln!(out, "{}", first9.join("\t"))?;
            } else {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 9 {
                    bail!("Variant line has fewer than 9 fields");
                }
                let mut kept = Vec::with_capacity(cut);
                kept.extend_from_slice(&fields[..9]);
                // Be robust if line has fewer sample cols than header suggested
                let end = cut.min(fields.len());
                if end > 9 {
                    kept.extend_from_slice(&fields[9..end]);
                }
                writeln!(out, "{}", kept.join("\t"))?;
            }
        } else {
            writeln!(out, "{line}")?;
        }
    }
    out.flush()?;
    Ok(())
}

fn open_maybe_gz(path: &str) -> Result<Box<dyn Read>> {
    if path.ends_with(".gz") {
        let f = File::open(path).with_context(|| format!("open {path}"))?;
        Ok(Box::new(MultiGzDecoder::new(f)))
    } else {
        let f = File::open(path).with_context(|| format!("open {path}"))?;
        Ok(Box::new(f))
    }
}

/// Return Some(cut_idx) where cut_idx is the 0-based column index at which duplicated
/// content (vs previous sample) begins, judged by the first `same` variant lines.
/// The index refers to absolute field index (including fixed 0..8), so:
/// - samples start at index 9
/// - we cut to keep [0..cut_idx)
/// If no duplication is detected, return None.
fn detect_cut_idx(path: &str, same: usize) -> Result<Option<usize>> {
    let reader = BufReader::new(open_maybe_gz(path)?);

    let mut chrom_fields: Option<Vec<String>> = None;
    let mut sample_count: usize = 0;

    // For columns >= 10 (0-based >= 9+1), keep counters of matches to their immediate left
    // We’ll allocate after seeing #CHROM.
    let mut counts: Vec<usize> = Vec::new();
    let mut active = false;
    let mut lines_seen: usize = 0;

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                let fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
                if fields.len() < 10 {
                    // No sample columns at all → nothing to cut
                    return Ok(Some(9)); // implies keep only 1..9 if caller chooses
                }
                sample_count = fields.len().saturating_sub(9);
                // counts for columns 10..end (indices 10..)
                // We count for each i >= 10 equality with i-1.
                if sample_count >= 2 {
                    counts = vec![0; sample_count - 1];
                }
                chrom_fields = Some(fields);
                active = true;
            }
            continue;
        }

        if !active {
            // Skipping any pre-header junk just in case
            continue;
        }

        // We only need to inspect the first `same` variant lines.
        if lines_seen >= same {
            break;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            // malformed line; treat as non-match for safety
            lines_seen += 1;
            continue;
        }
        // Compare samples i (absolute idx 9 + k) with previous sample (9 + k - 1)
        // counts[k-1] corresponds to absolute index (9 + k)
        if sample_count >= 2 {
            let total_samples_in_line = fields.len().saturating_sub(9);
            let max_k = std::cmp::min(sample_count, total_samples_in_line);
            // k is 1..max_k-1 → compare k with k-1
            for k in 1..max_k {
                let a = fields[9 + k];
                let b = fields[9 + k - 1];
                if a == b {
                    counts[k - 1] += 1;
                }
            }
        }
        lines_seen += 1;
    }

    if counts.is_empty() {
        // No or single sample column only → nothing to cut
        return Ok(None);
    }

    // Earliest absolute column index where counts[idx-10] == same
    // counts[0] corresponds to absolute index 10 (2nd sample)
    for (j, &c) in counts.iter().enumerate() {
        if c >= same {
            let absolute = 10 + j; // because j==0 → col 10
            return Ok(Some(absolute));
        }
    }
    Ok(None)
}
