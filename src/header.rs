
pub fn header_run(
    vcf_in: &str,
    reference_tsv: &str,
    threads: Option<usize>,
    output: Option<&str>,
    ignore: u8,
) -> Result<String, Box<dyn std::error::Error>> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
        println!("[info] Rayon thread pool set to {n} threads");
    }

    // Output path default: ./<basename-without-.vcf>.withheader.vcf (also handle .vcf.gz)
    let out_path = if let Some(o) = output { o.to_string() } else {
        let fname = std::path::Path::new(vcf_in)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
        let stem = if fname.ends_with(".vcf.gz") {
            fname.trim_end_matches(".vcf.gz")
        } else if fname.ends_with(".vcf") {
            fname.trim_end_matches(".vcf")
        } else {
            fname
        };
        format!("./{}.withheader.vcf", stem)
    };

    println!("[info] [header] --vcf {vcf_in}");
    println!("[info] [header] --reference {reference_tsv}");
    println!("[info] [header] --output {out_path}");

    // Reader supports plain text and .gz
    let infile = File::open(vcf_in)?;
    let ext = Path::new(vcf_in).extension().and_then(|s| s.to_str()).unwrap_or("");
    let reader: Box<dyn BufRead> = if ext.eq_ignore_ascii_case("gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(infile)))
    } else {
        Box::new(BufReader::new(infile))
    };

    // Temp spool for body
    let tmp_path = format!("{}.spool.tmp", out_path);
    let tmpf = File::create(&tmp_path)?;
    let tmpw = BufWriter::new(tmpf);

    // Read header + blocks; spool body as we go
    let block_cap = 100_000; // lines per block
    let (pre_header, column_header, blocks) = read_blocks_and_spool(reader, tmpw, block_cap)?;

    // Parallel inference
    let (inferred_info, inferred_fmt, first_data) = infer_from_blocks_parallel(blocks);

    // Contigs
    let contigs = parse_reference_tsv(reference_tsv)?;
    if contigs.is_empty() {
        eprintln!("[warn] reference.tsv produced no contigs; emitting header without lengths.");
    }

    // Track existing INFO/FORMAT/FILTER
    let mut existing_info: BTreeSet<String> = BTreeSet::new();
    let mut existing_format: BTreeSet<String> = BTreeSet::new();
    for l in &pre_header {
        if let Some(id) = l.strip_prefix("##INFO=<ID=").and_then(|x| x.split(',').next()) {
            existing_info.insert(id.to_string());
        } else if let Some(id) = l.strip_prefix("##FORMAT=<ID=").and_then(|x| x.split(',').next()) {
            existing_format.insert(id.to_string());
        }
    }

    // Final header
    let mut new_header: Vec<String> = Vec::new();
    let fileformat_line = pre_header
        .iter()
        .find(|l| l.starts_with("##fileformat="))
        .cloned()
        .unwrap_or_else(|| "##fileformat=VCFv4.2".to_string());
    println!("[info] Emitting fileformat line: {fileformat_line}");
    new_header.push(fileformat_line);
    new_header.push("##source=gfa2bin-aligner/header".to_string());

    // Apply ignore rules and aggregate by normalized contig ID (max length per ID)
    let mut norm_map: BTreeMap<String, u64> = BTreeMap::new();
    for (k, len) in &contigs {
        if let Some(id) = apply_ignore_rules(k, ignore) {
            let entry = norm_map.entry(id).or_insert(0);
            if *len > *entry { *entry = *len; }
        }
    }
    for (id, len) in norm_map {
        if len > 0 {
            new_header.push(format!("##contig=<ID={},length={}>", id, len));
        } else {
            new_header.push(format!("##contig=<ID={}>", id));
        }
    }

    for l in &pre_header {
        if l.starts_with("##INFO=<ID=") || l.starts_with("##FORMAT=<ID=") || l.starts_with("##FILTER=<") {
            new_header.push(l.clone());
        }
    }

    for (k, ks) in inferred_info {
        // Always add INFO fields present in the data, even GT (which will be ignored by downstream tools if present in INFO).
        if !existing_info.contains(&k) {
            // If no type/description can be inferred, allow empty description.
            // Use a flag to indicate that empty description is allowed.
            new_header.push(infer_info_def_with_empty(&k, &ks, true));
        }
    }
    for (k, (kind, card)) in inferred_fmt {
        if !existing_format.contains(&k) {
            new_header.push(infer_format_def(&k, Some(kind), Some(card)));
        }
    }

    if let Some(ch) = column_header {
        new_header.push(ch);
    } else {
        new_header.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT".to_string());
        if first_data.is_none() {
            eprintln!("[warn] No column header and no data line found early; downstream tools may reject this file.");
        }
    }

    // Write header, then append spooled body
    {
        let mut out = BufWriter::new(File::create(&out_path)?);
        for l in new_header { writeln!(out, "{}", l)?; }
        out.flush()?;
    }
    {
        use std::io::copy;
        let mut out_app = std::fs::OpenOptions::new().append(true).open(&out_path)?;
        let mut tmp_r = BufReader::new(File::open(&tmp_path)?);
        copy(&mut tmp_r, &mut out_app)?;
    }

    let _ = std::fs::remove_file(&tmp_path);

    println!("[info] Header synthesis complete → {out_path}");
    println!("[note] Streaming + parallel inference. Record-body normalization is not performed.");
    Ok(out_path)
}
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use clap::ArgMatches;
use std::path::Path;
use flate2::read::MultiGzDecoder;

type Contigs = BTreeMap<String, u64>;

fn parse_reference_tsv(p: &str) -> io::Result<Contigs> {
    let f = BufReader::new(File::open(p)?);
    let mut contigs: Contigs = BTreeMap::new();
    for (i, line) in f.lines().enumerate() {
        let l = line?;
        if l.trim().is_empty() { continue; }
        if i == 0 && l.starts_with("node\tstart\tend\tpath") { continue; }
        let mut it = l.split('\t');
        let _node = it.next().unwrap_or_default();
        let _start = it.next().unwrap_or_default();
        let end = it.next().unwrap_or("0").parse::<u64>().unwrap_or(0);
        let path = it.next().unwrap_or("").to_string();
        if path.is_empty() { continue; }
        contigs.entry(path).and_modify(|m| *m = (*m).max(end)).or_insert(end);
    }
    Ok(contigs)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ValKind { Int, Float, Stringy }

#[derive(Clone, Debug)]
struct KeyStats {
    // INFO
    seen_as_flag: bool,
    all_int: bool,
    any_float: bool,
    // Number inference
    all_singleton: bool,    // always one value
    matches_a: usize,       // how many lines length == #ALT
    matches_r: usize,       // how many lines length == #ALT+1
    samples: usize,         // how many lines we saw this key in
}

impl Default for KeyStats {
    fn default() -> Self {
        Self {
            seen_as_flag: false,
            all_int: true,
            any_float: false,
            all_singleton: true,
            matches_a: 0,
            matches_r: 0,
            samples: 0,
        }
    }
}

fn classify_value_token(tok: &str) -> ValKind {
    if tok.is_empty() { return ValKind::Stringy; }
    // Try integer
    if let Ok(_) = tok.parse::<i64>() { return ValKind::Int; }
    // Try float
    if let Ok(_) = tok.parse::<f64>() { return ValKind::Float; }
    ValKind::Stringy
}

fn infer_info_def_with_empty(id: &str, ks: &KeyStats, allow_empty: bool) -> String {
    // Type
    let typ = if ks.seen_as_flag {
        "Flag".to_string()
    } else if ks.all_int && !ks.any_float {
        "Integer".to_string()
    } else if ks.any_float {
        "Float".to_string()
    } else {
        "String".to_string()
    };

    // Number
    let number = if ks.seen_as_flag {
        "0".to_string()
    } else if ks.matches_a * 2 >= ks.samples && ks.samples > 0 {
        "A".to_string()
    } else if ks.matches_r * 2 >= ks.samples && ks.samples > 0 {
        "R".to_string()
    } else if ks.all_singleton {
        "1".to_string()
    } else {
        ".".to_string()
    };

    // Description
    let desc = if allow_empty && ks.samples == 0 {
        ""
    } else if allow_empty && typ == "" {
        ""
    } else {
        "Inferred from body"
    };
    // If allow_empty and desc is empty, emit empty Description
    if allow_empty && desc.is_empty() {
        format!("##INFO=<ID={id},Number={number},Type={typ},Description=\"\">")
    } else {
        format!("##INFO=<ID={id},Number={number},Type={typ},Description=\"{desc}\">")
    }
}

fn infer_format_def(id: &str, example_kind: Option<ValKind>, example_card: Option<usize>) -> String {
    let (number, _typ) = match (example_kind, example_card) {
        (Some(_), Some(0)) => ("0".to_string(), "String".to_string()),
        (Some(_), Some(1)) => ("1".to_string(), "String".to_string()),
        (Some(_), Some(_n)) => (".".to_string(), "String".to_string()),
        (Some(_), _none) => (".".to_string(), "String".to_string()),
        (_none, _) => (".".to_string(), "String".to_string()),
    };
    let (number, typ) = match example_kind {
        Some(ValKind::Int)   => (number, "Integer".to_string()),
        Some(ValKind::Float) => (number, "Float".to_string()),
        _ => (number, "String".to_string()),
    };
    format!("##FORMAT=<ID={id},Number={number},Type={typ},Description=\"Inferred from FORMAT column\">")
}

// --- Streaming + parallel inference helpers ---

#[inline]
fn merge_keystats(mut a: KeyStats, b: &KeyStats) -> KeyStats {
    a.seen_as_flag |= b.seen_as_flag;
    a.all_int &= b.all_int;
    a.any_float |= b.any_float;
    a.all_singleton &= b.all_singleton;
    a.matches_a += b.matches_a;
    a.matches_r += b.matches_r;
    a.samples += b.samples;
    a
}

#[inline]
fn merge_info_maps(mut a: BTreeMap<String, KeyStats>, b: BTreeMap<String, KeyStats>) -> BTreeMap<String, KeyStats> {
    for (k, v) in b {
        a.entry(k).and_modify(|x| *x = merge_keystats(std::mem::take(x), &v)).or_insert(v);
    }
    a
}

#[inline]
fn merge_format_maps(mut a: BTreeMap<String, (ValKind, usize)>, b: BTreeMap<String, (ValKind, usize)>) -> BTreeMap<String, (ValKind, usize)> {
    for (k, v) in b {
        a.entry(k).or_insert(v);
    }
    a
}

struct Block { lines: Vec<String> }

/// Read pre-header and column header, then read the body in blocks, while spooling all body lines to `spool_writer`.
fn read_blocks_and_spool<R: BufRead, W: Write>(mut reader: R, mut spool_writer: W, block_cap: usize)
    -> io::Result<(Vec<String>, Option<String>, Vec<Block>)>
{
    let mut pre_header: Vec<String> = Vec::new();
    let mut column_header: Option<String> = None;

    // Consume header lines
    loop {
        let mut buf = String::new();
        let n = reader.read_line(&mut buf)?;
        if n == 0 { break; }
        if buf.starts_with("##") {
            pre_header.push(buf.trim_end().to_string());
        } else if buf.starts_with("#CHROM\t") || buf.starts_with("#CHROM ") {
            column_header = Some(buf.trim_end().to_string());
            break;
        } else {
            // No proper header; treat buf as first data line
            let mut blk = Block { lines: Vec::with_capacity(block_cap) };
            blk.lines.push(buf);
            let mut blocks = Vec::new();
            loop {
                while blk.lines.len() < block_cap {
                    let mut l = String::new();
                    let n2 = reader.read_line(&mut l)?;
                    if n2 == 0 { break; }
                    blk.lines.push(l);
                }
                for l in &blk.lines { spool_writer.write_all(l.as_bytes())?; }
                blocks.push(blk);
                let mut peek = String::new();
                let n3 = reader.read_line(&mut peek)?;
                if n3 == 0 { break; }
                blk = Block { lines: Vec::with_capacity(block_cap) };
                blk.lines.push(peek);
            }
            return Ok((pre_header, column_header, blocks));
        }
    }

    // Normal path: we had header; now read the body in blocks
    let mut blocks: Vec<Block> = Vec::new();
    let mut blk = Block { lines: Vec::with_capacity(block_cap) };
    loop {
        let mut l = String::new();
        let n = reader.read_line(&mut l)?;
        if n == 0 {
            if !blk.lines.is_empty() {
                for x in &blk.lines { spool_writer.write_all(x.as_bytes())?; }
                blocks.push(blk);
            }
            break;
        }
        blk.lines.push(l);
        if blk.lines.len() == block_cap {
            for x in &blk.lines { spool_writer.write_all(x.as_bytes())?; }
            blocks.push(blk);
            blk = Block { lines: Vec::with_capacity(block_cap) };
        }
    }

    Ok((pre_header, column_header, blocks))
}

fn infer_from_blocks_parallel(blocks: Vec<Block>)
    -> (BTreeMap<String, KeyStats>, BTreeMap<String, (ValKind, usize)>, Option<String>)
{
    use rayon::prelude::*;

    let results = blocks.into_par_iter().map(|blk| {
        let mut info_map: BTreeMap<String, KeyStats> = BTreeMap::new();
        let mut fmt_map: BTreeMap<String, (ValKind, usize)> = BTreeMap::new();
        let mut first_data: Option<String> = None;

        for line in &blk.lines {
            let trimmed = line.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
            if first_data.is_none() { first_data = Some(trimmed.to_string()); }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 8 { continue; }
            let alt_ct = fields.get(4)
                .map(|alts| alts.split(',').filter(|x| !x.is_empty()).count())
                .unwrap_or(0);

            // INFO
            if let Some(info) = fields.get(7) {
                for item in info.split(';') {
                    if item.is_empty() { continue; }
                    if let Some((k, v)) = item.split_once('=') {
                        let ks = info_map.entry(k.to_string()).or_default();
                        ks.samples += 1;
                        if v.is_empty() { continue; }
                        let vals: Vec<&str> = v.split(',').collect();
                        ks.all_singleton &= vals.len() == 1;
                        if vals.len() == alt_ct { ks.matches_a += 1; }
                        if vals.len() == alt_ct + 1 { ks.matches_r += 1; }
                        for vv in &vals {
                            match classify_value_token(vv) {
                                ValKind::Int => { /* keep all_int */ }
                                ValKind::Float => { ks.any_float = true; ks.all_int = false; }
                                ValKind::Stringy => { ks.all_int = false; }
                            }
                        }
                    } else {
                        // Flag (no '=')
                        let ks = info_map.entry(item.to_string()).or_default();
                        ks.samples += 1;
                        ks.seen_as_flag = true;
                    }
                }
            }

            // FORMAT
            if fields.len() >= 10 {
                if let Some(fmt) = fields.get(8) {
                    for key in fmt.split(':') {
                        if key.is_empty() { continue; }
                        if let Some(sample) = fields.get(9) {
                            if let Some((pos, _)) = fmt.split(':').enumerate().find(|(_, k)| *k == key) {
                                let sample_tok = sample.split(':').nth(pos).unwrap_or("");
                                let card = sample_tok.split(',').filter(|x| !x.is_empty()).count();
                                let kind = if card == 0 { ValKind::Stringy } else {
                                    classify_value_token(sample_tok.split(',').next().unwrap_or(""))
                                };
                                fmt_map.entry(key.to_string()).or_insert((kind, card));
                            }
                        }
                    }
                }
            }
        }

        (info_map, fmt_map, first_data)
    }).collect::<Vec<_>>();

    let mut final_info: BTreeMap<String, KeyStats> = BTreeMap::new();
    let mut final_fmt: BTreeMap<String, (ValKind, usize)> = BTreeMap::new();
    let mut any_first: Option<String> = None;

    for (im, fm, fd) in results {
        final_info = merge_info_maps(final_info, im);
        final_fmt = merge_format_maps(final_fmt, fm);
        if any_first.is_none() { any_first = fd; }
    }

    (final_info, final_fmt, any_first)
}

pub fn header_main(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let vcf_in = matches.get_one::<String>("vcf").unwrap();
    let reference_tsv = matches.get_one::<String>("reference").unwrap();
    let threads = matches.get_one::<String>("threads").and_then(|s| s.parse::<usize>().ok());
    let out_opt = matches.get_one::<String>("output").map(|s| s.as_str());
    let ignore: u8 = matches.get_one::<String>("ignore")
        .and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
    let _ = header_run(vcf_in, reference_tsv, threads, out_opt, ignore)?;
    Ok(())
}

/// Helper function to apply ignore rules to contig names.
/// Apply ignore rules 0..=5:
/// 0: Keep as-is (no checks)
/// 1: Keep only if string (case-insensitively) contains "chr"
/// 2: As 1, and drop if token after "chr" is not [digits|X|Y|M]
/// 3: As 2, and drop if there is any suffix after the token (e.g., "chr12_random")
/// 4: Keep only standard human set {1..22, X, Y, M}. Normalize to "chr{TOKEN}" (uppercase token),
///    stripping any extra context, e.g. "GRCh38.chr12_random" -> "chr12".
/// 5: Same as 4, but output only "{TOKEN}" without "chr" prefix, e.g. "12", "X", "Y", "M".
fn apply_ignore_rules(raw: &str, level: u8) -> Option<String> {
    match level {
        0 => Some(raw.to_string()),
        1 => {
            let contains_chr = raw.to_ascii_lowercase().contains("chr");
            if contains_chr { Some(raw.to_string()) } else { None }
        }
        2 => {
            let (tok_opt, _has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr { return None; }
            match tok_opt {
                Some(t) => {
                    let ok = matches!(t.as_str(), "X" | "Y" | "M") || t.chars().all(|c| c.is_ascii_digit());
                    if ok { Some(raw.to_string()) } else { None }
                }
                none => none,
            }
        }
        3 => {
            let (tok_opt, has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr { return None; }
            if tok_opt.is_none() { return None; }
            if has_suffix { return None; }
            Some(raw.to_string())
        }
        4 => {
            let (tok_opt, _has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr { return None; }
            let t = tok_opt?;
            if !is_std_human_chr_token(&t) { return None; }
            Some(format!("chr{}", t.to_ascii_uppercase()))
        }
        5 => {
            let (tok_opt, _has_suffix, found_chr) = extract_chr_token(raw);
            if !found_chr { return None; }
            let t = tok_opt?;
            if !is_std_human_chr_token(&t) { return None; }
            Some(t.to_ascii_uppercase())
        }
        _ => Some(raw.to_string()),
    }
}

/// Extracts the token after "chr" (case-insensitive) in the string.
/// Returns (Some(token), has_suffix, found_chr)
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

/// Returns true if token is a standard human chromosome (1-22, X, Y, M).
fn is_std_human_chr_token(tok: &str) -> bool {
    match tok.to_ascii_uppercase().as_str() {
        "X" | "Y" | "M" => true,
        n => n.parse::<u8>().map(|v| (1..=22).contains(&v)).unwrap_or(false),
    }
}