use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{File, rename};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

#[cfg(feature = "rayon")]
use rayon::ThreadPoolBuilder;
#[cfg(feature = "rayon")]
use rayon::prelude::*;

const MAP_CHUNK_LINES: usize = 100_000;

fn norm_chr(raw: &str) -> Option<String> {
    // strip common prefixes and normalize MT->M
    let mut c = raw.trim();
    if let Some(s) = c.strip_prefix("chr") {
        c = s;
    }
    if let Some(s) = c.strip_prefix("CHR") {
        c = s;
    }
    let c = match c {
        "Mt" | "MT" | "mt" => "M",
        other => other,
    };
    // Only allow 1..22, X, Y, M
    if c == "X" || c == "Y" || c == "M" {
        return Some(c.to_string());
    }
    if let Ok(n) = c.parse::<u32>() {
        if (1..=22).contains(&n) {
            return Some(n.to_string());
        }
    }
    None
}

fn make_key(chrom: &str, pos: &str, ref_allele: &str, alt_allele: &str) -> Option<String> {
    let chr = norm_chr(chrom)?;
    Some(format!("{}:{}:{}:{}", chr, pos, ref_allele, alt_allele))
}

fn is_probable_header_tsv(line: &str) -> bool {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 2 {
        return false;
    }
    let second = cols[1];
    second.eq_ignore_ascii_case("snp")
        || second.eq_ignore_ascii_case("variant")
        || second.eq_ignore_ascii_case("node")
}

fn try_normalize_token(token: &str) -> Option<String> {
    // Accept forms like (chr)10:12910:G:A
    let t = token.trim();
    if t.is_empty() {
        return None;
    }
    // Some FastQTL SNPs are like 10:12910:G:A, others can be rsIDs; normalize only when it's the former
    let parts: Vec<&str> = t.split(':').collect();
    if parts.len() >= 4 {
        let chrom = parts[0];
        let pos = parts[1];
        let r#ref = parts[2];
        let alt = parts[3];
        return make_key(chrom, pos, r#ref, alt);
    }
    None
}

fn load_map_from_tsv(path: &Path) -> std::io::Result<HashMap<String, String>> {
    let f = File::open(path)?;
    let r = BufReader::new(f);
    let mut map = HashMap::with_capacity(1_000_000);
    for line in r.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        if let (Some(k), Some(v)) = (it.next(), it.next()) {
            map.insert(k.to_string(), v.to_string());
        }
    }
    Ok(map)
}

pub fn run_rename(vcf_path: &str, qtl_path: &str, threads: usize) -> Result<(), Box<dyn Error>> {
    let t0 = Instant::now();
    eprintln!(
        "[INFO] rename: start vcf='{}' qtl='{}' threads={}",
        vcf_path, qtl_path, threads
    );

    // Check for --map-mode flag in args
    let args: Vec<String> = std::env::args().collect();
    let mut map_mode = "ids-only".to_string();
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--map-mode" && i + 1 < args.len() {
            map_mode = args[i + 1].clone();
            i += 1;
        }
        i += 1;
    }
    eprintln!("[INFO] Map mode: {}", map_mode);

    #[cfg(feature = "rayon")]
    {
        if threads > 0 {
            let _ = ThreadPoolBuilder::new().num_threads(threads).build_global();
        }
    }

    // Decide target directory (place map next to the QTL file)
    let qtl_dir = Path::new(qtl_path)
        .parent()
        .unwrap_or_else(|| Path::new("."));
    let map_final = qtl_dir.join("map.tsv");
    let map_tmp = qtl_dir.join("map.tsv.tmp");

    // If a stable map exists and no tmp in progress, reuse it. Otherwise build from VCF and write atomically.
    let id_key_map: HashMap<String, String> = if map_final.exists() && !map_tmp.exists() {
        eprintln!("[INFO] detected existing map: {}", map_final.display());
        let m = load_map_from_tsv(&map_final)?;
        eprintln!("[INFO] loaded existing map entries: {}", m.len());
        m
    } else {
        eprintln!("[INFO] building key→normalized map from VCF ...");
        let mut map: HashMap<String, String> = HashMap::with_capacity(1_000_000);
        let mut pos_map: HashMap<String, String> = HashMap::with_capacity(1_000_000);
        let mut pos_dupe: HashSet<String> = HashSet::with_capacity(64);

        // Create temp map file and write as we go
        let tmp_path = qtl_dir.join("map.tsv.tmp");
        eprintln!("[INFO] creating temp map at {}", tmp_path.display());
        let tmp_file = File::create(&tmp_path)?;
        let mut map_writer = BufWriter::new(tmp_file);
        let mut map_lines_written: usize = 0;
        let mut first_entries: usize = 0;

        let mut vcf_body_lines: usize = 0;
        let mut vcf_ids_added: usize = 0;
        let mut vcf_keys_added: usize = 0;
        let mut vcf_pos_added: usize = 0;
        let mut vcf_pos_ambiguous: usize = 0;

        {
            let vf = File::open(vcf_path)?;
            let mut vr = BufReader::new(vf);
            let mut line = String::new();
            let mut chunk: Vec<String> = Vec::with_capacity(MAP_CHUNK_LINES);
            let mut chunk_idx: usize = 0;
            let ids_only = map_mode == "ids-only";

            let mut process_chunk = |chunk_in: &Vec<String>| -> Result<(), Box<dyn Error>> {
                chunk_idx += 1;

                #[cfg(feature = "rayon")]
                let triples: Vec<(
                    Option<(String, String)>,
                    Option<(String, String)>,
                    Option<(String, String)>,
                )> = chunk_in
                    .par_iter()
                    .map(|l| {
                        let cols: Vec<&str> = l.split('\t').collect();
                        if cols.len() < 5 {
                            return (None, None, None);
                        }
                        let chrom = cols[0];
                        let pos = cols[1].trim();
                        let id = cols.get(2).copied().unwrap_or(".");
                        let r#ref = cols[3];
                        let alt_all = cols[4];
                        let alt = alt_all.split(',').next().unwrap_or(alt_all);
                        if let Some(norm) = make_key(chrom, pos, r#ref, alt) {
                            let id_pair = if id != "." && !id.is_empty() {
                                Some((id.to_string(), norm.clone()))
                            } else {
                                None
                            };
                            let key_pair = Some((norm.clone(), norm.clone()));
                            let pos_pair = Some((pos.to_string(), norm.clone()));
                            (id_pair, key_pair, pos_pair)
                        } else {
                            (None, None, None)
                        }
                    })
                    .collect();

                #[cfg(not(feature = "rayon"))]
                let triples: Vec<(
                    Option<(String, String)>,
                    Option<(String, String)>,
                    Option<(String, String)>,
                )> = chunk_in
                    .iter()
                    .map(|l| {
                        let cols: Vec<&str> = l.split('\t').collect();
                        if cols.len() < 5 {
                            return (None, None, None);
                        }
                        let chrom = cols[0];
                        let pos = cols[1].trim();
                        let id = cols.get(2).copied().unwrap_or(".");
                        let r#ref = cols[3];
                        let alt_all = cols[4];
                        let alt = alt_all.split(',').next().unwrap_or(alt_all);
                        if let Some(norm) = make_key(chrom, pos, r#ref, alt) {
                            let id_pair = if id != "." && !id.is_empty() {
                                Some((id.to_string(), norm.clone()))
                            } else {
                                None
                            };
                            let key_pair = Some((norm.clone(), norm.clone()));
                            let pos_pair = Some((pos.to_string(), norm.clone()));
                            (id_pair, key_pair, pos_pair)
                        } else {
                            (None, None, None)
                        }
                    })
                    .collect();

                let mut local_ids_added = 0usize;
                let mut local_keys_added = 0usize;
                let mut local_pos_added = 0usize;
                let mut local_pos_amb = 0usize;

                for (id_opt, key_opt, pos_opt) in triples {
                    // Only store id→normalized mapping by default
                    if ids_only {
                        // Only id→normalized
                        if let Some((idk, v)) = id_opt {
                            let id_is_numeric =
                                !idk.is_empty() && idk.chars().all(|c| c.is_ascii_digit());
                            if id_is_numeric {
                                if map.insert(idk.clone(), v.clone()).is_none() {
                                    writeln!(map_writer, "{}\t{}", idk, v)?;
                                    map_lines_written += 1;
                                    local_ids_added += 1;
                                    vcf_ids_added += 1;
                                    if first_entries < 5 {
                                        eprintln!("[INFO] Writing map entry: {} -> {}", idk, v);
                                        first_entries += 1;
                                    }
                                }
                            }
                        }
                    } else {
                        // Store all three mappings (id, normalized, pos)
                        if let Some((k, v)) = key_opt.clone() {
                            if map.insert(k, v).is_none() {
                                local_keys_added += 1;
                                vcf_keys_added += 1;
                            }
                        }
                        if let Some((idk, v)) = id_opt.clone() {
                            if map.insert(idk.clone(), v.clone()).is_none() {
                                writeln!(map_writer, "{}\t{}", idk, v)?;
                                map_lines_written += 1;
                                local_ids_added += 1;
                                vcf_ids_added += 1;
                                if first_entries < 5 {
                                    eprintln!("[INFO] Writing map entry: {} -> {}", idk, v);
                                    first_entries += 1;
                                }
                            }
                        }
                        if let Some((posk, normv)) = pos_opt.clone() {
                            if !pos_dupe.contains(&posk) {
                                if let Some(prev) = pos_map.get(&posk) {
                                    if prev != &normv {
                                        pos_map.remove(&posk);
                                        pos_dupe.insert(posk);
                                        local_pos_amb += 1;
                                        vcf_pos_ambiguous += 1;
                                    }
                                } else {
                                    pos_map.insert(posk, normv);
                                    local_pos_added += 1;
                                    vcf_pos_added += 1;
                                }
                            }
                        }
                    }
                }

                vcf_body_lines += chunk_in.len();
                eprintln!(
                    "[INFO] map-chunk #{}, lines={}, ids+={}, keys+={}, pos_unique+={}, pos_amb+={}",
                    chunk_idx,
                    chunk_in.len(),
                    local_ids_added,
                    local_keys_added,
                    local_pos_added,
                    local_pos_amb
                );
                Ok(())
            };

            loop {
                line.clear();
                let n = vr.read_line(&mut line)?;
                if n == 0 {
                    break;
                }
                if line.starts_with('#') {
                    continue;
                }
                let trimmed = line.trim_end_matches(['\n', '\r']).to_string();
                if trimmed.is_empty() {
                    continue;
                }
                chunk.push(trimmed);
                if chunk.len() >= MAP_CHUNK_LINES {
                    process_chunk(&chunk)?;
                    chunk = Vec::with_capacity(MAP_CHUNK_LINES);
                    eprintln!("[INFO] map-chunk dropped to release memory");
                }
            }
            if !chunk.is_empty() {
                process_chunk(&chunk)?;
                eprintln!("[INFO] final map-chunk dropped to release memory");
            }

            // Do not write POS-only entries to the map file
            map_writer.flush()?;
            eprintln!("[INFO] temp map lines written: {}", map_lines_written);
        }
        eprintln!(
            "[INFO] VCF pass: lines={} ids_mapped={} keys_mapped={} pos_unique={} pos_ambiguous={}",
            vcf_body_lines, vcf_ids_added, vcf_keys_added, vcf_pos_added, vcf_pos_ambiguous
        );
        if map_mode == "ids-only" {
            eprintln!(
                "[INFO] map.tsv policy: only numeric VCF IDs are written as keys; normalized and POS keys are not stored"
            );
        } else {
            eprintln!("[INFO] map.tsv policy: all mappings (id, normalized, pos) are stored");
        }

        // Merge pos_map into map for a unified lookup (only if not ids-only)
        if map_mode != "ids-only" {
            for (k, v) in pos_map.into_iter() {
                map.entry(k).or_insert(v);
            }
        }

        eprintln!("[INFO] Map build complete: {} entries", map.len());
        // Atomic rename temp -> final
        let final_map_path = qtl_dir.join("map.tsv");
        rename(&tmp_path, &final_map_path)?;
        eprintln!("[INFO] wrote map to {}", final_map_path.display());
        eprintln!("[INFO] map finalized; starting QTL phase");
        map
    };
    eprintln!("[INFO] active map size: {} entries", id_key_map.len());
    eprintln!(
        "[INFO] map build mode: {}",
        if map_final.exists() && !map_tmp.exists() {
            "reused existing map (no temp file)"
        } else {
            "rebuilt from VCF (temp existed until rename)"
        }
    );

    // Prepare output (streaming; we only reach here AFTER the map is finalized)
    let out_path = format!("{}.renamed.tsv", qtl_path);
    let mut qtl_reader = BufReader::new(File::open(qtl_path)?);
    let writer = Arc::new(Mutex::new(BufWriter::new(File::create(&out_path)?)));
    let replaced_ctr = Arc::new(AtomicUsize::new(0));
    let unchanged_ctr = Arc::new(AtomicUsize::new(0));

    // Detect and write header unchanged, if present
    let mut wrote_header = false;
    let mut first_line = String::new();
    if qtl_reader.read_line(&mut first_line)? > 0 {
        let trimmed = first_line.trim_end_matches(['\n', '\r']).to_string();
        if is_probable_header_tsv(&trimmed) {
            wrote_header = true;
            eprintln!("[INFO] QTL header detected; preserving first line as-is");
            let mut w = writer.lock().unwrap();
            writeln!(w, "{}", trimmed)?;
        } else if !trimmed.is_empty() {
            // process this first data row immediately as a single-item chunk
            #[cfg(feature = "rayon")]
            {
                let w = Arc::clone(&writer);
                let map_ref = &id_key_map;
                let rep = Arc::clone(&replaced_ctr);
                let unc = Arc::clone(&unchanged_ctr);
                [trimmed].par_iter().for_each(|row| {
                    let out = replace_col2_with_map(row, map_ref);
                    let mut wlock = w.lock().unwrap();
                    writeln!(wlock, "{}", out.line).unwrap();
                    if out.changed {
                        rep.fetch_add(1, Ordering::Relaxed);
                    } else {
                        unc.fetch_add(1, Ordering::Relaxed);
                    }
                });
            }
            #[cfg(not(feature = "rayon"))]
            {
                let out = replace_col2_with_map(&trimmed, &id_key_map);
                let mut wlock = writer.lock().unwrap();
                writeln!(wlock, "{}", out.line)?;
                if out.changed {
                    replaced_ctr.fetch_add(1, Ordering::Relaxed);
                } else {
                    unchanged_ctr.fetch_add(1, Ordering::Relaxed);
                }
            }
        }
    }

    // QTL replacement streaming with par_bridge
    use std::io::BufRead;
    let mut _total_rows: usize = 0;
    eprintln!("[INFO] Starting QTL replacement stream...");
    #[cfg(feature = "rayon")]
    {
        use rayon::iter::ParallelBridge;
        let w = Arc::clone(&writer);
        let map_ref = &id_key_map;
        let rep = Arc::clone(&replaced_ctr);
        let unc = Arc::clone(&unchanged_ctr);
        let lines = qtl_reader.lines();
        lines.par_bridge().for_each(|res| {
            if let Ok(line) = res {
                if line.trim().is_empty() {
                    return;
                }
                let replaced = replace_col2_with_map(&line, map_ref);
                let mut lock = w.lock().unwrap();
                let _ = writeln!(lock, "{}", replaced.line);
                if replaced.changed {
                    rep.fetch_add(1, Ordering::Relaxed);
                } else {
                    unc.fetch_add(1, Ordering::Relaxed);
                }
            }
        });
        _total_rows = replaced_ctr.load(Ordering::Relaxed) + unchanged_ctr.load(Ordering::Relaxed);
    }
    #[cfg(not(feature = "rayon"))]
    {
        for l in qtl_reader.lines() {
            let s = l?;
            if s.trim().is_empty() {
                continue;
            }
            let out = replace_col2_with_map(&s, &id_key_map);
            let mut wlock = writer.lock().unwrap();
            writeln!(wlock, "{}", out.line)?;
            if out.changed {
                replaced_ctr.fetch_add(1, Ordering::Relaxed);
            } else {
                unchanged_ctr.fetch_add(1, Ordering::Relaxed);
            }
            total_rows += 1;
        }
    }

    // Finish
    writer.lock().unwrap().flush()?;
    eprintln!(
        "[INFO] Replacement done: {} lines processed",
        replaced_ctr.load(Ordering::Relaxed) + unchanged_ctr.load(Ordering::Relaxed)
    );
    eprintln!(
        "[INFO] QTL replaced: rows={} replaced={} unchanged={} header_written={} output={}",
        replaced_ctr.load(Ordering::Relaxed) + unchanged_ctr.load(Ordering::Relaxed),
        replaced_ctr.load(Ordering::Relaxed),
        unchanged_ctr.load(Ordering::Relaxed),
        wrote_header,
        out_path
    );
    eprintln!("[INFO] total elapsed: {:.2?}", t0.elapsed());
    Ok(())
}

struct ReplaceOut {
    line: String,
    changed: bool,
}

fn replace_col2_with_map(line: &str, id_key_map: &HashMap<String, String>) -> ReplaceOut {
    let mut cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 2 {
        return ReplaceOut {
            line: line.to_string(),
            changed: false,
        };
    }
    let col2 = cols[1].trim();

    // 1) Try exact token in map
    if let Some(norm) = id_key_map.get(col2) {
        cols[1] = norm;
        return ReplaceOut {
            line: cols.join("\t"),
            changed: true,
        };
    }

    // 2) Try normalized token (e.g. strip 'chr', normalize MT/M)
    if let Some(norm2) = try_normalize_token(col2) {
        if let Some(norm) = id_key_map.get(&norm2) {
            cols[1] = norm;
            return ReplaceOut {
                line: cols.join("\t"),
                changed: true,
            };
        }
    }

    ReplaceOut {
        line: line.to_string(),
        changed: false,
    }
}
