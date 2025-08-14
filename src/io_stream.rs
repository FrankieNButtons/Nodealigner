use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::OnceLock;
use gfa_reader::Gfa;

static REF_NODE2PATH: OnceLock<HashMap<u64, String>> = OnceLock::new();
static REF_NODE2START: OnceLock<HashMap<u64, u64>> = OnceLock::new();
static REF_NODE2SEQ: OnceLock<HashMap<u64, String>> = OnceLock::new();

// #[allow(dead_code)]
// pub struct VcfRecord {
//     pub chrom: String,
//     pub pos: u64,
//     pub rest: Vec<String>,
//     pub raw: String, // Full original line, for output
// }

// /// Read VCF file, returning (header_lines, records)
// #[allow(dead_code)]
// pub fn read_vcf_records(
//     path: &str,
// ) -> Result<(Vec<String>, Vec<VcfRecord>), Box<dyn std::error::Error>> {
//     let f = File::open(path)?;
//     let reader = BufReader::new(f);

//     let mut headers = Vec::new();
//     let mut records = Vec::new();

//     for line in reader.lines() {
//         let line = line?;
//         if line.starts_with('#') {
//             headers.push(line);
//         } else {
//             let fields: Vec<&str> = line.split('\t').collect();
//             if fields.len() < 2 {
//                 continue;
//             }
//             let chrom = fields[0].to_string();
//             let pos: u64 = fields[1].parse().unwrap_or(0);
//             let rest = fields[2..].iter().map(|s| s.to_string()).collect();
//             records.push(VcfRecord {
//                 chrom,
//                 pos,
//                 rest,
//                 raw: line,
//             });
//         }
//     }
//     Ok((headers, records))
// }

// /// Write VCF file
// #[allow(dead_code)]
// pub fn write_vcf_records(
//     path: &str,
//     headers: &[String],
//     records: &[VcfRecord],
// ) -> Result<(), Box<dyn std::error::Error>> {
//     let f = File::create(path)?;
//     let mut writer = BufWriter::new(f);

//     for h in headers {
//         writeln!(writer, "{h}")?;
//     }
//     for rec in records {
//         writeln!(writer, "{}", rec.raw)?;
//     }
//     Ok(())
// }

/// Streaming stats for VCF transform
#[derive(Default, Debug)]
pub struct StreamStats {
    pub total: u64,
    pub replaced: u64,
    pub skipped: u64,
    pub unmapped: u64,
    pub replaced_chrom: u64,
    pub replaced_pos: u64,
    pub replaced_id: u64,
    pub replaced_ref: u64,
    pub missing_start: u64,
    pub missing_seq: u64,
    pub used_ref_map: u64,
    pub used_aln_map: u64,
}

#[derive(Debug, Clone, Default)]
pub struct AlnInfo {
    pub path: String,
    pub distance: i64,
    pub position: u64,
}

/// Return true if the raw CHROM field should be skipped entirely (substring match)
pub fn should_skip_chrom(chrom: &str, skip: &HashSet<String>) -> bool {
    if skip.is_empty() {
        return false;
    }
    // case-sensitive substring match by default
    for k in skip {
        if !k.is_empty() && chrom.contains(k) {
            return true;
        }
    }
    false
}

/// Parse a node id from a CHROM string. Accepts either a plain integer (e.g., "1234")
/// or a suffix-digit pattern (e.g., "node_1234"). Returns None if no digits found.
fn parse_node_id_from_chrom(chrom: &str) -> Option<u64> {
    // fast path: pure integer
    if let Ok(v) = chrom.parse::<u64>() {
        return Some(v);
    }
    // fallback: extract last continuous digit run
    let mut digits = String::new();
    for c in chrom.chars().rev() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else if !digits.is_empty() {
            break;
        }
    }
    if digits.is_empty() {
        return None;
    }
    let digits_fwd: String = digits.chars().rev().collect();
    digits_fwd.parse::<u64>().ok()
}

/// Case-insensitive search for "chr" and extract a canonical token:
/// - Looks for the *last* occurrence of "chr" (case-insensitive) so
///   strings like "GRCh38.chr12_random" work.
/// - Captures immediately following token consisting of [0-9]+ or X/Y/M (case-insensitive).
/// - Returns (token_upper, has_suffix_after_token, found_chr)
fn extract_chr_token(raw: &str) -> (Option<String>, bool, bool) {
    let lower = raw.to_ascii_lowercase();
    if let Some(idx) = lower.rfind("chr") {
        let tail = &raw[idx + 3..];
        // Walk tail to capture token consisting of digits or X/Y/M (case-insensitive)
        let mut token = String::new();
        for ch in tail.chars() {
            let c = ch.to_ascii_uppercase();
            if c.is_ascii_digit() || c == 'X' || c == 'Y' || c == 'M' {
                token.push(c);
            } else {
                break;
            }
        }
        if token.is_empty() {
            return (None, !tail.is_empty(), true);
        }
        // Whether there is anything *after* the token and before hitting end
        let has_suffix = tail.len() > token.len();
        return (Some(token), has_suffix, true);
    }
    (None, false, false)
}

/// Return whether a token is a standard human chromosome label:
/// 1..22 or X/Y/M (token must be uppercase, digits allowed).
fn is_std_human_chr_token(token_upper: &str) -> bool {
    match token_upper {
        "X" | "Y" | "M" => true,
        _ => {
            token_upper.chars().all(|c| c.is_ascii_digit())
                && token_upper
                    .parse::<u32>()
                    .map(|n| (1..=22).contains(&n))
                    .unwrap_or(false)
        }
    }
}

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
                    // digits or X/Y/M (not necessarily standard range at this level)
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
            Some(format!("chr{}", t))
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
            Some(t)
        }

        _ => Some(raw.to_string()),
    }
}

/// Replacement strategy (new):
/// - 丢弃含有 `--skip` 关键词（子串匹配）的记录。
/// - 解析 node_id：优先从 CHROM 提取纯数字；否则使用原始 POS 值。
/// - 若能在 alignment.tsv 找到该 node：
///     * CHROM := 规范化后的 path 名（按 ignore_level 过滤/归一化）。
///     * ID    := 原始 POS 字符串（未替换前）。
///     * POS   := distance + position + 1（alignment.tsv）。
///     * REF   := 来自 GFA 的段序列；若未提供 GFA，则回退 reference.tsv 的 seq。
/// - 若 alignment.tsv 缺失该 node：
///     * 仍可用 reference.tsv 的 path/start/seq（若存在）进行回退；否则按 ignore 规则仅重写/保留 CHROM。
/// - 若提供 GFA，则序列优先来自 GFA，且 reference.tsv 的序列将被忽略。
pub fn stream_replace_chrom_to_tmp(
    vcf_path: &str,
    tmp_out_path: &str,
    node2aln: &HashMap<u64, AlnInfo>,
    skip: &HashSet<String>,
    ignore_level: u8, // 0..=5
    gfa: Option<&Gfa<u32, (), ()>>, // 如果提供GFA，优先用其序列并可忽略reference.tsv
) -> Result<StreamStats, Box<dyn std::error::Error>> {
    let f_in = File::open(vcf_path)?;
    let reader = BufReader::new(f_in);

    let f_out = File::create(tmp_out_path)?;
    let mut writer = BufWriter::new(f_out);

    let mut stats = StreamStats::default();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            writeln!(writer, "{line}")?;
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }

        let raw_chrom = fields[0];
        let skip_now = should_skip_chrom(raw_chrom, skip);

        let chrom_node_id_opt = parse_node_id_from_chrom(raw_chrom);
        let mut pos_node_id_opt: Option<u64> = None;
        if chrom_node_id_opt.is_none() {
            if let Ok(pos_as_u64) = fields[1].trim().parse::<u64>() {
                pos_node_id_opt = Some(pos_as_u64);
            }
        }
        let node_id_opt = chrom_node_id_opt.or(pos_node_id_opt);

        if skip_now {
            stats.skipped += 1;
            continue;
        }

        stats.total += 1;
        let mut wrote = false;

        if let Some(node_id) = node_id_opt {
            let path_from_ref = REF_NODE2PATH.get().and_then(|m| m.get(&node_id));
            let seq_from_ref = REF_NODE2SEQ.get().and_then(|m| m.get(&node_id));
            let aln_info = node2aln.get(&node_id);

            // 选路径：alignment.tsv 优先；否则回退 reference.tsv
            let chosen_path_opt: Option<&str> = aln_info
                .map(|a| a.path.as_str())
                .or(path_from_ref.map(|s| s.as_str()));

            if let Some(path_val) = chosen_path_opt {
                if let Some(norm_chr) = apply_ignore_rules(path_val, ignore_level) {
                    let mut out_fields: Vec<String> = fields.iter().map(|s| s.to_string()).collect();

                    // 保存原始 POS 到 ID（如果有第三列）
                    let orig_pos = out_fields.get(1).cloned().unwrap_or_default();
                    if out_fields.len() >= 3 {
                        out_fields[2] = orig_pos.clone();
                        stats.replaced_id += 1;
                    }
                    // 1) CHROM 由规范化的 path 名得到
                    out_fields[0] = norm_chr.clone();
                    stats.replaced_chrom += 1;

                    // 2) REF：优先从 GFA 取段序列；若无 GFA，则回退到 reference.tsv 的 seq
                    let ref_set = false;
                    if let Some(g) = gfa {
                        let nid_u32 = node_id as u32;
                        if nid_u32 as usize >= g.get_index_low() && nid_u32 as usize <= g.get_index_high() {
                            let seq = g.get_sequence_by_id(&nid_u32);
                            if out_fields.len() >= 4 {
                                out_fields[3] = seq.to_string();
                                stats.replaced_ref += 1;
                            }
                        }
                    }
                    if !ref_set {
                        if let Some(seq) = seq_from_ref {
                            if out_fields.len() >= 4 {
                                out_fields[3] = seq.to_string();
                                stats.replaced_ref += 1;
                            }
                        } else {
                            stats.missing_seq += 1;
                        }
                    }

                    // 3) POS：来自 alignment.tsv 的 distance + position + 1；若缺失则回退到原始 POS 或 reference start
                    let mut pos_set = false;
                    if let Some(a) = aln_info {
                        // POS calculation rules:
                        // - If distance is extremely large (> 1_000_000_000), use position directly.
                        // - Otherwise, compute (distance + 1) as unsigned offset and add to position.
                        //   Special case: distance == -1 => offset = 0 (same node), so POS = position.
                        //   Any distance < -1 is clamped to behave like -1 to avoid u64 wrap.
                        let new_pos: u64 = if a.distance > 1_000_000_000 {
                            a.position
                        } else {
                            let offset = if a.distance >= -1 {
                                (a.distance + 1) as u64
                            } else {
                                0u64
                            };
                            offset.saturating_add(a.position)
                        };
                        if out_fields.len() >= 2 {
                            out_fields[1] = new_pos.to_string();
                            stats.replaced_pos += 1;
                            pos_set = true;
                        }
                        stats.used_aln_map += 1;
                    }
                    if !pos_set {
                        if let Some(start_val) = REF_NODE2START.get().and_then(|m| m.get(&node_id)).copied() {
                            if out_fields.len() >= 2 {
                                out_fields[1] = start_val.to_string();
                                stats.replaced_pos += 1;
                                pos_set = true;
                            }
                            stats.used_ref_map += 1;
                        }
                    }
                    if !pos_set { stats.missing_start += 1; }

                    writeln!(writer, "{}", out_fields.join("\t"))?;
                    stats.replaced += 1;
                    wrote = true;
                } else {
                    stats.skipped += 1;
                    wrote = true;
                }
            }
        }

        if !wrote {
            match apply_ignore_rules(fields[0], ignore_level) {
                Some(norm_chr) => {
                    let mut out_line = String::new();
                    out_line.push_str(&norm_chr);
                    for f in &fields[1..] {
                        out_line.push('\t');
                        out_line.push_str(f);
                    }
                    writeln!(writer, "{out_line}")?;
                    stats.unmapped += 1;
                }
                _none => {
                    stats.skipped += 1;
                }
            }
        }
    }

    eprintln!(
        "[summary] total={} replaced={} unmapped={} skipped={} chrom={} pos={} id={} ref={} missing_start={} missing_seq={} via_ref={} via_aln={}",
        stats.total,
        stats.replaced,
        stats.unmapped,
        stats.skipped,
        stats.replaced_chrom,
        stats.replaced_pos,
        stats.replaced_id,
        stats.replaced_ref,
        stats.missing_start,
        stats.missing_seq,
        stats.used_ref_map,
        stats.used_aln_map
    );

    writer.flush()?;
    Ok(stats)
}

/// 读取 alignment.tsv（首行为表头或数据）。
/// 需要列：node、path，以及可选列 distance、position；
/// - 自动探测列名（不区分大小写）。
/// - 若无表头，则按默认索引解析：node(0), distance(1), position(2), path(4)。
pub fn read_alignment_tsv(path: &str) -> Result<HashMap<u64, AlnInfo>, Box<dyn std::error::Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);

    let mut map: HashMap<u64, AlnInfo> = HashMap::new();

    // Column indices (discovered from header if present)
    let mut idx_node: Option<usize> = None;
    let mut idx_path: Option<usize> = None;
    let mut idx_distance: Option<usize> = None;
    let mut idx_position: Option<usize> = None;

    let mut header_seen = false;

    for (_line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        // Detect header by typical tokens in the first non-empty line
        if !header_seen {
            let l0 = fields[0].to_ascii_lowercase();
            if l0.contains("node") || l0.contains("id") {
                header_seen = true;
                for (i, col) in fields.iter().enumerate() {
                    let c = col.to_ascii_lowercase();
                    if idx_node.is_none() && (c == "node" || c == "id" || c == "segment" || c == "seg") { idx_node = Some(i); }
                    if idx_path.is_none() && (c == "path" || c == "chrom" || c == "name") { idx_path = Some(i); }
                    if idx_distance.is_none() && (c == "distance" || c == "dist" || c == "offset") { idx_distance = Some(i); }
                    if idx_position.is_none() && (c == "position" || c == "pos") { idx_position = Some(i); }
                }
                // Fall through to next line for data
                continue;
            }
        }

        // If we got here and we still don't have indices, assume default positional layout
        let node_i = idx_node.unwrap_or(0);
        let dist_i = idx_distance.unwrap_or(1);
        let pos_i  = idx_position.unwrap_or(2);
        let path_i = idx_path.unwrap_or_else(|| fields.len().saturating_sub(1).max(4));

        if fields.len() <= node_i { continue; }
        let node_raw = fields[node_i].trim();
        if node_raw.is_empty() { continue; }
        let node: u64 = match node_raw.parse() { Ok(v) => v, Err(_) => continue };

        let mut info = AlnInfo::default();
        // path
        if fields.len() > path_i { info.path = fields[path_i].trim().to_string(); }
        // distance
        if fields.len() > dist_i { info.distance = fields[dist_i].trim().parse().unwrap_or(0); }
        // position
        if fields.len() > pos_i { info.position = fields[pos_i].trim().parse().unwrap_or(0); }

        map.insert(node, info);
    }

    Ok(map)
}

/// 读取 reference.tsv（extract产物），支持4列和6列格式：
/// 4列: node, start, end, path
/// 6列: node, start, end, seq, length, path
/// - 返回 node->path（用于向后兼容）
/// - 同时写入全局 OnceLock：node->start, node->seq, node->path（供流式替换使用）
pub fn read_reference_tsv(path: &str) -> Result<HashMap<u64, String>, Box<dyn std::error::Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);

    let mut map_path: HashMap<u64, String> = HashMap::new();
    let mut map_start: HashMap<u64, u64> = HashMap::new();
    let mut map_seq: HashMap<u64, String> = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        // 跳过表头（以 "node\t" 开头）
        if i == 0 && trimmed.to_ascii_lowercase().starts_with("node\t") {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }

        // node 在首列
        let node: u64 = match fields[0].trim().parse::<u64>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // start 在第2列（索引1）
        if let Ok(s) = fields.get(1).unwrap_or(&"").trim().parse::<u64>() {
            map_start.insert(node, s);
        }

        // path 在最后一列
        let path_col = fields.len() - 1;
        let path_val = fields[path_col].trim().to_string();
        if !path_val.is_empty() {
            map_path.insert(node, path_val);
        }

        // seq 在6列表时的第4列（索引3）
        if fields.len() >= 6 {
            let seq_val = fields[3].trim().to_string();
            if !seq_val.is_empty() {
                map_seq.insert(node, seq_val);
            }
        }
    }

    // 写入 OnceLock（若已存在则保持原样）
    let _ = REF_NODE2PATH.set(map_path.clone());
    let _ = REF_NODE2START.set(map_start);
    let _ = REF_NODE2SEQ.set(map_seq);

    Ok(map_path)
}
