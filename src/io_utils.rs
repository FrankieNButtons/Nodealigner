use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

#[allow(dead_code)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub rest: Vec<String>,
    pub raw: String, // Full original line, for output
}

/// Read VCF file, returning (header_lines, records)
pub fn read_vcf_records(
    path: &str,
) -> Result<(Vec<String>, Vec<VcfRecord>), Box<dyn std::error::Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);

    let mut headers = Vec::new();
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            headers.push(line);
        } else {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                continue;
            }
            let chrom = fields[0].to_string();
            let pos: u64 = fields[1].parse().unwrap_or(0);
            let rest = fields[2..].iter().map(|s| s.to_string()).collect();
            records.push(VcfRecord {
                chrom,
                pos,
                rest,
                raw: line,
            });
        }
    }
    Ok((headers, records))
}

/// Write VCF file
pub fn write_vcf_records(
    path: &str,
    headers: &[String],
    records: &[VcfRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    let f = File::create(path)?;
    let mut writer = BufWriter::new(f);

    for h in headers {
        writeln!(writer, "{h}")?;
    }
    for rec in records {
        writeln!(writer, "{}", rec.raw)?;
    }
    Ok(())
}


/// 读取 alignment.tsv，建立 node->path（首行为表头自动跳过，首列node，末列path，自动trim）
pub fn read_tsv_node_path(
    path: &str,
) -> Result<HashMap<u64, String>, Box<dyn std::error::Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    let mut map = HashMap::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        if idx == 0 || line.trim().is_empty() { continue; } // 跳过表头
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 { continue; }
        // debug打印node=1相关内容
        if fields[0].trim() == "1" {
            // println!("[debug][node==1] idx={idx}, fields={:?}", fields);
        }
        let node: u64 = fields[0].trim().parse().unwrap_or(0);
        let path = fields[4].trim().to_string();
        map.insert(node, path);
    }
    // println!("[debug] node2path contains key 1? {}", map.contains_key(&1));
    Ok(map)
}

/// 读取 reference.tsv（extract产物），建立 node->path（跳表头，node在首列，path在第四列）
pub fn read_reference_tsv(
    path: &str,
) -> Result<HashMap<u64, String>, Box<dyn std::error::Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    let mut map = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 || line.trim().is_empty() { continue; } // 跳表头和空行
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 { continue; }
        let node: u64 = fields[0].trim().parse().unwrap_or(0);
        let path = fields[3].trim().to_string();
        map.insert(node, path);
    }
    Ok(map)
}