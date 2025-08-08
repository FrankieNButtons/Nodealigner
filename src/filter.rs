use crate::io_utils::VcfRecord;
use std::collections::{HashMap, HashSet};
use rayon::prelude::*;

/// Parallel filtering of records by keywords in path (skip if any keyword is present)
pub fn filter_records_by_path(
    records: Vec<VcfRecord>,
    node2path: &HashMap<u64, String>,
    skip_keywords: &HashSet<String>,
) -> Vec<VcfRecord> {
    records
        .into_par_iter()
        .filter(|rec| {
            if let Some(path) = node2path.get(&rec.pos) {
                !skip_keywords.iter().any(|kw| path.contains(kw))
            } else {
                true
            }
        })
        .collect()
}