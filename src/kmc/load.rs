//! Loads a stranded KMC database into myloasm's split k-mer table.

use super::reader::KmcReader;
use crate::seq_parse::split_kmer_passes_count_thresholds;
use crate::types::{Kmer64, BYTE_TO_SEQ};
use std::path::Path;

/// Reads `db` (built by `myloasm-kmc`) into `(canonical k-mer, [reverse count, forward count])`
/// pairs, the exact table `seq_parse::read_to_split_kmers` produces:
///
/// * the canonical k-mer is the one whose split form (middle base masked) is smaller. For every
///   k-mer that is not a split-palindrome this coincides with KMC's canonical form, but the choice
///   is recomputed here rather than assumed;
/// * split-palindromes (`split(kmer) == split(revcomp(kmer))`) are dropped, as in the in-memory path;
/// * the count thresholds are re-applied, in case the database was built with looser ones.
pub fn split_kmers_from_db(db: &Path, k: usize) -> Result<Vec<(Kmer64, [u32; 2])>, String> {
    let mut reader = KmcReader::open(db)
        .map_err(|e| format!("could not open KMC database {}: {e}", db.display()))?;
    let info = reader.info().clone();
    if !info.is_stranded() {
        return Err(format!(
            "{} was not built with per-strand counters; build it with myloasm-kmc (or kmc -sc)",
            db.display()
        ));
    }
    if info.kmer_length as usize != k {
        return Err(format!(
            "{} was built with k = {} but myloasm is running with k = {k}",
            db.display(),
            info.kmer_length
        ));
    }
    log::info!(
        "Loading {} k-mers from {} (k = {}, database thresholds: total >= {}, per strand >= {})",
        info.total_kmers,
        db.display(),
        info.kmer_length,
        info.min_count,
        info.min_count_per_strand
    );

    let split_mask: u64 = !(3u64 << (k - 1));
    let mut kmers = Vec::with_capacity(info.total_kmers as usize);
    let mut palindromes = 0u64;
    let mut below_thresholds = 0u64;
    while let Some(record) = reader
        .next_kmer()
        .map_err(|e| format!("error reading {}: {e}", db.display()))?
    {
        let forward = encode(record.kmer);
        let reverse = reverse_complement(forward, k);
        let (split_f, split_r) = (forward & split_mask, reverse & split_mask);
        if split_f == split_r {
            palindromes += 1;
            continue;
        }
        // counts[1]: seen as the canonical k-mer itself, counts[0]: seen as its reverse complement.
        let as_itself = saturate(record.count_fwd);
        let as_revcomp = saturate(record.count_rev);
        let (canonical, counts) = if split_f < split_r {
            (forward, [as_revcomp, as_itself])
        } else {
            (reverse, [as_itself, as_revcomp])
        };
        if !split_kmer_passes_count_thresholds(counts) {
            below_thresholds += 1;
            continue;
        }
        kmers.push((canonical, counts));
    }
    log::info!(
        "Loaded {} k-mers ({} split-palindromes and {} k-mers below the count thresholds skipped)",
        kmers.len(),
        palindromes,
        below_thresholds
    );
    Ok(kmers)
}

/// Deletes both files of a KMC database. Returns whether any file was removed; missing files are
/// not an error, other failures are logged.
pub fn remove_db(db: &Path) -> bool {
    let mut removed = false;
    for extension in [".kmc_pre", ".kmc_suf"] {
        let mut path = db.as_os_str().to_owned();
        path.push(extension);
        match std::fs::remove_file(&path) {
            Ok(()) => removed = true,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => log::warn!("Could not remove {}: {}", Path::new(&path).display(), e),
        }
    }
    removed
}

/// Packs an ACGT string into 2 bits per base, first base in the most significant position.
#[inline]
fn encode(kmer: &[u8]) -> u64 {
    kmer.iter()
        .fold(0u64, |acc, &b| (acc << 2) | BYTE_TO_SEQ[b as usize] as u64)
}

#[inline]
fn reverse_complement(mut kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (kmer & 3));
        kmer >>= 2;
    }
    rc
}

#[inline]
fn saturate(count: u64) -> u32 {
    count.min(u32::MAX as u64) as u32
}
