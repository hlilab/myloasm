//! On-disk k-mer counting with the external `myloasm-kmc-v1` program.
//!
//! myloasm's default counter (`seq_parse::read_to_split_kmers`) holds every k-mer in memory. For
//! very large datasets this module instead runs `myloasm-kmc-v1` (a separately installed wrapper
//! around a KMC fork with per-strand counters, https://github.com/bluenote-1577/myloasm-kmc),
//! which counts canonical k-mers on disk with one counter per strand, and then streams the
//! resulting database back into the same split k-mer table:
//!
//! * [`runner`] locates and runs `myloasm-kmc-v1`;
//! * [`reader`] parses KMC databases, including the per-strand counter layout;
//! * [`load`] turns a stranded database into `(k-mer, [reverse count, forward count])` pairs,
//!   applying the same canonical form, palindrome rule and count thresholds as the in-memory path.
//!
//! Both counters give identical tables for the same reads (checked by
//! `tests/kmer_counting_equivalence.rs` whenever the binary is installed), with one exception:
//! KMC skips k-mers containing non-ACGT symbols, while the in-memory counter reads them as `A`.

pub mod load;
pub mod reader;
pub mod runner;

use crate::cli::Cli;
use crate::types::Kmer64;
use std::path::{Path, PathBuf};

/// Base name of the database that [`split_kmers_via_disk_count`] writes under `binary_temp/`.
pub const DISK_COUNT_DB_NAME: &str = "kmer_counts_kmc";

/// Where `--kmc` writes its database. The database can be large, so it is kept only until the
/// SNPmer checkpoint has been written (see [`remove_disk_count_db`]); until then a run that dies
/// can be restarted from it with `--kmc-stranded-db`.
pub fn disk_count_db_path(output_dir: &Path) -> PathBuf {
    output_dir.join("binary_temp").join(DISK_COUNT_DB_NAME)
}

/// Deletes the database written by `--kmc`, if it exists.
pub fn remove_disk_count_db(output_dir: &Path) {
    if load::remove_db(&disk_count_db_path(output_dir)) {
        log::info!(
            "Removed the k-mer count database; it is not needed once SNPmers are extracted."
        );
    }
}

/// Counts the split k-mers of `args.input_files` on disk and loads them.
pub fn split_kmers_via_disk_count(
    args: &Cli,
    output_dir: &Path,
) -> Result<Vec<(Kmer64, [u32; 2])>, String> {
    let db = disk_count_db_path(output_dir);
    let tmp_dir = output_dir.join("kmc_tmp");
    std::fs::create_dir_all(&tmp_dir)
        .map_err(|e| format!("could not create {}: {e}", tmp_dir.display()))?;

    let config = runner::DiskCountConfig {
        input_files: &args.input_files,
        kmer_size: args.kmer_size,
        threads: args.threads,
        ram_gb: args.kmc_ram,
        output_db: &db,
        tmp_dir: &tmp_dir,
        binary: None,
    };
    let counted = runner::count_kmers_on_disk(&config);
    // KMC deletes its own temporary files; the directory (and anything left after a failure) goes too.
    if let Err(e) = std::fs::remove_dir_all(&tmp_dir) {
        log::warn!("Could not remove {}: {}", tmp_dir.display(), e);
    }
    counted?;

    load::split_kmers_from_db(&db, args.kmer_size)
}
