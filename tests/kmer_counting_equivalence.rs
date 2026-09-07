//! The on-disk k-mer counter (`myloasm-kmc-v1`, an external program wrapping a KMC fork) must
//! produce exactly the split k-mer table that myloasm's in-memory counter produces: the same
//! canonical k-mers with the same `[reverse count, forward count]` pairs, after the same
//! palindrome rule, count thresholds and middle-base quality filter.
//!
//! The binary is found through `MYLOASM_KMC_BIN`, or next to the test executable / on `PATH`
//! like myloasm itself does. When it is not installed these tests print a notice and pass
//! without checking anything.

use myloasm::cli::Cli;
use myloasm::kmc;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const K: usize = 21;
const THREADS: usize = 4;

type SplitKmerTable = Vec<(u64, [u32; 2])>;

/// `bloom_filter_gb == 0.0` disables the bloom-filter prepass and counts every k-mer exactly.
fn count_in_memory(reads: &str, bloom_filter_gb: f64) -> SplitKmerTable {
    let args = Cli {
        input_files: vec![reads.to_string()],
        kmer_size: K,
        threads: THREADS,
        bloom_filter_size: Some(bloom_filter_gb),
        ..Default::default()
    };
    let mut table = myloasm::seq_parse::read_to_split_kmers(K, THREADS, &args);
    table.sort_unstable();
    table
}

/// The counter binary, or `None` (after printing why) if it is not installed.
fn counter_binary() -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("MYLOASM_KMC_BIN") {
        return Some(PathBuf::from(path));
    }
    match kmc::runner::locate_binary() {
        Ok(path) => Some(path),
        Err(message) => {
            eprintln!("SKIPPED: {message}");
            None
        }
    }
}

fn count_on_disk(binary: &Path, reads: &str) -> SplitKmerTable {
    let work_dir = TempDir::new().unwrap();
    let db = work_dir.path().join("counts");
    let tmp_dir = work_dir.path().join("kmc_tmp");
    std::fs::create_dir(&tmp_dir).unwrap();

    let config = kmc::runner::DiskCountConfig {
        input_files: &[reads.to_string()],
        kmer_size: K,
        threads: THREADS,
        ram_gb: 2,
        output_db: &db,
        tmp_dir: &tmp_dir,
        binary: Some(binary),
    };
    kmc::runner::count_kmers_on_disk(&config).expect("myloasm-kmc-v1 should succeed");

    let mut table = kmc::load::split_kmers_from_db(&db, K).unwrap();
    table.sort_unstable();
    table
}

fn assert_same_table(in_memory: &SplitKmerTable, on_disk: &SplitKmerTable) {
    assert!(
        in_memory.len() > 10_000,
        "test data too small to be meaningful"
    );
    assert_eq!(in_memory.len(), on_disk.len(), "different number of k-mers");
    for (mem, disk) in in_memory.iter().zip(on_disk) {
        assert_eq!(mem, disk, "k-mer or per-strand counts differ");
    }
}

#[test]
fn fastq_with_base_qualities() {
    // Real nanopore reads: exercises the middle-base quality filter on both sides.
    let Some(binary) = counter_binary() else {
        return;
    };
    let reads = "tests/reads/40kb_plas.fq";
    assert_same_table(
        &count_in_memory(reads, 0.05),
        &count_on_disk(&binary, reads),
    );
}

#[test]
fn fastq_without_bloom_filter_prepass() {
    let Some(binary) = counter_binary() else {
        return;
    };
    let reads = "tests/reads/40kb_plas.fq";
    assert_same_table(&count_in_memory(reads, 0.0), &count_on_disk(&binary, reads));
}

#[test]
fn fasta_without_qualities() {
    let Some(binary) = counter_binary() else {
        return;
    };
    let reads = "tests/reads/40kb_plas.fa";
    assert_same_table(
        &count_in_memory(reads, 0.05),
        &count_on_disk(&binary, reads),
    );
}

#[test]
fn quality_filter_actually_removes_kmers() {
    // Sanity check that the FASTQ tests are not vacuous: the same reads without their qualities
    // keep more k-mers, on both sides.
    let Some(binary) = counter_binary() else {
        return;
    };
    let fastq = "tests/reads/40kb_plas.fq";
    let work_dir = TempDir::new().unwrap();
    let fasta = work_dir.path().join("reads.fa");
    fastq_to_fasta(fastq, &fasta);
    let fasta = fasta.to_str().unwrap();

    let with_qualities = count_in_memory(fastq, 0.0);
    let without_qualities = count_in_memory(fasta, 0.0);
    assert!(with_qualities.len() < without_qualities.len());
    assert_same_table(&without_qualities, &count_on_disk(&binary, fasta));
}

fn fastq_to_fasta(fastq: &str, fasta: &Path) {
    let text = std::fs::read_to_string(fastq).unwrap();
    let mut out = String::new();
    for record in text.lines().collect::<Vec<_>>().chunks(4) {
        out.push('>');
        out.push_str(&record[0][1..]);
        out.push('\n');
        out.push_str(record[1]);
        out.push('\n');
    }
    std::fs::write(fasta, out).unwrap();
}
