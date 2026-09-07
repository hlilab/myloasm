//! Runs the external `myloasm-kmc-v1` program, which wraps a KMC fork with per-strand counters.
//!
//! The counter is a separate, separately installed executable: KMC is GPL-3 and needs a C++14
//! toolchain, neither of which myloasm itself takes on. The `v1` in the name is the interface
//! version (command-line flags and database format) that this module was written against; a
//! binary with a different suffix is not used.

use crate::constants::{
    MID_BASE_THRESHOLD_INITIAL, MIN_SPLIT_KMER_COUNT_PER_STRAND, MIN_SPLIT_KMER_COUNT_TOTAL,
};
use std::collections::VecDeque;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub const BINARY_NAME: &str = "myloasm-kmc-v1";

/// How many trailing lines of the counter's output are quoted when it fails.
const ERROR_CONTEXT_LINES: usize = 20;

pub struct DiskCountConfig<'a> {
    pub input_files: &'a [String],
    pub kmer_size: usize,
    pub threads: usize,
    pub ram_gb: usize,
    /// Database base path to write (`.kmc_pre` / `.kmc_suf` get appended).
    pub output_db: &'a Path,
    /// Existing directory for KMC's temporary files.
    pub tmp_dir: &'a Path,
    /// Explicit path to the counter; `None` looks next to the running executable, then on `PATH`.
    pub binary: Option<&'a Path>,
}

/// Counts the k-mers of `config.input_files` into `config.output_db` with the count thresholds
/// and middle-base quality filter of the in-memory counter. The counter's own log lines are
/// forwarded at debug level.
pub fn count_kmers_on_disk(config: &DiskCountConfig) -> Result<(), String> {
    let binary = match config.binary {
        Some(path) => path.to_path_buf(),
        None => locate_binary()?,
    };

    let mut command = Command::new(&binary);
    command
        .arg("--output")
        .arg(config.output_db)
        .arg("--tmp-dir")
        .arg(config.tmp_dir)
        .arg("--kmer-size")
        .arg(config.kmer_size.to_string())
        .arg("--threads")
        .arg(config.threads.max(1).to_string())
        .arg("--ram-gb")
        .arg(config.ram_gb.max(1).to_string())
        .arg("--min-count")
        .arg(MIN_SPLIT_KMER_COUNT_TOTAL.to_string())
        .arg("--min-count-per-strand")
        .arg(MIN_SPLIT_KMER_COUNT_PER_STRAND.to_string())
        .arg("--min-mid-quality")
        .arg(MID_BASE_THRESHOLD_INITIAL.to_string())
        .arg("--")
        .args(config.input_files)
        .stdin(Stdio::null())
        .stdout(Stdio::inherit())
        .stderr(Stdio::piped());

    log::info!("Counting k-mers on disk with {}", binary.display());
    log::debug!("{:?}", command);
    let mut child = command
        .spawn()
        .map_err(|e| format!("could not start {}: {e}", binary.display()))?;

    let stderr = child.stderr.take().expect("stderr is piped");
    let mut last_lines = VecDeque::with_capacity(ERROR_CONTEXT_LINES);
    for line in BufReader::new(stderr).lines() {
        let line = line.map_err(|e| format!("reading {BINARY_NAME} output: {e}"))?;
        log::debug!("[{BINARY_NAME}] {line}");
        if last_lines.len() == ERROR_CONTEXT_LINES {
            last_lines.pop_front();
        }
        last_lines.push_back(line);
    }

    let status = child
        .wait()
        .map_err(|e| format!("waiting for {BINARY_NAME}: {e}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!(
            "{BINARY_NAME} failed ({status}). Its last output was:\n{}",
            Vec::from(last_lines).join("\n")
        ))
    }
}

/// Finds `myloasm-kmc-v1` next to the current executable, then on `PATH`.
pub fn locate_binary() -> Result<PathBuf, String> {
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let sibling = dir.join(BINARY_NAME);
            if sibling.is_file() {
                return Ok(sibling);
            }
        }
    }
    if let Some(path) = std::env::var_os("PATH") {
        for dir in std::env::split_paths(&path) {
            let candidate = dir.join(BINARY_NAME);
            if candidate.is_file() {
                return Ok(candidate);
            }
        }
    }
    Err(format!(
        "could not find the `{BINARY_NAME}` executable next to myloasm or on PATH. \
         --kmc needs it: install it from https://github.com/bluenote-1577/myloasm-kmc \
         (`cargo install --path .` there) and make sure it is on PATH."
    ))
}
