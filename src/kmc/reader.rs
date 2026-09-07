//! Sequential reader for KMC databases, including the per-strand counter layout written by the
//! vendored KMC fork (`kmc -sc`, format versions 0x100 / 0x300). Pure Rust, no dependencies.
//! Kept in sync with the reader in the KMC fork repository (kmc_rust_disk_reader/src/lib.rs).

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// KMC database file information
#[derive(Debug, Clone)]
pub struct KmcFileInfo {
    pub kmer_length: u32,
    pub mode: u32,
    pub counter_size: u32,
    pub lut_prefix_length: u32,
    pub signature_len: u32,
    pub min_count: u32,
    pub max_count: u64,
    pub total_kmers: u64,
    pub both_strands: bool,
    pub kmc_version: u32,
    /// Number of counters stored per k-mer: 2 for a database built with stranded
    /// counters (`kmc -sc`), 1 otherwise.
    pub n_counters: u32,
    /// Minimal per-strand count the database was built with (stranded counters only).
    pub min_count_per_strand: u32,
}

impl KmcFileInfo {
    /// True if the database stores a separate counter for each strand.
    #[inline]
    pub fn is_stranded(&self) -> bool {
        self.n_counters == 2
    }
}

/// A kmer with its count (uses borrowed slice for zero-copy)
#[derive(Debug)]
pub struct KmerRecord<'a> {
    pub kmer: &'a [u8],
    /// Total number of occurrences (`count_fwd + count_rev`).
    pub count: u64,
    /// Times the k-mer was read as itself. Equals `count` for a database without
    /// stranded counters.
    pub count_fwd: u64,
    /// Times the k-mer was read as its reverse complement. Always 0 for a database
    /// without stranded counters.
    pub count_rev: u64,
}

/// KMC database reader for sequential listing
pub struct KmcReader {
    info: KmcFileInfo,

    // Prefix file state
    prefix_buf: Vec<u64>,
    prefix_pos: usize,
    prefix_mask: u64,

    // Suffix file state - single buffer, no BufReader overhead
    suffix_file: File,
    suffix_buf: Vec<u8>,
    suffix_buf_pos: usize,
    suffix_buf_len: usize,

    // Iteration state
    suffix_number: u64,

    // Derived values
    suffix_size: u32,
    suffix_rec_size: u32,
    byte_alignment: u8,
    no_of_rows: u32,

    // Pre-allocated buffers for hot path
    kmer_data: Vec<u64>,
    kmer_str_buf: Vec<u8>,
}

const SUFFIX_BUF_SIZE: usize = 1 << 25; // 32MB buffer (matches C++)
const CHAR_CODES: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl KmcReader {
    /// Open a KMC database for sequential reading
    pub fn open<P: AsRef<Path>>(db_path: P) -> io::Result<Self> {
        let db_path = db_path.as_ref();
        let db_str = db_path.to_string_lossy();

        // Remove .kmc_pre or .kmc_suf extension if present
        let base_path = db_str
            .strip_suffix(".kmc_pre")
            .or_else(|| db_str.strip_suffix(".kmc_suf"))
            .unwrap_or(&db_str);

        let pre_path = format!("{}.kmc_pre", base_path);
        let suf_path = format!("{}.kmc_suf", base_path);

        // Open and read prefix file
        let mut pre_file = File::open(&pre_path)?;
        let pre_size = pre_file.metadata()?.len();

        // Verify markers
        let mut marker = [0u8; 4];
        pre_file.read_exact(&mut marker)?;
        if &marker != b"KMCP" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid KMCP header",
            ));
        }

        pre_file.seek(SeekFrom::End(-4))?;
        pre_file.read_exact(&mut marker)?;
        if &marker != b"KMCP" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid KMCP trailer",
            ));
        }

        // Read KMC version
        pre_file.seek(SeekFrom::End(-12))?;
        let mut buf4 = [0u8; 4];
        pre_file.read_exact(&mut buf4)?;
        let kmc_version = u32::from_le_bytes(buf4);

        // 0x100 / 0x300 are the KMC1 / KMC2 layouts whose suffix records carry a pair
        // of per-strand counters (forward first, then reverse).
        if kmc_version != 0 && kmc_version != 0x200 && kmc_version != 0x100 && kmc_version != 0x300
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Unsupported KMC version: 0x{:x}", kmc_version),
            ));
        }
        let is_kmc2 = kmc_version == 0x200 || kmc_version == 0x300;

        // Read header offset
        pre_file.seek(SeekFrom::End(-8))?;
        let mut header_offset_byte = [0u8; 1];
        pre_file.read_exact(&mut header_offset_byte)?;
        let header_offset = header_offset_byte[0] as i64;

        // Read metadata
        pre_file.seek(SeekFrom::End(-(header_offset + 8)))?;

        let info = if is_kmc2 {
            Self::read_kmc2_header(&mut pre_file, kmc_version)?
        } else {
            Self::read_kmc1_header(&mut pre_file, kmc_version)?
        };

        // Calculate derived values
        let byte_alignment = if info.kmer_length % 4 != 0 {
            4 - (info.kmer_length % 4) as u8
        } else {
            0
        };

        let suffix_size = (info.kmer_length - info.lut_prefix_length) / 4;
        let suffix_rec_size = suffix_size + info.counter_size * info.n_counters;

        let no_of_rows = {
            let total_symbols = info.kmer_length + byte_alignment as u32;
            if total_symbols % 32 != 0 {
                total_symbols / 32 + 1
            } else {
                total_symbols / 32
            }
        };

        // Read prefix LUT
        let prefix_mask = (1u64 << (2 * info.lut_prefix_length)) - 1;

        let (prefix_buf, _signature_map) = if is_kmc2 {
            Self::read_kmc2_lut(&mut pre_file, pre_size, &info, header_offset)?
        } else {
            Self::read_kmc1_lut(&mut pre_file, &info)?
        };

        // Open suffix file
        let mut suf_file = File::open(&suf_path)?;

        // Verify suffix file markers
        suf_file.read_exact(&mut marker)?;
        if &marker != b"KMCS" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid KMCS header",
            ));
        }

        suf_file.seek(SeekFrom::End(-4))?;
        suf_file.read_exact(&mut marker)?;
        if &marker != b"KMCS" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid KMCS trailer",
            ));
        }

        // Position at start of suffix data
        suf_file.seek(SeekFrom::Start(4))?;

        // Pre-allocate buffers
        let kmer_data = vec![0u64; no_of_rows as usize];
        let kmer_str_buf = vec![0u8; info.kmer_length as usize];

        Ok(KmcReader {
            info,
            prefix_buf,
            prefix_pos: 1,
            prefix_mask,
            suffix_file: suf_file,
            suffix_buf: vec![0u8; SUFFIX_BUF_SIZE],
            suffix_buf_pos: 0,
            suffix_buf_len: 0,
            suffix_number: 0,
            suffix_size,
            suffix_rec_size,
            byte_alignment,
            no_of_rows,
            kmer_data,
            kmer_str_buf,
        })
    }

    fn read_kmc1_header(file: &mut File, kmc_version: u32) -> io::Result<KmcFileInfo> {
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        file.read_exact(&mut buf4)?;
        let kmer_length = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let mode = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let counter_size = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let lut_prefix_length = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let min_count = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let max_count_lo = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf8)?;
        let total_kmers = u64::from_le_bytes(buf8);

        let mut both_strands_byte = [0u8; 1];
        file.read_exact(&mut both_strands_byte)?;
        let both_strands = both_strands_byte[0] == 0;

        // counter mode: 1 means a pair of per-strand counters, then 2 reserved bytes
        let mut counter_mode_byte = [0u8; 1];
        file.read_exact(&mut counter_mode_byte)?;
        let n_counters = if kmc_version == 0x100 {
            if counter_mode_byte[0] != 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "database claims per-strand counters but the counter mode says otherwise",
                ));
            }
            2
        } else {
            1
        };
        file.seek(SeekFrom::Current(2))?;

        file.read_exact(&mut buf4)?;
        let max_count_hi = u32::from_le_bytes(buf4);
        let max_count = (max_count_hi as u64) << 32 | max_count_lo as u64;

        let min_count_per_strand = if kmc_version == 0x100 {
            file.read_exact(&mut buf4)?;
            u32::from_le_bytes(buf4)
        } else {
            0
        };

        Ok(KmcFileInfo {
            kmer_length,
            mode,
            counter_size,
            lut_prefix_length,
            signature_len: 0,
            min_count,
            max_count,
            total_kmers,
            both_strands,
            kmc_version,
            n_counters,
            min_count_per_strand,
        })
    }

    fn read_kmc2_header(file: &mut File, kmc_version: u32) -> io::Result<KmcFileInfo> {
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        file.read_exact(&mut buf4)?;
        let kmer_length = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let mode = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let counter_size = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let lut_prefix_length = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let signature_len = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let min_count = u32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let max_count = u32::from_le_bytes(buf4) as u64;

        file.read_exact(&mut buf8)?;
        let total_kmers = u64::from_le_bytes(buf8);

        let mut both_strands_byte = [0u8; 1];
        file.read_exact(&mut both_strands_byte)?;
        let both_strands = both_strands_byte[0] == 0;

        let (n_counters, min_count_per_strand) = if kmc_version == 0x300 {
            // counter mode: 1 means a pair of per-strand counters, then the per-strand cutoff
            let mut counter_mode_byte = [0u8; 1];
            file.read_exact(&mut counter_mode_byte)?;
            if counter_mode_byte[0] != 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "database claims per-strand counters but the counter mode says otherwise",
                ));
            }
            file.read_exact(&mut buf4)?;
            (2, u32::from_le_bytes(buf4))
        } else {
            (1, 0)
        };

        Ok(KmcFileInfo {
            kmer_length,
            mode,
            counter_size,
            lut_prefix_length,
            signature_len,
            min_count,
            max_count,
            total_kmers,
            both_strands,
            kmc_version,
            n_counters,
            min_count_per_strand,
        })
    }

    fn read_kmc1_lut(
        file: &mut File,
        info: &KmcFileInfo,
    ) -> io::Result<(Vec<u64>, Option<Vec<u32>>)> {
        let lut_size = (1u64 << (2 * info.lut_prefix_length)) + 1;

        file.seek(SeekFrom::Start(4))?;

        // Read entire LUT at once as bytes, then convert
        let byte_size = lut_size as usize * 8;
        let mut raw_buf = vec![0u8; byte_size];
        file.read_exact(&mut raw_buf)?;

        // Convert bytes to u64s (little-endian)
        let mut prefix_buf = vec![0u64; lut_size as usize];
        for i in 0..lut_size as usize {
            let offset = i * 8;
            prefix_buf[i] = u64::from_le_bytes(raw_buf[offset..offset + 8].try_into().unwrap());
        }

        prefix_buf[lut_size as usize - 1] = info.total_kmers;

        Ok((prefix_buf, None))
    }

    fn read_kmc2_lut(
        file: &mut File,
        file_size: u64,
        info: &KmcFileInfo,
        header_offset: i64,
    ) -> io::Result<(Vec<u64>, Option<Vec<u32>>)> {
        let signature_map_size = (1u32 << (2 * info.signature_len)) + 1;

        let content_size = file_size - 8;
        let lut_area_size =
            content_size as i64 - 4 - (signature_map_size as i64 * 4) - header_offset;

        let last_data_index = lut_area_size as u64 / 8;

        file.seek(SeekFrom::Start(4))?;

        // Read entire LUT at once as bytes
        let byte_size = last_data_index as usize * 8;
        let mut raw_buf = vec![0u8; byte_size];
        file.read_exact(&mut raw_buf)?;

        // Convert bytes to u64s (little-endian)
        let mut prefix_buf = vec![0u64; (last_data_index + 1) as usize];
        for i in 0..last_data_index as usize {
            let offset = i * 8;
            prefix_buf[i] = u64::from_le_bytes(raw_buf[offset..offset + 8].try_into().unwrap());
        }

        prefix_buf[last_data_index as usize] = info.total_kmers;

        Ok((prefix_buf, None))
    }

    /// Get database information
    #[inline]
    pub fn info(&self) -> &KmcFileInfo {
        &self.info
    }

    /// Get the current prefix for a suffix number
    #[inline]
    fn get_prefix(&mut self) -> u64 {
        while self.prefix_pos < self.prefix_buf.len() {
            if self.suffix_number != self.prefix_buf[self.prefix_pos] {
                break;
            }
            self.prefix_pos += 1;
        }
        ((self.prefix_pos - 1) as u64) & self.prefix_mask
    }

    /// Ensure we have enough bytes in the buffer, refill if needed
    #[inline]
    fn ensure_buffer(&mut self, need: usize) -> io::Result<()> {
        if self.suffix_buf_pos + need > self.suffix_buf_len {
            // Move remaining bytes to start
            let remaining = self.suffix_buf_len - self.suffix_buf_pos;
            if remaining > 0 {
                self.suffix_buf
                    .copy_within(self.suffix_buf_pos..self.suffix_buf_len, 0);
            }
            self.suffix_buf_pos = 0;
            self.suffix_buf_len = remaining;

            // Read more data - loop to fill buffer (read() may return less than requested)
            while self.suffix_buf_len < self.suffix_buf.len() {
                let read = self
                    .suffix_file
                    .read(&mut self.suffix_buf[self.suffix_buf_len..])?;
                if read == 0 {
                    break; // EOF
                }
                self.suffix_buf_len += read;
            }
        }
        Ok(())
    }

    /// Read the next kmer into internal buffer, return borrowed reference
    #[inline]
    pub fn next_kmer(&mut self) -> io::Result<Option<KmerRecord<'_>>> {
        if self.suffix_number >= self.info.total_kmers {
            return Ok(None);
        }

        // Get prefix
        let prefix = self.get_prefix();

        // Reset kmer_data
        for i in 0..self.no_of_rows as usize {
            self.kmer_data[i] = 0;
        }

        // Place prefix
        let off = 64 - (self.info.lut_prefix_length * 2) - (self.byte_alignment as u32 * 2);
        self.kmer_data[0] = prefix << off;

        // Ensure we have enough data for suffix + counter
        self.ensure_buffer(self.suffix_rec_size as usize)?;

        // Read suffix bytes directly from buffer
        let mut row_index = 0usize;
        let mut bit_off = off as i32 - 8;

        for i in 0..self.suffix_size as usize {
            let suf_byte = self.suffix_buf[self.suffix_buf_pos + i] as u64;
            self.kmer_data[row_index] |= suf_byte << bit_off;

            if bit_off == 0 {
                bit_off = 56;
                row_index += 1;
            } else {
                bit_off -= 8;
            }
        }
        self.suffix_buf_pos += self.suffix_size as usize;

        // Read counter(s): a stranded database stores the forward one first, then the reverse one
        let (count_fwd, count_rev) = if self.info.counter_size == 0 {
            (1u64, 0u64)
        } else {
            let counter_size = self.info.counter_size as usize;
            let mut counts = [0u64; 2];
            for (c, slot) in counts
                .iter_mut()
                .enumerate()
                .take(self.info.n_counters as usize)
            {
                let base = self.suffix_buf_pos + c * counter_size;
                let mut value = self.suffix_buf[base] as u64;
                for b in 1..counter_size {
                    value |= (self.suffix_buf[base + b] as u64) << (8 * b);
                }
                *slot = value;
            }
            self.suffix_buf_pos += counter_size * self.info.n_counters as usize;
            (counts[0], counts[1])
        };

        // Convert to string in pre-allocated buffer
        self.kmer_data_to_buf();

        self.suffix_number += 1;

        Ok(Some(KmerRecord {
            kmer: &self.kmer_str_buf,
            count: count_fwd + count_rev,
            count_fwd,
            count_rev,
        }))
    }

    /// Convert kmer_data to string in pre-allocated buffer
    #[inline]
    fn kmer_data_to_buf(&mut self) {
        let mut temp_byte_alignment = self.byte_alignment;
        let mut pos = 0usize;
        let kmer_len = self.info.kmer_length as usize;

        for row_counter in 0..self.no_of_rows as usize {
            let word = self.kmer_data[row_counter];
            let bytes = word.to_be_bytes();

            let mut i = 0u32;
            let mut byte_idx = 0usize;

            while i < 32 && byte_idx < 8 && pos < kmer_len {
                let byte = bytes[byte_idx];

                // Symbol 1 (top 2 bits)
                if i == 0 && temp_byte_alignment > 0 {
                    temp_byte_alignment -= 1;
                } else if pos < kmer_len {
                    self.kmer_str_buf[pos] = CHAR_CODES[((byte >> 6) & 0x03) as usize];
                    pos += 1;
                }

                // Symbol 2
                if i == 0 && temp_byte_alignment > 0 {
                    temp_byte_alignment -= 1;
                } else if pos < kmer_len {
                    self.kmer_str_buf[pos] = CHAR_CODES[((byte >> 4) & 0x03) as usize];
                    pos += 1;
                }

                // Symbol 3
                if i == 0 && temp_byte_alignment > 0 {
                    temp_byte_alignment -= 1;
                } else if pos < kmer_len {
                    self.kmer_str_buf[pos] = CHAR_CODES[((byte >> 2) & 0x03) as usize];
                    pos += 1;
                }

                // Symbol 4
                if pos < kmer_len {
                    self.kmer_str_buf[pos] = CHAR_CODES[(byte & 0x03) as usize];
                    pos += 1;
                }

                i += 4;
                byte_idx += 1;
            }
        }
    }

    /// Write all kmers directly to a writer (most efficient for dumping).
    /// A stranded database gets two counter columns: forward, then reverse.
    pub fn write_all<W: Write>(&mut self, mut writer: W) -> io::Result<u64> {
        let mut count = 0u64;
        let mut line_buf = Vec::with_capacity(self.info.kmer_length as usize + 48);
        let stranded = self.info.is_stranded();

        while let Some(record) = self.next_kmer()? {
            line_buf.clear();
            line_buf.extend_from_slice(record.kmer);
            line_buf.push(b'\t');
            if stranded {
                push_u64(&mut line_buf, record.count_fwd);
                line_buf.push(b'\t');
                push_u64(&mut line_buf, record.count_rev);
            } else {
                push_u64(&mut line_buf, record.count);
            }
            line_buf.push(b'\n');

            writer.write_all(&line_buf)?;
            count += 1;
        }
        Ok(count)
    }
}

/// Append the decimal representation of `num` to `out`.
#[inline]
pub fn push_u64(out: &mut Vec<u8>, mut num: u64) {
    let start = out.len();
    if num == 0 {
        out.push(b'0');
        return;
    }
    while num > 0 {
        out.push(b'0' + (num % 10) as u8);
        num /= 10;
    }
    out[start..].reverse();
}

/// Convenience function to open a database
pub fn read_kmers<P: AsRef<Path>>(db_path: P) -> io::Result<KmcReader> {
    KmcReader::open(db_path)
}
