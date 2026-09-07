use crate::cli::Cli;
use crate::constants::{MIN_SPLIT_KMER_COUNT_PER_STRAND, MIN_SPLIT_KMER_COUNT_TOTAL};
use crate::seeding;
use crate::utils::*;
use crossbeam_channel::bounded;
use fastbloom::BloomFilter;
use fxhash::FxBuildHasher;
use fxhash::FxHashMap;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;

#[inline]
pub fn quickhash(kmer: u64) -> u64 {
    //fxhash
    //return fxhash::hash64(&kmer);
    return seeding::mm_hash64(kmer);
}

/// Whether a split k-mer with `[reverse count, forward count]` is kept for SNPmer detection.
#[inline]
pub fn split_kmer_passes_count_thresholds(counts: [u32; 2]) -> bool {
    counts[0] >= MIN_SPLIT_KMER_COUNT_PER_STRAND
        && counts[1] >= MIN_SPLIT_KMER_COUNT_PER_STRAND
        && counts[0].saturating_add(counts[1]) >= MIN_SPLIT_KMER_COUNT_TOTAL
}

/// Counts all split k-mers of the input reads in memory and returns those passing the count
/// thresholds as `(canonical k-mer, [reverse count, forward count])`. `--kmc` does
/// the same through `kmc::split_kmers_via_disk_count`.
pub fn read_to_split_kmers(k: usize, threads: usize, args: &Cli) -> Vec<(u64, [u32; 2])> {
    let start = std::time::Instant::now();
    let bf_vec_maps = first_iteration(k, threads, args);
    log::info!(
        "Finished with bloom filter processing in {:?}. Round 2 - Start k-mer counting...",
        start.elapsed()
    );
    log_memory_usage(true, "Memory usage after bloom filter processing");

    let start = std::time::Instant::now();
    let vec_maps = second_iteration(k, threads, args, bf_vec_maps);
    let map_size_raw = vec_maps.iter().map(|x| x.len()).sum::<usize>();
    log::info!(
        "Finished with k-mer counting in {:?}. Total kmers after bloom filter: {}",
        start.elapsed(),
        map_size_raw
    );
    log_memory_usage(
        true,
        "Memory usage after second round of k-mer counting processing",
    );

    let mut vectorized_map = vec![];
    for map in vec_maps.into_iter() {
        for (kmer, counts) in map.into_iter() {
            if split_kmer_passes_count_thresholds(counts) {
                vectorized_map.push((kmer, counts));
            }
        }
    }

    let map_size_retain = vectorized_map.len();
    log::info!(
        "Removed {} kmers with counts < 1 in both strands and <= 3 multiplicity.",
        map_size_raw - map_size_retain
    );
    if map_size_retain < map_size_raw / 1000 {
        log::warn!("Less than 0.1% of kmers have counts > 1 in both strands and > 2 multiplicity. This may indicate a problem with the input data or very low coverage.");
    }
    log::debug!("Final Hashmap len after vectorization: {}", map_size_retain);
    log_memory_usage(
        false,
        "Memory usage after second round of k-mer counting processing",
    );

    return vectorized_map;
}

fn estimate_bf_size(args: &Cli) -> f64 {
    let mut est_bases = 0.;
    let mut is_gzipped = false;
    for fq_file in args.input_files.iter() {
        if fq_file.contains(".gz") || fq_file.contains(".gzip") || fq_file.contains(".bz") {
            is_gzipped = true;
        }

        // Fail if file can not be read
        if Path::new(fq_file).metadata().is_err() {
            log::error!(
                "Unable to read file: {}. Please check the file path and permissions.",
                fq_file
            );
            std::process::exit(1);
        }

        let metadata = std::fs::metadata(fq_file).expect("Unable to read file metadata");
        log::debug!(
            "File: {}, size (Gbytes): {}",
            fq_file,
            metadata.len() as f64 / 1_000_000_000.
        );
        est_bases += metadata.len() as f64 / 1_000_000_000.;
        if is_gzipped {
            est_bases *= 1.5; //rough estimate of compression ratio
        } else {
            est_bases /= 2.;
        }
    }
    let bf_size = (est_bases / 2.0).min(200.0).max(2.0);
    return bf_size;
}

fn first_iteration(k: usize, threads: usize, args: &Cli) -> Vec<FxHashMap<u64, [u32; 2]>> {
    //Topology is
    //      A-SEND: tx_head, , B-REC: rx_head1, rx_head2...
    // |   |  ...
    // B   B  ...  B-SEND: txs[0...], txs2[0...],... C-REC: rxs
    // | x | x | ...
    // C   C  ...
    let hm_size = threads;
    let mask = !(1 << 63);
    let bf_size;
    if let Some(bf_size_manual) = args.bloom_filter_size {
        bf_size = bf_size_manual;
    } else {
        bf_size = estimate_bf_size(args);
        log::info!("Using automatic bloom filter size: {:.2} Gbytes", bf_size);
    }

    let aggressive_bloom = args.aggressive_bloom;
    let mut bf_vec_maps: Vec<FxHashMap<u64, [u32; 2]>> = vec![FxHashMap::default(); hm_size];
    if bf_size > 0. {
        let num_b = threads / 10 + 1;
        let counter = Arc::new(Mutex::new(0));
        let read_lengths = Arc::new(Mutex::new(vec![]));
        let mut rxs = vec![];
        let mut txs_vecs = vec![vec![]; num_b];
        for _ in 0..threads {
            //let (tx, rx) = unbounded();
            let (tx, rx) = bounded(500);
            for i in 1..num_b {
                txs_vecs[i].push(tx.clone());
            }
            txs_vecs[0].push(tx);
            rxs.push(rx);
        }

        //let (tx_head, rx_head1) = unbounded();
        let (tx_head, rx_head1) = bounded(500);
        let mut rx_heads = vec![];
        for _ in 1..num_b {
            let rx_head2 = rx_head1.clone();
            rx_heads.push(rx_head2);
        }
        rx_heads.push(rx_head1);

        assert!(txs_vecs.len() == rx_heads.len());

        let fq_files = args.input_files.clone();
        //A: Get k-mers
        thread::spawn(move || {
            for fq_file in fq_files {
                let bufreader = BufReader::new(std::fs::File::open(fq_file).expect("valid path"));
                let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
                while let Some(record) = reader.next() {
                    let rec = record.expect("Error reading record");
                    let seq = rec.seq().to_vec();
                    let qualities = rec.qual().map(Vec::from);
                    tx_head.send((seq, qualities)).unwrap();
                }
            }
            drop(tx_head);
            log::debug!("Finished reading all reads.");
        });
        //B: Process kmers and send to hash maps
        for (rx_head, txs) in rx_heads.into_iter().zip(txs_vecs.into_iter()) {
            let clone_counter = Arc::clone(&counter);
            let clone_read_lengths = Arc::clone(&read_lengths);
            thread::spawn(move || {
                let mut split_kmer_buf: Vec<u64> = Vec::new();
                let mut vec_and_canon: Vec<Vec<u64>> = vec![Vec::new(); hm_size];
                loop {
                    match rx_head.recv() {
                        Ok((seq, qualities)) => {
                            {
                                let mut read_lengths = clone_read_lengths.lock().unwrap();
                                read_lengths.push(seq.len());
                            }
                            seeding::split_kmer_mid(
                                &seq,
                                qualities.as_deref(),
                                k,
                                &mut split_kmer_buf,
                            );
                            for &kmer_i_and_canon in split_kmer_buf.iter() {
                                let kmer = kmer_i_and_canon & mask;
                                let hash = quickhash(kmer) % hm_size as u64;
                                vec_and_canon[hash as usize].push(kmer_i_and_canon);
                            }

                            for (i, bucket) in vec_and_canon.iter_mut().enumerate() {
                                if !bucket.is_empty() {
                                    let cap = bucket.capacity();
                                    let to_send =
                                        std::mem::replace(bucket, Vec::with_capacity(cap));
                                    txs[i].send(to_send).unwrap();
                                }
                            }

                            {
                                let mut counter = clone_counter.lock().unwrap();
                                *counter += 1;
                                if *counter % 100000 == 0 {
                                    log::info!("Processed {} reads.", counter);
                                }
                                if *counter % 1_000_000 == 0 {
                                    log_memory_usage(
                                        false,
                                        &format!(
                                            "Processed {} reads for bloom filter stage",
                                            *counter
                                        ),
                                    );
                                }
                            }
                        }
                        Err(_) => {
                            break;
                        }
                    }
                }
                for tx in txs {
                    drop(tx);
                }
            });
        }

        //C: Update bloom filter
        let mut handles = Vec::new();
        for rx in rxs.into_iter() {
            handles.push(thread::spawn(move || {
                let mut filter_canonical = BloomFilter::with_num_bits(
                    (bf_size * 4. * 1_000_000_000. / threads as f64) as usize,
                )
                .hasher(FxBuildHasher::default())
                .expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut filter_noncanonical = BloomFilter::with_num_bits(
                    (bf_size * 4. * 1_000_000_000. / threads as f64) as usize,
                )
                .hasher(FxBuildHasher::default())
                .expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut map: FxHashMap<u64, [u32; 2]> = FxHashMap::default();
                loop {
                    match rx.recv() {
                        Ok(msg) => {
                            let kmer_vecs = msg;
                            for kmer_i_canon in kmer_vecs {
                                let canonical = kmer_i_canon >> 63;
                                let kmer = kmer_i_canon & mask;
                                let kmer_canon = kmer | (1 << 63);
                                if canonical == 1 {
                                    let already_present_canon =
                                        filter_canonical.insert(&kmer_canon);
                                    let already_present_noncanon =
                                        filter_noncanonical.contains(&kmer);
                                    if aggressive_bloom {
                                        if already_present_noncanon && already_present_canon {
                                            map.insert(kmer, [0, 0]);
                                        }
                                    } else {
                                        if already_present_noncanon {
                                            map.insert(kmer, [0, 0]);
                                        }
                                    }
                                } else {
                                    let already_present_noncanon =
                                        filter_noncanonical.insert(&kmer);
                                    let already_present_canon =
                                        filter_canonical.contains(&kmer_canon);
                                    if aggressive_bloom {
                                        if already_present_noncanon && already_present_canon {
                                            map.insert(kmer, [0, 0]);
                                        }
                                    } else {
                                        if already_present_canon {
                                            map.insert(kmer, [0, 0]);
                                        }
                                    }
                                }
                            }
                        }
                        Err(_) => {
                            log::trace!("Thread finished.");
                            break;
                        }
                    }
                }

                map.shrink_to_fit();
                map
            }));
        }

        for (map_ind, handle) in handles.into_iter().enumerate() {
            bf_vec_maps[map_ind] = handle.join().unwrap();
            bf_vec_maps[map_ind].shrink_to_fit();
        }

        let read_lengths = Arc::try_unwrap(read_lengths).unwrap().into_inner().unwrap();
        log::info!(
            "Read lengths - {}",
            get_nx_from_vec(&read_lengths, &[10, 50, 90])
        );
        log::info!(
            "Total bases - {} million bp",
            (read_lengths.iter().map(|x| *x as usize).sum::<usize>() as f64 / 1_000_000.).round()
        );
    }

    return bf_vec_maps;
}

fn second_iteration(
    k: usize,
    threads: usize,
    args: &Cli,
    bf_vec_maps: Vec<FxHashMap<u64, [u32; 2]>>,
) -> Vec<FxHashMap<u64, [u32; 2]>> {
    let bf_size = args.bloom_filter_size;
    let mask = !(1 << 63);
    let mut vec_maps: Vec<FxHashMap<u64, [u32; 2]>> = vec![FxHashMap::default(); threads];

    let num_b = threads / 10 + 1;
    let counter = Arc::new(Mutex::new(0));
    let mut rxs = vec![];
    let mut txs_vecs = vec![vec![]; num_b];
    for _ in 0..threads {
        //let (tx, rx) = unbounded();
        let (tx, rx) = bounded(500);
        for i in 1..num_b {
            txs_vecs[i].push(tx.clone());
        }
        txs_vecs[0].push(tx);
        rxs.push(rx);
    }

    //let (tx_head, rx_head1) = unbounded();
    let (tx_head, rx_head1) = bounded(500);
    let mut rx_heads = vec![];
    for _ in 1..num_b {
        let rx_head2 = rx_head1.clone();
        rx_heads.push(rx_head2);
    }
    rx_heads.push(rx_head1);
    let bf_size = if let Some(bf_size_manual) = bf_size {
        bf_size_manual
    } else {
        estimate_bf_size(args)
    };

    assert!(txs_vecs.len() == rx_heads.len());

    let fq_files = args.input_files.clone();
    thread::spawn(move || {
        for fq_file in fq_files {
            let bufreader = BufReader::new(std::fs::File::open(fq_file).expect("valid path"));
            let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                let seq = rec.seq().to_vec();
                let qualities = rec.qual().map(Vec::from);
                tx_head.send((seq, qualities)).unwrap();
            }
        }
        drop(tx_head);
        log::debug!("Finished reading all reads.");
    });

    //B: Process kmers and send to hash maps
    for (rx_head, txs) in rx_heads.into_iter().zip(txs_vecs.into_iter()) {
        let clone_counter = Arc::clone(&counter);
        thread::spawn(move || {
            let mut split_kmer_buf: Vec<u64> = Vec::new();
            let mut vec_and_canon: Vec<Vec<u64>> = vec![Vec::new(); threads];
            loop {
                match rx_head.recv() {
                    Ok((seq, qualities)) => {
                        seeding::split_kmer_mid(&seq, qualities.as_deref(), k, &mut split_kmer_buf);
                        for &kmer_i_and_canon in split_kmer_buf.iter() {
                            let kmer = kmer_i_and_canon & mask;
                            let hash = quickhash(kmer) % threads as u64;
                            vec_and_canon[hash as usize].push(kmer_i_and_canon);
                        }

                        for (i, bucket) in vec_and_canon.iter_mut().enumerate() {
                            if !bucket.is_empty() {
                                let cap = bucket.capacity();
                                let to_send = std::mem::replace(bucket, Vec::with_capacity(cap));
                                txs[i].send(to_send).unwrap();
                            }
                        }

                        {
                            let mut counter = clone_counter.lock().unwrap();
                            *counter += 1;
                            if *counter % 100000 == 0 {
                                log::debug!("Processed {} reads.", counter);
                            }
                        }
                    }
                    Err(_) => {
                        break;
                    }
                }
            }
            for tx in txs {
                drop(tx);
            }
        });
    }

    let mut handles = Vec::new();
    for (rx, my_map) in rxs.into_iter().zip(bf_vec_maps.into_iter()) {
        handles.push(thread::spawn(move || {
            let mut my_map = my_map;
            loop {
                match rx.recv() {
                    Ok(msg) => {
                        let vec_and_canon = msg;
                        if bf_size > 0. {
                            for kmer_and_canon in vec_and_canon.into_iter() {
                                let kmer = kmer_and_canon & mask;
                                let canon = kmer_and_canon >> 63;
                                if let Some(val) = my_map.get_mut(&kmer) {
                                    val[canon as usize] += 1;
                                }
                            }
                        } else {
                            for kmer_and_canon in vec_and_canon.into_iter() {
                                let kmer = kmer_and_canon & mask;
                                let canon = kmer_and_canon >> 63;
                                let val = my_map.entry(kmer).or_insert([0, 0]);
                                val[canon as usize] += 1
                            }
                        }
                    }
                    Err(_) => {
                        log::trace!("Thread finished.");
                        break;
                    }
                }
            }
            my_map.shrink_to_fit();
            my_map
        }));
    }

    for (i, handle) in handles.into_iter().enumerate() {
        vec_maps[i] = handle.join().unwrap();
    }

    return vec_maps;
}

// Intuitively I think this helps, not sure if we want to use it TODO
pub fn quality_pool(qualities: Vec<u8>) -> Vec<u8> {
    let pool_width = 5;
    let mut pool = Vec::new();
    for i in 0..qualities.len() {
        if i > pool_width / 2 && i < qualities.len() - pool_width / 2 {
            let mut min = 255;
            for j in i - pool_width / 2..i + pool_width / 2 {
                if qualities[j] < min {
                    min = qualities[j];
                }
            }
            pool.push(min);
        } else {
            pool.push(qualities[i]);
        }
    }
    return pool;
}
