/// Integration test: run the full pipeline through write_reads and validate
/// the produced R1/R2 fastq output against the golden files.
///
/// Golden parameters:
///   c=0.5, n=10000, u=300, sd=100, l=150, abundance=1.0
///   Default adapters (single-element pools → deterministic random.choice)
///   Basic mode (no quality scores), filler_base='G'
///
/// ORDERING NOTE: Python's process_df uses `df.sort_values(by=['length'])` with
/// the default `kind='quicksort'`, which is NOT stable for equal-length fragments.
/// We therefore compare the MULTISET of reads (sequence → count) rather than
/// requiring byte-exact line ordering. This is the correct semantic: two fastq
/// files are equivalent iff they contain the same reads in the same quantities.
use readsynth::digest::digest_fasta;
use readsynth::iupac::{check_for_palindrome, iupac_motifs, regex_bp_len};
use readsynth::prob_copies::{apply_probabilities, prob_to_counts};
use readsynth::size_select::size_select_gaussian;
use readsynth::types::MotifConfig;
use readsynth::write_reads::{write_reads, ReadConfig};
use std::collections::HashMap;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn build_hhai_msei_config() -> MotifConfig {
    let m1_raw = check_for_palindrome(&["GCG/C".to_string()]);
    let m2_raw = check_for_palindrome(&["T/TAA".to_string()]);
    let motif_dt1 = iupac_motifs(&m1_raw);
    let motif_dt2 = iupac_motifs(&m2_raw);
    let mut motif_dt = motif_dt1.clone();
    for (k, v) in &motif_dt2 {
        motif_dt.insert(k.clone(), *v);
    }
    let motif_len: HashMap<String, usize> =
        motif_dt.keys().map(|p| (p.clone(), regex_bp_len(p))).collect();
    MotifConfig { motif_dt, motif_dt1, motif_dt2, motif_len }
}

/// Run the full pipeline and return the R1 and R2 fastq content as byte vectors.
fn run_pipeline() -> (Vec<u8>, Vec<u8>) {
    let root = repo_root();
    let fasta = root.join("test/test_project4/genomes/Escherichia_coli_plasmid.fasta");

    let max_len = 900usize;
    let n       = 10_000usize;
    let u       = 300usize;
    let sd      = 100usize;
    let x       = u;
    let l1      = 150usize;

    let config = build_hhai_msei_config();
    let fragments = digest_fasta(&fasta, &config, max_len).expect("digest_fasta failed");
    let (count_frags, len_freqs) = apply_probabilities(fragments, 0.5, 1.0, max_len);
    let ss = size_select_gaussian(&len_freqs, n, u, sd, x, max_len).expect("size_select failed");
    let final_frags = prob_to_counts(count_frags, &ss.fragment_comps, ss.adjustment);

    let rc = ReadConfig::default_adapters(l1, None, 'G');
    let mut r1_buf: Vec<u8> = Vec::new();
    let mut r2_buf: Vec<u8> = Vec::new();
    write_reads(&final_frags, &mut r1_buf, &mut r2_buf, "Escherichia_coli_plasmid.fasta", &rc)
        .expect("write_reads failed");

    (r1_buf, r2_buf)
}

/// Parse a fastq byte buffer into a map of (seq_line → count).
/// Only the sequence lines (every 4th line starting at index 1) are counted.
fn seq_counts(fastq_bytes: &[u8]) -> HashMap<String, usize> {
    let text = std::str::from_utf8(fastq_bytes).expect("non-UTF8 fastq");
    let mut map: HashMap<String, usize> = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        if i % 4 == 1 {
            *map.entry(line.to_string()).or_insert(0) += 1;
        }
    }
    map
}

/// Count total lines in a fastq byte buffer.
fn line_count(fastq_bytes: &[u8]) -> usize {
    std::str::from_utf8(fastq_bytes)
        .expect("non-UTF8 fastq")
        .lines()
        .count()
}

#[test]
fn write_reads_r1_matches_golden() {
    let root = repo_root();
    let golden_r1 = root.join("test/test_project4/golden_output/sim_metagenome_R1.fastq");
    assert!(golden_r1.exists(), "Golden R1 not found: {}", golden_r1.display());

    let (r1_buf, _) = run_pipeline();
    let golden_bytes = std::fs::read(&golden_r1).expect("cannot read golden R1");

    // 1. Total read count must match (line count / 4).
    let rust_lines   = line_count(&r1_buf);
    let golden_lines = line_count(&golden_bytes);
    assert_eq!(
        rust_lines, golden_lines,
        "R1 line count mismatch: rust={rust_lines}, golden={golden_lines}"
    );

    // 2. Multiset of read sequences must match exactly.
    let rust_seqs   = seq_counts(&r1_buf);
    let golden_seqs = seq_counts(&golden_bytes);

    // Check every rust sequence appears in golden with the same count.
    let mut mismatches = 0usize;
    for (seq, &rc) in &rust_seqs {
        let gc = *golden_seqs.get(seq).unwrap_or(&0);
        if rc != gc {
            eprintln!("R1 count mismatch: seq={}...  rust={rc}, golden={gc}", &seq[..seq.len().min(40)]);
            mismatches += 1;
        }
    }
    // Check for sequences in golden not present in rust.
    for (seq, &gc) in &golden_seqs {
        if !rust_seqs.contains_key(seq) {
            eprintln!("R1 missing in rust: seq={}...  golden_count={gc}", &seq[..seq.len().min(40)]);
            mismatches += 1;
        }
    }

    assert_eq!(mismatches, 0, "{mismatches} R1 sequence(s) have wrong counts (see stderr)");

    println!(
        "write_reads R1: {} reads, {} unique sequences — multiset matches golden",
        rust_lines / 4,
        rust_seqs.len()
    );
}

#[test]
fn write_reads_r2_matches_golden() {
    let root = repo_root();
    let golden_r2 = root.join("test/test_project4/golden_output/sim_metagenome_R2.fastq");

    if !golden_r2.exists() {
        eprintln!("Skipping R2 test: golden R2 not found");
        return;
    }

    let (_, r2_buf) = run_pipeline();
    let golden_bytes = std::fs::read(&golden_r2).expect("cannot read golden R2");

    let rust_lines   = line_count(&r2_buf);
    let golden_lines = line_count(&golden_bytes);
    assert_eq!(rust_lines, golden_lines, "R2 line count mismatch");

    let rust_seqs   = seq_counts(&r2_buf);
    let golden_seqs = seq_counts(&golden_bytes);

    let mut mismatches = 0usize;
    for (seq, &rc) in &rust_seqs {
        let gc = *golden_seqs.get(seq).unwrap_or(&0);
        if rc != gc {
            eprintln!("R2 count mismatch: seq={}...  rust={rc}, golden={gc}", &seq[..seq.len().min(40)]);
            mismatches += 1;
        }
    }
    for (seq, &gc) in &golden_seqs {
        if !rust_seqs.contains_key(seq) {
            eprintln!("R2 missing in rust: seq={}...  golden_count={gc}", &seq[..seq.len().min(40)]);
            mismatches += 1;
        }
    }

    assert_eq!(mismatches, 0, "{mismatches} R2 sequence(s) have wrong counts");

    println!(
        "write_reads R2: {} reads, {} unique sequences — multiset matches golden",
        rust_lines / 4,
        rust_seqs.len()
    );
}

/// Verify the per-fragment counts (how many reads each fragment gets) match the
/// golden individual_counts CSV, independent of ordering.
#[test]
fn write_reads_per_fragment_counts_match_golden() {
    let root = repo_root();
    let golden_counts_csv = root.join(
        "test/test_project4/golden_output/individual_counts/counts_Escherichia_coli_plasmid.fasta.csv",
    );
    assert!(golden_counts_csv.exists(), "Golden counts CSV not found");

    let max_len = 900usize;
    let n       = 10_000usize;
    let u       = 300usize;
    let sd      = 100usize;
    let x       = u;
    let config  = build_hhai_msei_config();

    let fragments = digest_fasta(
        &repo_root().join("test/test_project4/genomes/Escherichia_coli_plasmid.fasta"),
        &config,
        max_len,
    ).expect("digest_fasta failed");
    let (count_frags, len_freqs) = apply_probabilities(fragments, 0.5, 1.0, max_len);
    let ss = size_select_gaussian(&len_freqs, n, u, sd, x, max_len).expect("size_select failed");
    let final_frags = prob_to_counts(count_frags, &ss.fragment_comps, ss.adjustment);

    // Index rust output by (start, end, reverse) → counts
    let rust_counts: HashMap<(usize, usize, u8), u64> = final_frags
        .iter()
        .map(|f| ((f.fragment.start, f.fragment.end, f.fragment.reverse), f.counts))
        .collect();

    // Load golden counts CSV
    let mut rdr = csv::Reader::from_path(&golden_counts_csv).expect("cannot open golden counts CSV");
    let headers = rdr.headers().expect("no headers").clone();
    let col = |name: &str| headers.iter().position(|h| h == name).unwrap();
    let i_start   = col("start");
    let i_end     = col("end");
    let i_rev     = col("reverse");
    let i_counts  = col("counts");

    let mut mismatches = 0;
    let mut total_golden = 0u64;
    for rec in rdr.records().map(|r| r.unwrap()) {
        let start:  usize = rec[i_start].parse().unwrap();
        let end:    usize = rec[i_end].parse().unwrap();
        let rev:    u8    = rec[i_rev].parse().unwrap();
        let golden_count: u64 = rec[i_counts].parse::<f64>().unwrap() as u64;
        total_golden += golden_count;

        let rust_count = *rust_counts.get(&(start, end, rev)).unwrap_or(&0);
        if rust_count != golden_count {
            eprintln!(
                "Count mismatch (start={start}, end={end}, rev={rev}): rust={rust_count}, golden={golden_count}"
            );
            mismatches += 1;
        }
    }

    assert_eq!(mismatches, 0, "{mismatches} fragment(s) have wrong read counts");
    println!(
        "per-fragment counts: {} fragments, total_golden={total_golden} — all match",
        final_frags.len()
    );
}
