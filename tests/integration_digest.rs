/// Integration test: digest the E. coli plasmid with the same parameters used
/// to generate the golden output, then compare row-for-row against the golden CSV.
///
/// Golden parameters (from test/test_project4/make_golden_standard.sh):
///   -m1 GCG/C  -m2 T/TAA  -c 0.5  -n 10000  -u 300  -sd 100  -l 150
///   max_len = u + 6*sd = 300 + 600 = 900
///
/// The golden raw_digest CSV was produced by process_df() in readsynth.py and lives at:
///   test/test_project4/golden_output/raw_digests/raw_digest_Escherichia_coli_plasmid.fasta.csv
use readsynth::digest::digest_fasta;
use readsynth::iupac::{check_for_palindrome, iupac_motifs, regex_bp_len};
use readsynth::types::{MotifConfig, ProcessedFragment};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build the MotifConfig for GCG/C + T/TAA, mirroring readsynth.py's setup block.
/// Palindrome check is applied (both motifs are palindromic so no extras are added).
fn build_hhai_msei_config() -> MotifConfig {
    let m1_raw = check_for_palindrome(&["GCG/C".to_string()]);
    let m2_raw = check_for_palindrome(&["T/TAA".to_string()]);

    let motif_dt1 = iupac_motifs(&m1_raw);
    let motif_dt2 = iupac_motifs(&m2_raw);

    // motif_dt = dt1 first, then dt2 (insertion order matters for trimming)
    let mut motif_dt = motif_dt1.clone();
    for (k, v) in &motif_dt2 {
        motif_dt.insert(k.clone(), *v);
    }

    let motif_len: HashMap<String, usize> = motif_dt
        .keys()
        .map(|p| (p.clone(), regex_bp_len(p)))
        .collect();

    MotifConfig { motif_dt, motif_dt1, motif_dt2, motif_len }
}

/// Return the repo root (two levels up from readsynth-rs/).
fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..")
}

/// Load the golden raw_digest CSV into a Vec<ProcessedFragment>.
fn load_golden(path: &Path) -> Vec<ProcessedFragment> {
    let mut rdr = csv::Reader::from_path(path)
        .unwrap_or_else(|e| panic!("Cannot open golden CSV {}: {e}", path.display()));

    rdr.deserialize::<ProcessedFragment>()
        .map(|r| r.expect("malformed golden CSV row"))
        .collect()
}

/// Sort key for stable comparison: (start, end, reverse).
fn sort_key(f: &ProcessedFragment) -> (usize, usize, u8) {
    (f.start, f.end, f.reverse)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn digest_ecoli_plasmid_matches_golden() {
    let root = repo_root();
    let fasta = root.join("test/test_project4/genomes/Escherichia_coli_plasmid.fasta");
    let golden_csv = root.join(
        "test/test_project4/golden_output/raw_digests/raw_digest_Escherichia_coli_plasmid.fasta.csv",
    );

    assert!(fasta.exists(), "FASTA not found: {}", fasta.display());
    assert!(golden_csv.exists(), "Golden CSV not found: {}", golden_csv.display());

    let config = build_hhai_msei_config();
    // max_len = u + 6*sd = 300 + 6*100 = 900
    let max_len = 900usize;

    let mut rust_frags = digest_fasta(&fasta, &config, max_len)
        .expect("digest_fasta failed");
    let mut golden_frags = load_golden(&golden_csv);

    // Sort both by (start, end, reverse) for order-independent comparison
    rust_frags.sort_by_key(sort_key);
    golden_frags.sort_by_key(sort_key);

    assert_eq!(
        rust_frags.len(),
        golden_frags.len(),
        "Fragment count mismatch: Rust={}, Golden={}",
        rust_frags.len(),
        golden_frags.len()
    );

    for (i, (r, g)) in rust_frags.iter().zip(golden_frags.iter()).enumerate() {
        assert_eq!(r.chrom,    g.chrom,    "row {i}: chrom mismatch");
        assert_eq!(r.seq,      g.seq,      "row {i}: seq mismatch (start={}, end={})", g.start, g.end);
        assert_eq!(r.start,    g.start,    "row {i}: start mismatch");
        assert_eq!(r.end,      g.end,      "row {i}: end mismatch");
        assert_eq!(r.m1,       g.m1,       "row {i}: m1 mismatch");
        assert_eq!(r.m2,       g.m2,       "row {i}: m2 mismatch");
        assert_eq!(r.internal, g.internal, "row {i}: internal mismatch");
        assert_eq!(r.reverse,  g.reverse,  "row {i}: reverse mismatch");
        assert_eq!(r.revc,     g.revc,     "row {i}: revc mismatch");
        assert_eq!(r.length,   g.length,   "row {i}: length mismatch");
    }
}

/// Smoke test on the tiny synthetic genome from test_project3
/// (two fragments, hand-verifiable).
#[test]
fn digest_fake_genome_1_smoke() {
    let root = repo_root();
    let fasta = root.join("test/test_project3/genomes/fake_genome_1.fasta");

    if !fasta.exists() {
        eprintln!("Skipping smoke test: {}", fasta.display());
        return;
    }

    let config = build_hhai_msei_config();
    let frags = digest_fasta(&fasta, &config, 900).expect("digest failed");

    // Two RE sites produce exactly 2 fragments (GCGC→TTAA forward, TTAA→GCGC reverse).
    // Both have 21 N's between motifs (genome is symmetric), so both trim to length 23.
    assert_eq!(frags.len(), 2, "Expected 2 fragments, got {}", frags.len());
    assert!(frags.iter().all(|f| f.length == 23), "Expected all lengths=23, got {:?}",
        frags.iter().map(|f| f.length).collect::<Vec<_>>());
}
