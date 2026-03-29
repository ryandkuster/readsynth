/// Integration test: run digest + apply_probabilities on the E. coli plasmid and compare
/// probability / adj_prob against the golden counts CSV.
///
/// Golden parameters: c=0.5, abundance=1.0 (single genome, comp=1.0)
/// Expected: probability = 0.25 * 0.5^internal, adj_prob = probability * 1.0
use readsynth::digest::digest_fasta;
use readsynth::iupac::{check_for_palindrome, iupac_motifs, regex_bp_len};
use readsynth::prob_copies::apply_probabilities;
use readsynth::types::{CountFragment, MotifConfig};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

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

/// Load probability / adj_prob from the golden counts CSV, keyed by (start, end, reverse).
/// Also returns the internal value so we can double-check the formula.
fn load_golden_probs(path: &Path) -> HashMap<(usize, usize, u8), (usize, f64, f64)> {
    let mut rdr = csv::Reader::from_path(path)
        .unwrap_or_else(|e| panic!("Cannot open golden counts CSV {}: {e}", path.display()));

    // Headers: (empty), chrom, seq, start, end, m1, m2, internal, reverse, revc, length,
    //          probability, adj_prob, counts
    let headers = rdr.headers().unwrap().clone();
    let col = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .unwrap_or_else(|| panic!("column '{name}' not found in golden counts CSV"))
    };
    let i_start    = col("start");
    let i_end      = col("end");
    let i_internal = col("internal");
    let i_reverse  = col("reverse");
    let i_prob     = col("probability");
    let i_adj_prob = col("adj_prob");

    let mut map = HashMap::new();
    for result in rdr.records() {
        let rec = result.unwrap();
        let start:    usize = rec[i_start].parse().unwrap();
        let end:      usize = rec[i_end].parse().unwrap();
        let internal: usize = rec[i_internal].parse().unwrap();
        let reverse:  u8    = rec[i_reverse].parse().unwrap();
        let prob:     f64   = rec[i_prob].parse().unwrap();
        let adj_prob: f64   = rec[i_adj_prob].parse().unwrap();
        map.insert((start, end, reverse), (internal, prob, adj_prob));
    }
    map
}

#[test]
fn prob_copies_ecoli_plasmid_matches_golden() {
    let root = repo_root();
    let fasta = root.join("test/test_project4/genomes/Escherichia_coli_plasmid.fasta");
    let golden_csv = root.join(
        "test/test_project4/golden_output/individual_counts/counts_Escherichia_coli_plasmid.fasta.csv",
    );

    assert!(fasta.exists(), "FASTA not found: {}", fasta.display());
    assert!(golden_csv.exists(), "Golden counts CSV not found: {}", golden_csv.display());

    let config  = build_hhai_msei_config();
    let max_len = 900usize;

    let fragments = digest_fasta(&fasta, &config, max_len).expect("digest_fasta failed");
    let (rust_counts, _len_freqs) = apply_probabilities(fragments, 0.5, 1.0, max_len);

    // Index rust output by (start, end, reverse)
    let rust_map: HashMap<(usize, usize, u8), &CountFragment> = rust_counts
        .iter()
        .map(|f| ((f.fragment.start, f.fragment.end, f.fragment.reverse), f))
        .collect();

    let golden = load_golden_probs(&golden_csv);

    // For every fragment in the golden, check our probability and adj_prob match.
    for ((start, end, reverse), (internal, g_prob, g_adj_prob)) in &golden {
        let r = rust_map
            .get(&(*start, *end, *reverse))
            .unwrap_or_else(|| {
                panic!(
                    "Rust output missing fragment (start={start}, end={end}, reverse={reverse})"
                )
            });

        let eps = 1e-12;
        assert!(
            (r.probability - g_prob).abs() < eps,
            "probability mismatch at (start={start}, end={end}, reverse={reverse}, internal={internal}): \
             rust={}, golden={g_prob}",
            r.probability
        );
        assert!(
            (r.adj_prob - g_adj_prob).abs() < eps,
            "adj_prob mismatch at (start={start}, end={end}, reverse={reverse}): \
             rust={}, golden={g_adj_prob}",
            r.adj_prob
        );
    }

    // Also verify every rust fragment that matches a golden key satisfies the formula
    // probability = 0.5² × 0.5^internal, just as a sanity check.
    for f in &rust_counts {
        let expected = 0.5f64.powi(2) * 0.5f64.powi(f.fragment.internal as i32);
        let bidirectional = rust_counts
            .iter()
            .filter(|g| {
                g.fragment.start == f.fragment.start && g.fragment.end == f.fragment.end
                    && g.fragment.reverse != f.fragment.reverse
            })
            .count()
            > 0;
        let expected = if bidirectional { expected / 2.0 } else { expected };
        let eps = 1e-12;
        assert!(
            (f.probability - expected).abs() < eps,
            "formula mismatch: got {}, expected {expected} (internal={}, bidirectional={bidirectional})",
            f.probability,
            f.fragment.internal
        );
    }

    println!(
        "prob_copies: {} rust fragments, {} golden rows — all probability/adj_prob match",
        rust_counts.len(),
        golden.len()
    );
}
