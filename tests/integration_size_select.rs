/// Integration test: run digest → apply_probabilities → size_select_gaussian and compare
/// per-length read counts against the golden combined.csv.
///
/// Golden parameters: c=0.5, n=10000, u=300, sd=100, x=u=300, max_len=900, abundance=1.0
/// The `counts` column in combined.csv = round(fragment_comps[L] * sum_prob[L] * adjustment)
use readsynth::digest::digest_fasta;
use readsynth::iupac::{check_for_palindrome, iupac_motifs, regex_bp_len};
use readsynth::prob_copies::apply_probabilities;
use readsynth::size_select::size_select_gaussian;
use readsynth::types::MotifConfig;
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

/// Load (length → (sum_prob, golden_counts)) from the combined.csv,
/// returning only rows where sum_prob > 0.
fn load_golden_combined(path: &Path) -> Vec<(usize, f64, f64)> {
    let mut rdr = csv::Reader::from_path(path)
        .unwrap_or_else(|e| panic!("Cannot open golden combined.csv {}: {e}", path.display()));

    // Header: (empty), length, sum_prob, name, counts_file, counts
    let headers = rdr.headers().unwrap().clone();
    let col = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .unwrap_or_else(|| panic!("column '{name}' not found in combined.csv"))
    };
    let i_length   = col("length");
    let i_sum_prob = col("sum_prob");
    let i_counts   = col("counts");

    rdr.records()
        .map(|r| r.unwrap())
        .filter_map(|rec| {
            let length:   usize = rec[i_length].parse().unwrap();
            let sum_prob: f64   = rec[i_sum_prob].parse().unwrap();
            let counts:   f64   = rec[i_counts].parse().unwrap();
            if sum_prob > 0.0 { Some((length, sum_prob, counts)) } else { None }
        })
        .collect()
}

#[test]
fn size_select_gaussian_ecoli_plasmid_matches_golden() {
    let root = repo_root();
    let fasta = root.join("test/test_project4/genomes/Escherichia_coli_plasmid.fasta");
    let golden_combined =
        root.join("test/test_project4/golden_output/combined.csv");

    assert!(fasta.exists(), "FASTA not found: {}", fasta.display());
    assert!(golden_combined.exists(), "Golden combined.csv not found: {}", golden_combined.display());

    let config  = build_hhai_msei_config();
    let max_len = 900usize;
    let n       = 10_000usize;
    let u       = 300usize;
    let sd      = 100usize;
    let x       = u; // default: x = u when -x not specified

    let fragments = digest_fasta(&fasta, &config, max_len).expect("digest_fasta failed");
    let (_counts, len_freqs) = apply_probabilities(fragments, 0.5, 1.0, max_len);

    let result = size_select_gaussian(&len_freqs, n, u, sd, x, max_len)
        .expect("size_select_gaussian failed");

    let golden = load_golden_combined(&golden_combined);

    let mut mismatches = 0usize;
    for (length, sum_prob, golden_counts) in &golden {
        let fc = *result.fragment_comps.get(length).unwrap_or(&0.0);
        // Python: round(fragment_comps[L] * sum_prob[L] * adjustment)
        // Python's round() uses banker's rounding for .5; Rust round() rounds away from zero.
        // For values far from .5 this makes no difference.
        let rust_counts = (fc * sum_prob * result.adjustment).round();
        if (rust_counts - golden_counts).abs() > 0.5 {
            eprintln!(
                "MISMATCH length={length}: rust={rust_counts}, golden={golden_counts} \
                 (fc={fc:.6}, sum_prob={sum_prob}, adj={:.6})",
                result.adjustment
            );
            mismatches += 1;
        }
    }

    assert_eq!(
        mismatches, 0,
        "{mismatches} length(s) had count mismatches (see stderr for details)"
    );

    let total_golden: f64 = golden.iter().map(|(_, _, c)| c).sum();
    println!(
        "size_select: {} non-zero lengths, total golden reads={total_golden}, adjustment={:.6}",
        golden.len(),
        result.adjustment
    );
}
