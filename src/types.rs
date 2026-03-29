use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

/// Restriction enzyme motif configuration.
///
/// Mirrors the Python args.motif_dt / motif_dt1 / motif_dt2 / motif_len namespace.
/// IndexMap is used (not HashMap) because the trimming in process_fragments applies
/// motifs in insertion order — GCGC before TTAA — and that order affects the result.
#[derive(Debug, Clone)]
pub struct MotifConfig {
    /// All motifs (m1 ∪ m2): compiled regex pattern → cut position.
    pub motif_dt: IndexMap<String, usize>,
    /// m1 (forward adapter) motifs: compiled regex pattern → cut position.
    pub motif_dt1: IndexMap<String, usize>,
    /// m2 (reverse adapter) motifs: compiled regex pattern → cut position.
    pub motif_dt2: IndexMap<String, usize>,
    /// regex pattern → actual bp length of the motif (bracket groups count as 1 bp).
    pub motif_len: HashMap<String, usize>,
}

/// A raw fragment produced by digest_seq/digest_frag, before process_fragments filtering.
/// Mirrors the intermediate DataFrame from Python's digest_genomes.main().
#[derive(Debug, Clone)]
pub struct RawFragment {
    pub chrom: String,
    pub seq: String,
    pub start: usize,
    pub end: usize,
    pub m1: String,
    pub m2: String,
    pub internal: usize,
}

/// A fully processed fragment — this is what gets written to raw_digest_*.csv.
/// Column order matches the Python golden output exactly.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ProcessedFragment {
    pub chrom: String,
    pub seq: String,
    pub start: usize,
    pub end: usize,
    pub m1: String,
    pub m2: String,
    pub internal: usize,
    pub reverse: u8,
    pub revc: String,
    pub length: usize,
}

/// One entry from the genome abundance CSV (-g flag).
#[derive(Debug, Clone)]
pub struct GenomeEntry {
    pub path: PathBuf,
    /// Normalised fractional abundance (sums to 1.0 across all genomes).
    pub abundance: f64,
}

/// Per-fragment probability and read-count data, written to counts_*.csv.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CountFragment {
    #[serde(flatten)]
    pub fragment: ProcessedFragment,
    pub probability: f64,
    pub adj_prob: f64,
    pub counts: u64,
}

/// Iso-enzyme (Type IIB) motif configuration.
/// Each recognition site has two cut positions instead of one.
#[derive(Debug, Clone)]
pub struct MotifConfigIso {
    /// regex pattern → [cut_pos1, cut_pos2] (positions in the stripped-slash sequence)
    pub motif_dt: IndexMap<String, [usize; 2]>,
    /// regex pattern → bp length of full recognition sequence (no slashes)
    pub motif_len: HashMap<String, usize>,
}

/// Global simulation configuration (mirrors the Python args namespace).
#[derive(Debug, Clone)]
pub struct Config {
    pub output_dir: PathBuf,
    pub n: usize,
    pub c: f64,
    pub l1: usize,
    pub l2: usize,
    pub max_len: usize,
    /// Gaussian mean for size selection (None when using low-pass or custom dist).
    pub u: Option<usize>,
    /// Gaussian SD for size selection.
    pub sd: Option<usize>,
    /// Intersection point for Gaussian/fragment distribution scaling.
    pub x: Option<usize>,
    /// Path to custom fragment-length JSON distribution.
    pub custom_dist: Option<PathBuf>,
    /// Low-pass max fragment length (disables Gaussian size selection).
    pub low_pass: Option<usize>,
    pub filler_base: char,
    pub test_mode: bool,
    pub motifs: MotifConfig,
    pub is_iso: bool,
}
