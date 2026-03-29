use anyhow::{bail, Result};
use rand::seq::SliceRandom;
use std::fs;
use std::io::Write;
use std::path::Path;

use crate::types::CountFragment;

/// One sequencing adapter: (top_strand, bottom_strand, identifier).
pub type Adapter = (String, String, String);

/// Configuration for fastq read output.
pub struct ReadConfig {
    /// Pool of R1-side adapters (randomly drawn each read).
    pub adapters1: Vec<Adapter>,
    /// Pool of R2-side adapters (randomly drawn each read).
    pub adapters2: Vec<Adapter>,
    /// Offset into the R1 ligated string where sequencing begins.
    pub a1s: usize,
    /// Offset into the R2 ligated string where sequencing begins.
    pub a2s: usize,
    /// R1 read length.
    pub l1: usize,
    /// R2 read length.
    pub l2: usize,
    /// Filler base when fragment + adapter is shorter than the read length.
    pub filler_base: char,
    /// Sampled R1 quality score strings. If None, quality is 'I' × l1.
    pub q1_scores: Option<Vec<String>>,
    /// Sampled R2 quality score strings. If None, quality is 'I' × l2.
    pub q2_scores: Option<Vec<String>>,
}

impl ReadConfig {
    /// Build the default adapter configuration used when no custom adapters are provided.
    ///
    /// Mirrors Python's create_adapters() + the a1s/a2s = 62/58 assignment in readsynth.py.
    pub fn default_adapters(l1: usize, l2: Option<usize>, filler_base: char) -> Self {
        let l2 = l2.unwrap_or(l1);
        ReadConfig {
            adapters1: vec![(
                "AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG".into(),
                "CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT".into(),
                "rs1".into(),
            )],
            adapters2: vec![(
                "CAAGCAGAAGACGGCATACGAGATGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG".into(),
                "CTGTCTCTTATACACATCTCCGAGCCCACGAGACATCTCGTATGCCGTCTTCTGCTTG".into(),
                "rs2".into(),
            )],
            a1s: 62,
            a2s: 58,
            l1,
            l2,
            filler_base,
            q1_scores: None,
            q2_scores: None,
        }
    }

    /// Build iso-enzyme default adapters, padding each adapter with N×cut_pos1
    /// overhangs to reflect the single-enzyme two-cut-site geometry.
    ///
    /// Mirrors Python's create_adapters_iso().
    pub fn iso_adapters(
        l1: usize,
        l2: Option<usize>,
        filler_base: char,
        m1_cut1: usize,
        m2_cut1: usize,
    ) -> Self {
        let l2 = l2.unwrap_or(l1);
        let n1 = "N".repeat(m1_cut1);
        let n2 = "N".repeat(m2_cut1);
        ReadConfig {
            adapters1: vec![(
                format!("AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG{n1}"),
                format!("{n2}CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATT"),
                "rs1".into(),
            )],
            adapters2: vec![(
                format!("CAAGCAGAAGACGGCATACGAGATGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG{n1}"),
                format!("{n2}CTGTCTCTTATACACATCTCCGAGCCCACGAGACATCTCGTATGCCGTCTTCTGCTTG"),
                "rs2".into(),
            )],
            a1s: 62,
            a2s: 58,
            l1,
            l2,
            filler_base,
            q1_scores: None,
            q2_scores: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Adapter + quality score file loading
// ---------------------------------------------------------------------------

/// Load a custom adapter file: whitespace-separated `top_seq  bot_seq  id` per line.
///
/// Mirrors Python's get_adapters().
pub fn load_adapters(path: &Path) -> Result<Vec<Adapter>> {
    let content = fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("cannot open adapter file {}: {e}", path.display()))?;
    let mut adapters: Vec<Adapter> = Vec::new();
    for (i, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            bail!(
                "adapter file line {}: expected 'top bot id', got: {:?}",
                i + 1,
                line
            );
        }
        adapters.push((parts[0].to_string(), parts[1].to_string(), parts[2].to_string()));
    }
    if adapters.is_empty() {
        bail!("no adapters found in {}", path.display());
    }
    Ok(adapters)
}

/// Load a quality score file: one Phred quality string per line (no header).
///
/// Mirrors Python's `pd.read_csv(args.q1, names=['score'], sep='\t')`.
pub fn load_q_scores(path: &Path) -> Result<Vec<String>> {
    let content = fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("cannot open q-score file {}: {e}", path.display()))?;
    let scores: Vec<String> = content
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    if scores.is_empty() {
        bail!("no quality score strings found in {}", path.display());
    }
    Ok(scores)
}

// ---------------------------------------------------------------------------
// Read writing
// ---------------------------------------------------------------------------

/// Write paired-end fastq reads for one genome's count fragments.
///
/// Mirrors Python's write_reads.read_writer_basic() (no q-scores) and
/// read_writer_samples() (with q-scores).
///
/// R1 ligated string: a1_top + seq  + a2_bot
/// R2 ligated string: a2_top + revc + a1_bot
/// Both are sliced to [a?s .. a?s+l?] and right-padded with filler_base.
///
/// When q-scores are provided, a score string is sampled per read and both
/// the sequence and score are truncated to min(seq_len, score_len), mirroring
/// Python's even_score_seq().
///
/// Header format: @{r_no}:{idx}:{full}:{length}:{a1_id}:{a2_id}:{genome_name}
pub fn write_reads(
    fragments: &[CountFragment],
    r1: &mut dyn Write,
    r2: &mut dyn Write,
    genome_name: &str,
    config: &ReadConfig,
) -> Result<()> {
    let mut rng = rand::thread_rng();
    let a1e = config.a1s + config.l1;
    let a2e = config.a2s + config.l2;
    let filler = config.filler_base;

    // Pre-build default quality strings (used when no q-score files are configured).
    let default_score1 = "I".repeat(config.l1);
    let default_score2 = "I".repeat(config.l2);

    let mut r_no: usize = 0;

    for (idx, frag) in fragments.iter().enumerate() {
        for _ in 0..frag.counts {
            let a1 = config.adapters1.choose(&mut rng).expect("adapters1 empty");
            let a2 = config.adapters2.choose(&mut rng).expect("adapters2 empty");

            let r1_lig = format!("{}{}{}", a1.0, frag.fragment.seq, a2.1);
            let r2_lig = format!("{}{}{}", a2.0, frag.fragment.revc, a1.1);
            let full = r1_lig.len();

            let mut r1_read = slice_pad(&r1_lig, config.a1s, a1e, config.l1, filler);
            let mut r2_read = slice_pad(&r2_lig, config.a2s, a2e, config.l2, filler);

            // Apply quality scores, truncating both seq and score to the shorter length.
            let score1 = match &config.q1_scores {
                Some(pool) => {
                    let s = pool.choose(&mut rng).expect("q1_scores empty");
                    let len = r1_read.len().min(s.len());
                    r1_read.truncate(len);
                    s[..len].to_string()
                }
                None => default_score1.clone(),
            };
            let score2 = match &config.q2_scores {
                Some(pool) => {
                    let s = pool.choose(&mut rng).expect("q2_scores empty");
                    let len = r2_read.len().min(s.len());
                    r2_read.truncate(len);
                    s[..len].to_string()
                }
                None => default_score2.clone(),
            };

            let header = format!(
                "@{r_no}:{idx}:{full}:{}:{}:{}:{genome_name}",
                frag.fragment.length, a1.2, a2.2
            );

            writeln!(r1, "{} 1\n{}\n+\n{}", header, r1_read, score1)?;
            writeln!(r2, "{} 2\n{}\n+\n{}", header, r2_read, score2)?;
            r_no += 1;
        }
    }

    Ok(())
}

/// Slice `s[start..end]` and right-pad to `len` with `filler`.
fn slice_pad(s: &str, start: usize, end: usize, len: usize, filler: char) -> String {
    let bytes = s.as_bytes();
    let slice_start = start.min(bytes.len());
    let slice_end = end.min(bytes.len());
    let taken = &bytes[slice_start..slice_end];
    let mut out = String::with_capacity(len);
    out.push_str(std::str::from_utf8(taken).unwrap());
    while out.len() < len {
        out.push(filler);
    }
    out
}
