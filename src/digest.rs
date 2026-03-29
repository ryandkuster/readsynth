use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::iupac::reverse_comp;
use crate::types::{MotifConfig, MotifConfigIso, ProcessedFragment, RawFragment};

// ---------------------------------------------------------------------------
// Regex helpers
// ---------------------------------------------------------------------------

/// Pre-compile one Regex per motif pattern. Panics on invalid patterns.
fn compile_regexes(motif_len: &HashMap<String, usize>) -> HashMap<String, Regex> {
    motif_len
        .keys()
        .map(|p| (p.clone(), Regex::new(p).expect("invalid motif regex")))
        .collect()
}

/// Find all (possibly overlapping) start positions of `re` within `text`.
///
/// Equivalent to Python's `[m.start() for m in re.finditer('(?=' + pat + ')', text)]`.
/// The lookahead trick in Python allows overlapping hits; here we replicate it by
/// advancing one byte at a time after each match.
fn overlapping_positions(re: &Regex, text: &str) -> Vec<usize> {
    let mut positions = Vec::new();
    let mut start = 0usize;
    while start < text.len() {
        match re.find(&text[start..]) {
            Some(m) => {
                let abs = start + m.start();
                positions.push(abs);
                start = abs + 1;
            }
            None => break,
        }
    }
    positions
}

// ---------------------------------------------------------------------------
// Internal-site counting
// ---------------------------------------------------------------------------

/// Count all RE recognition sites in `subseq` across all motifs.
/// Mirrors Python's internal_sites().
fn count_internal_sites(subseq: &str, motif_regexes: &HashMap<String, Regex>) -> usize {
    motif_regexes
        .values()
        .map(|re| overlapping_positions(re, subseq).len())
        .sum()
}

// ---------------------------------------------------------------------------
// Digest: per-chromosome
// ---------------------------------------------------------------------------

/// Search one chromosome for all valid RE fragments.
/// Mirrors Python's digest_seq().
fn digest_seq(
    chrom: &str,
    seq: &str,
    config: &MotifConfig,
    motif_regexes: &HashMap<String, Regex>,
    max_len: usize,
) -> Vec<RawFragment> {
    let mut fragments = Vec::new();

    // Iterate over motif_dt (IndexMap) to get deterministic insertion order,
    // matching Python's `for motif1 in args.motif_len.keys()` which preserves
    // the GCGC-before-TTAA order from the ordered motif_dt dict.
    for motif1 in config.motif_dt.keys() {
        let re1 = &motif_regexes[motif1];
        for start in overlapping_positions(re1, seq) {
            let window_end = (start + max_len).min(seq.len());
            let fragment = &seq[start..window_end];
            fragments.extend(digest_frag(chrom, fragment, motif1, start, config, motif_regexes));
        }
    }

    fragments
}

/// From a fragment window (starting at m1), find all valid m2 endpoints.
/// Mirrors Python's digest_frag().
fn digest_frag(
    chrom: &str,
    fragment: &str,
    motif1: &str,
    f_start: usize,
    config: &MotifConfig,
    motif_regexes: &HashMap<String, Regex>,
) -> Vec<RawFragment> {
    let mut result = Vec::new();

    // Search starts at position 1 (skip first char) so we don't re-match motif1
    // at its own start position. Mirrors Python's `fragment[1:]`.
    if fragment.len() < 2 {
        return result;
    }
    let search_region = &fragment[1..];

    // Use motif_dt (IndexMap) for deterministic insertion order, matching Python.
    for motif2 in config.motif_dt.keys() {
        let re2 = &motif_regexes[motif2];
        let m2_len = config.motif_len[motif2];

        for rel_pos in overlapping_positions(re2, search_region) {
            // rel_pos is relative to search_region (= fragment[1:])
            // end_rel is the m2 start position relative to fragment
            let end_rel = rel_pos + 1;
            let seq_end = end_rel + m2_len;

            if seq_end > fragment.len() {
                continue;
            }

            // Internal region: fragment[1 .. seq_end - 1]
            // Mirrors Python: fragment[1 : end + motif_len[m2] - 1]
            let internal_region = if seq_end > 1 {
                &fragment[1..seq_end - 1]
            } else {
                ""
            };
            let internals = count_internal_sites(internal_region, motif_regexes);

            result.push(RawFragment {
                chrom: chrom.to_string(),
                seq: fragment[..seq_end].to_string(),
                start: f_start,
                end: f_start + end_rel,
                m1: motif1.to_string(),
                m2: motif2.to_string(),
                internal: internals,
            });
        }
    }

    result
}

// ---------------------------------------------------------------------------
// FASTA reading
// ---------------------------------------------------------------------------

/// Open a FASTA file, handling both plain text and gzip (.gz) transparently.
fn open_fasta(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    if path.extension().map_or(false, |e| e == "gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Read a FASTA file and return (chromosome_name, sequence) pairs.
fn read_fasta(path: &Path) -> Result<Vec<(String, String)>> {
    let reader = open_fasta(path)?;
    let mut chromosomes: Vec<(String, String)> = Vec::new();
    let mut current_chrom = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                chromosomes.push((current_chrom.clone(), current_seq.clone()));
                current_seq.clear();
            }
            // Take first word after '>', strip leading '>'
            current_chrom = line
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
        } else {
            current_seq.push_str(line.trim().to_ascii_uppercase().as_str());
        }
    }

    if !current_seq.is_empty() {
        chromosomes.push((current_chrom, current_seq));
    }

    Ok(chromosomes)
}

// ---------------------------------------------------------------------------
// process_fragments — equivalent of Python's process_df()
// ---------------------------------------------------------------------------

/// Replace IUPAC ambiguity codes with 'N', leaving A/C/G/T/N unchanged.
/// Mirrors the regex substitution `df['seq'].str.replace('[RYSWKMBDHV]', 'N')`.
fn to_n(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' => 'N',
            _ => c,
        })
        .collect()
}

/// Compute the trimmed fragment length without permanently modifying seq/revc.
///
/// This mirrors the length-calculation block in Python's process_df(), which:
///   1. saves seq/revc as backups
///   2. applies overhang trimming in motif_dt iteration order
///   3. measures len(seq)
///   4. restores seq/revc from backups
///
/// The trimming order follows config.motif_dt (IndexMap, preserving insertion order).
/// For the standard test case this is GCGC first, TTAA second — that order is critical.
///
/// `seq` here is already the (potentially swapped) sequence as it will appear in the
/// output CSV, i.e. the reverse-complement for reverse=1 fragments.
fn compute_trimmed_length(
    seq: &str,
    m1: &str,
    m2: &str,
    reverse: u8,
    config: &MotifConfig,
) -> usize {
    let mut trimmed = seq.to_string();

    for (mot, &front) in &config.motif_dt {
        let bp_len = config.motif_len[mot];
        let back = bp_len - front; // bases on the 3' side of the cut

        // --- m1 trimming ---
        if mot == m1 {
            if reverse == 0 {
                // Forward: trim `front` bases from the 5' end of seq
                // Mirrors: df.loc[(m1==mot)&(rev==0), 'seq'] = seq.str[front:]
                trimmed = trimmed[front.min(trimmed.len())..].to_string();
            } else {
                // Reverse: trim `back` bases from the 3' end
                // Mirrors: df.loc[(m1==mot)&(rev==1), 'seq'] = seq.str[:-back]
                if back > 0 {
                    let new_len = trimmed.len().saturating_sub(back);
                    trimmed = trimmed[..new_len].to_string();
                }
                // back==0 → Python: seq[:-0] == seq[:0] == ""; guard avoids that footgun.
            }
        }

        // --- m2 trimming ---
        if mot == m2 {
            if reverse == 0 {
                // Forward: trim `back` bases from the 3' end
                // Mirrors: df.loc[(m2==mot)&(rev==0), 'seq'] = seq.str[:-back]
                if back > 0 {
                    let new_len = trimmed.len().saturating_sub(back);
                    trimmed = trimmed[..new_len].to_string();
                }
            } else {
                // Reverse: trim `front` bases from the 5' end
                // Mirrors: df.loc[(m2==mot)&(rev==1), 'seq'] = seq.str[front:]
                trimmed = trimmed[front.min(trimmed.len())..].to_string();
            }
        }
        // When m1 == m2 (homoenzyme), both blocks execute in the same iteration,
        // chaining: the m2 block sees the result of the m1 block.
    }

    trimmed.len()
}

/// Build a ProcessedFragment from a RawFragment plus its orientation.
fn make_processed(
    frag: &RawFragment,
    seq: &str,
    revc: &str,
    reverse: u8,
    config: &MotifConfig,
) -> ProcessedFragment {
    let length = compute_trimmed_length(seq, &frag.m1, &frag.m2, reverse, config);
    ProcessedFragment {
        chrom: frag.chrom.clone(),
        seq: seq.to_string(),
        start: frag.start,
        end: frag.end,
        m1: frag.m1.clone(),
        m2: frag.m2.clone(),
        internal: frag.internal,
        reverse,
        revc: revc.to_string(),
        length,
    }
}

/// Filter and orient raw fragments; compute trimmed lengths; sort by length.
/// Mirrors Python's process_df() in readsynth.py.
///
/// Key behaviours replicated:
///   - forward/reverse classification via motif_dt1 / motif_dt2 membership
///   - removal of fragments where motifs overlap (seq.len <= m1_len + m2_len)
///   - IUPAC → N conversion
///   - reverse-complement computation
///   - bidirectional fragments duplicated into both orientations
///   - seq/revc swapped for reverse-strand fragments
///   - trimmed length computed (seq/revc in output remain full/untruncated)
///   - sorted ascending by trimmed length
pub fn process_fragments(raw: Vec<RawFragment>, config: &MotifConfig) -> Vec<ProcessedFragment> {
    let mut result: Vec<ProcessedFragment> = Vec::new();

    for frag in raw {
        let is_forward = config.motif_dt1.contains_key(&frag.m1)
            && config.motif_dt2.contains_key(&frag.m2);
        let is_reverse = config.motif_dt2.contains_key(&frag.m1)
            && config.motif_dt1.contains_key(&frag.m2);

        if !is_forward && !is_reverse {
            continue;
        }

        // Drop fragments whose motifs overlap (can't form a valid insert)
        let min_len = config.motif_len[&frag.m1] + config.motif_len[&frag.m2];
        if frag.seq.len() <= min_len {
            continue;
        }

        // Convert IUPAC ambiguity codes to N, then compute reverse complement
        let seq_n = to_n(&frag.seq);
        let revc = reverse_comp(&seq_n);

        if is_forward && is_reverse {
            // Bidirectional: emit both orientations.
            // Forward copy keeps seq as-is.
            result.push(make_processed(&frag, &seq_n, &revc, 0, config));
            // Reverse copy swaps seq ↔ revc (seq becomes the reverse complement).
            result.push(make_processed(&frag, &revc, &seq_n, 1, config));
        } else if is_forward {
            result.push(make_processed(&frag, &seq_n, &revc, 0, config));
        } else {
            // Reverse-only: swap seq and revc before storing.
            result.push(make_processed(&frag, &revc, &seq_n, 1, config));
        }
    }

    // Sort ascending by trimmed length, matching Python's df.sort_values('length')
    result.sort_by_key(|f| f.length);
    result
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Read a FASTA file, digest all chromosomes with the given motif config,
/// apply process_fragments, and return the sorted processed fragments.
///
/// This is the combined equivalent of Python's:
///   digest_genomes.main(args)  →  process_df(df, digest_file, args)
pub fn digest_fasta(
    path: &Path,
    config: &MotifConfig,
    max_len: usize,
) -> Result<Vec<ProcessedFragment>> {
    let motif_regexes = compile_regexes(&config.motif_len);
    let chromosomes = read_fasta(path)?;

    let mut raw: Vec<RawFragment> = Vec::new();
    for (chrom, seq) in &chromosomes {
        raw.extend(digest_seq(chrom, seq, config, &motif_regexes, max_len));
    }

    Ok(process_fragments(raw, config))
}

// ---------------------------------------------------------------------------
// CSV output
// ---------------------------------------------------------------------------

/// Write processed fragments to a CSV file, matching the Python column order:
/// chrom, seq, start, end, m1, m2, internal, reverse, revc, length
pub fn write_raw_digest(
    fragments: &[ProcessedFragment],
    path: &Path,
) -> Result<()> {
    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record(["chrom","seq","start","end","m1","m2","internal","reverse","revc","length"])?;
    for f in fragments {
        wtr.write_record(&[
            &f.chrom,
            &f.seq,
            &f.start.to_string(),
            &f.end.to_string(),
            &f.m1,
            &f.m2,
            &f.internal.to_string(),
            &f.reverse.to_string(),
            &f.revc,
            &f.length.to_string(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Iso-enzyme digest
// ---------------------------------------------------------------------------

/// Read a FASTA file and digest all chromosomes using an iso-enzyme (Type IIB) motif.
///
/// Each recognition site produces one fragment spanning from cut_pos1 to cut_pos2
/// within the match. Both forward and RC motif variants are searched.
///
/// Mirrors Python's digest_genomes_iso.main() combined with process_df_iso().
pub fn digest_fasta_iso(
    path: &Path,
    config: &MotifConfigIso,
    max_len: usize,
) -> Result<Vec<ProcessedFragment>> {
    let regexes: HashMap<String, Regex> = config
        .motif_len
        .keys()
        .map(|p| (p.clone(), Regex::new(p).expect("invalid iso motif regex")))
        .collect();

    let chromosomes = read_fasta(path)?;
    let motif_keys: Vec<&String> = config.motif_dt.keys().collect();

    let mut fragments: Vec<ProcessedFragment> = Vec::new();
    let mut genome_offset: usize = 0;

    for (chrom, seq) in &chromosomes {
        for (mi, motif) in motif_keys.iter().enumerate() {
            let re = &regexes[*motif];
            let cuts = config.motif_dt[*motif]; // [cut1, cut2]
            let ml = config.motif_len[*motif];
            let is_rc = mi == 1; // second entry is the RC variant

            for pos in overlapping_positions(re, seq) {
                let match_end = pos + ml;
                if match_end > seq.len() {
                    continue;
                }
                let cut1 = cuts[0];
                let cut2 = cuts[1].min(ml);
                if cut2 <= cut1 {
                    continue;
                }
                let length = cut2 - cut1;
                if length == 0 || length > max_len {
                    continue;
                }

                let seq_n = to_n(&seq[pos + cut1..pos + cut2]);
                let revc = reverse_comp(&seq_n);

                // RC motif: swap seq ↔ revc, mirrors process_df_iso swap.
                let (final_seq, final_revc, reverse) = if is_rc {
                    (revc, seq_n, 1u8)
                } else {
                    (seq_n, revc, 0u8)
                };

                // Python stores start = match_start, end = match_start + motif_len
                // (the full recognition-sequence coordinates), not the trimmed cut
                // coordinates. The overlap-detection algorithm in prob_n_copies_iso
                // expects these full-sequence coordinates.
                fragments.push(ProcessedFragment {
                    chrom: chrom.clone(),
                    seq: final_seq,
                    start: genome_offset + pos,
                    end: genome_offset + pos + ml,
                    m1: motif.to_string(),
                    m2: motif.to_string(),
                    internal: 0,
                    reverse,
                    revc: final_revc,
                    length,
                });
            }
        }
        genome_offset += seq.len();
    }

    // Python's process_df_iso sorts by 'start' (genomic position).
    // The calculate_prob_iso cluster-detection algorithm requires this ordering.
    fragments.sort_by_key(|f| f.start);
    Ok(fragments)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::iupac::{iupac_motifs, regex_bp_len};
    use indexmap::IndexMap;
    use std::collections::HashMap;

    /// Build a MotifConfig for the palindrome homoenzyme case (m1 == m2 == GCGC).
    fn hhai_hhai_config() -> MotifConfig {
        let motif_dt1 = iupac_motifs(&["GCG/C".to_string()]);
        let motif_dt2 = iupac_motifs(&["GCG/C".to_string()]);
        let mut motif_dt = motif_dt1.clone();
        for (k, v) in &motif_dt2 {
            motif_dt.insert(k.clone(), *v);
        }
        let motif_len: HashMap<_, _> = motif_dt
            .keys()
            .map(|p| (p.clone(), regex_bp_len(p)))
            .collect();
        MotifConfig { motif_dt, motif_dt1, motif_dt2, motif_len }
    }

    /// Build the MotifConfig for GCG/C (HhaI) + T/TAA (MseI) — the golden test case.
    fn hhai_msei_config() -> MotifConfig {
        let motif_dt1 = iupac_motifs(&["GCG/C".to_string()]);
        let motif_dt2 = iupac_motifs(&["T/TAA".to_string()]);
        let mut motif_dt = motif_dt1.clone();
        for (k, v) in &motif_dt2 {
            motif_dt.insert(k.clone(), *v);
        }
        let motif_len: HashMap<_, _> = motif_dt
            .keys()
            .map(|p| (p.clone(), regex_bp_len(p)))
            .collect();
        MotifConfig { motif_dt, motif_dt1, motif_dt2, motif_len }
    }

    // --- overlapping_positions ---

    #[test]
    fn test_overlapping_single() {
        // mirrors Python test_dg_digest_seq_1 sequence: GCGC at positions 10 and 24
        let re = Regex::new("GCGC").unwrap();
        let positions = overlapping_positions(&re, "AAAAAAAAAAGCGCAAAAAAAAAAGCGCAAAAAAAAAA");
        assert!(positions.contains(&10), "expected hit at 10, got {:?}", positions);
        assert!(positions.contains(&24), "expected hit at 24, got {:?}", positions);
    }

    #[test]
    fn test_overlapping_nested() {
        // GCGCGC has two overlapping GCGC hits at 0 and 2
        let re = Regex::new("GCGC").unwrap();
        let positions = overlapping_positions(&re, "GCGCGC");
        assert_eq!(positions, vec![0, 2]);
    }

    // --- digest_seq (mirrors test_dg_digest_seq_* in Python) ---

    fn simple_config(motif_strs: &[&str]) -> (MotifConfig, HashMap<String, Regex>) {
        let mut motif_dt = IndexMap::new();
        let mut motif_len = HashMap::new();
        for m in motif_strs {
            let (regex, cut) = crate::iupac::motif_to_regex(m);
            let bp = regex_bp_len(&regex);
            motif_dt.insert(regex.clone(), cut);
            motif_len.insert(regex, bp);
        }
        let regexes = compile_regexes(&motif_len);
        let config = MotifConfig {
            motif_dt1: motif_dt.clone(),
            motif_dt2: motif_dt.clone(),
            motif_dt,
            motif_len,
        };
        (config, regexes)
    }

    #[test]
    fn test_digest_seq_single_motif_pair() {
        // Python test_dg_digest_seq_1
        let seq = "AAAAAAAAAAGCGCAAAAAAAAAAGCGCAAAAAAAAAA";
        let (config, regexes) = simple_config(&["GCGC/"]);
        // Use raw motif to build config with correct cut position
        let mut motif_dt = IndexMap::new();
        motif_dt.insert("GCGC".to_string(), 4usize); // cut after all 4 bases
        let mut motif_len = HashMap::new();
        motif_len.insert("GCGC".to_string(), 4usize);
        let config = MotifConfig {
            motif_dt1: motif_dt.clone(),
            motif_dt2: motif_dt.clone(),
            motif_dt,
            motif_len,
        };
        let regexes = compile_regexes(&config.motif_len);
        let frags = digest_seq("chromname", seq, &config, &regexes, 200);
        // Should find one fragment: GCGCAAAAAAAAAAGCGC at 10..24
        assert_eq!(frags.len(), 1);
        let f = &frags[0];
        assert_eq!(f.seq, "GCGCAAAAAAAAAAGCGC");
        assert_eq!(f.start, 10);
        assert_eq!(f.end, 24);
        assert_eq!(f.internal, 0);
    }

    #[test]
    fn test_digest_seq_two_motifs_max_limit() {
        // Python test_dg_digest_seq_3: max=20 limits which fragments are found
        let seq = "AAAAAAAAAAGCGCAAAAAAAAAAGCGCAAAAAAAAAAGCGCAAAAAAAAAA";
        let mut motif_dt = IndexMap::new();
        motif_dt.insert("GCGC".to_string(), 4usize);
        let mut motif_len = HashMap::new();
        motif_len.insert("GCGC".to_string(), 4usize);
        let config = MotifConfig {
            motif_dt1: motif_dt.clone(),
            motif_dt2: motif_dt.clone(),
            motif_dt,
            motif_len,
        };
        let regexes = compile_regexes(&config.motif_len);
        let frags = digest_seq("chromname", seq, &config, &regexes, 20);
        assert_eq!(frags.len(), 2);
        let seqs: Vec<&str> = frags.iter().map(|f| f.seq.as_str()).collect();
        assert!(seqs.contains(&"GCGCAAAAAAAAAAGCGC"));
    }

    // --- trimmed length (golden test case verification) ---

    #[test]
    fn test_trimmed_length_reverse_fragment() {
        // Golden row 0: GCGCTACCGTGATTAA, m1=TTAA, m2=GCGC, reverse=1 → length=10
        let config = hhai_msei_config();
        let seq = "GCGCTACCGTGATTAA";
        // For reverse=1, m1=TTAA comes from motif_dt2, m2=GCGC from motif_dt1
        // But in process_fragments, frag.m1/m2 are the ORIGINAL (pre-swap) labels.
        // After swap: seq becomes old-revc, which starts with GCGC for this fragment.
        // m1 and m2 labels are preserved from the original digest.
        let length = compute_trimmed_length(seq, "TTAA", "GCGC", 1, &config);
        assert_eq!(length, 10);
    }

    #[test]
    fn test_trimmed_length_forward_fragment() {
        // Golden row 1: GCGCTTTTGCACGTTTAA, m1=GCGC, m2=TTAA, reverse=0 → length=12
        let config = hhai_msei_config();
        let seq = "GCGCTTTTGCACGTTTAA";
        let length = compute_trimmed_length(seq, "GCGC", "TTAA", 0, &config);
        assert_eq!(length, 12);
    }
}
