use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::path::Path;

use crate::types::{CountFragment, ProcessedFragment};

/// Compute per-fragment probability given digestion efficiency `c` and
/// the number of internal RE sites within the fragment.
///
/// Probability = c² × (1 - c)^internal_sites
/// Mirrors Python's copies_dict() + apply_approach() in prob_n_copies.py.
pub fn fragment_probability(c: f64, internal: usize) -> f64 {
    c.powi(2) * (1.0 - c).powi(internal as i32)
}

/// For fragments that share the same (start, end) coordinates (bidirectional),
/// halve their probability — they compete for the same ligation event.
/// Mirrors Python's bidirectional_weights().
pub fn apply_bidirectional_weights(fragments: &mut [CountFragment]) {
    // Find (start, end) pairs that appear more than once
    let mut coord_counts: HashMap<(usize, usize), usize> = HashMap::new();
    for f in fragments.iter() {
        *coord_counts.entry((f.fragment.start, f.fragment.end)).or_insert(0) += 1;
    }
    for f in fragments.iter_mut() {
        if coord_counts[&(f.fragment.start, f.fragment.end)] > 1 {
            f.probability /= 2.0;
        }
    }
}

/// Assign probabilities and adjusted probabilities to all fragments for one genome.
/// Returns the per-length summed adj_prob dictionary used by size_selection.
///
/// Mirrors Python's prob_n_copies.main().
pub fn apply_probabilities(
    fragments: Vec<ProcessedFragment>,
    c: f64,
    abundance: f64,
    max_len: usize,
) -> (Vec<CountFragment>, HashMap<usize, f64>) {
    let mut counts: Vec<CountFragment> = fragments
        .into_iter()
        .map(|f| {
            let prob = fragment_probability(c, f.internal);
            CountFragment {
                fragment: f,
                probability: prob,
                adj_prob: 0.0, // filled in below after bidirectional adjustment
                counts: 0,
            }
        })
        .collect();

    // Remove zero-probability fragments
    counts.retain(|f| f.probability > 0.0);

    apply_bidirectional_weights(&mut counts);

    // Adjusted probability = probability × relative abundance
    for f in counts.iter_mut() {
        f.adj_prob = f.probability * abundance;
    }

    // Build length → sum(adj_prob) dictionary for size selection
    let mut len_freqs: HashMap<usize, f64> = (0..=max_len).map(|i| (i, 0.0)).collect();
    for f in &counts {
        if f.fragment.length <= max_len {
            *len_freqs.get_mut(&f.fragment.length).unwrap() += f.adj_prob;
        }
    }

    (counts, len_freqs)
}

/// Assign integer read counts to each fragment from size-selection results,
/// and drop fragments that receive zero reads.
///
/// Mirrors Python's prob_to_counts():
///   counts = round(fragment_comps[length] * adj_prob * adjustment)
///   keep only counts > 0
pub fn prob_to_counts(
    fragments: Vec<CountFragment>,
    fragment_comps: &HashMap<usize, f64>,
    adjustment: f64,
) -> Vec<CountFragment> {
    fragments
        .into_iter()
        .filter_map(|mut f| {
            let comp = *fragment_comps.get(&f.fragment.length).unwrap_or(&0.0);
            let raw = (comp * f.adj_prob * adjustment).round();
            if raw <= 0.0 {
                return None;
            }
            f.counts = raw as u64;
            Some(f)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Iso-enzyme probability calculation
// ---------------------------------------------------------------------------

/// For each fragment find which other fragments overlap it spatially.
/// Mirrors Python's prob_n_copies_iso.find_overlaps().
fn find_overlaps_iso(fragments: &[ProcessedFragment]) -> Vec<Option<Vec<usize>>> {
    let mut overlaps = vec![None; fragments.len()];
    for i in 0..fragments.len() {
        let (si, ei) = (fragments[i].start, fragments[i].end);
        let mut hits: Vec<usize> = (0..fragments.len())
            .filter(|&j| {
                j != i && {
                    let (sj, ej) = (fragments[j].start, fragments[j].end);
                    (sj > si && sj < ei) || (sj < si && ej > si)
                }
            })
            .collect();
        if !hits.is_empty() {
            overlaps[i] = Some(hits);
        }
    }
    overlaps
}

/// Build all non-exclusive fragment combination events for one overlap cluster.
/// Mirrors Python's prob_n_copies_iso.count_events().
fn count_events_iso(linked: &[Vec<usize>]) -> (usize, usize, Vec<Vec<usize>>) {
    let first = linked[0][0].saturating_sub(1);
    let last = first + linked.len() - 1;

    // For each fragment index, list the indices it is *compatible* with (doesn't overlap).
    let mut compat: HashMap<usize, Vec<usize>> = HashMap::new();
    for (idx1, _) in linked.iter().enumerate() {
        let current = first + idx1;
        for (idx2, link) in linked.iter().enumerate() {
            if idx1 != idx2 && !link.contains(&current) {
                compat.entry(current).or_default().push(first + idx2);
            }
        }
    }

    let mut events: Vec<Vec<usize>> = Vec::new();
    for k in first..=last {
        let mut event = vec![k];
        match compat.get(&k) {
            None => {
                events.push(event);
            }
            Some(v) => {
                for &i in v {
                    event.push(i);
                    if let Some(inner) = compat.get(&i) {
                        for &j in inner {
                            if v.contains(&j) {
                                event.push(j);
                            }
                        }
                    }
                    events.push(event.clone());
                    event = vec![k];
                }
            }
        }
    }

    // Deduplicate (sort each event then use a HashSet).
    let unique: HashSet<Vec<usize>> = events
        .into_iter()
        .map(|mut e| { e.sort(); e })
        .collect();
    let events: Vec<Vec<usize>> = unique.into_iter().collect();

    (first, last, events)
}

/// Assign per-fragment probabilities from the event list.
/// Mirrors Python's prob_n_copies_iso.apply_probs().
fn apply_probs_iso(first: usize, last: usize, events: &[Vec<usize>]) -> Vec<f64> {
    let denom = events.len() as f64;
    (first..=last)
        .map(|i| events.iter().filter(|e| e.contains(&i)).count() as f64 / denom)
        .collect()
}

/// Walk through the overlap list, cluster linked fragments, and compute probabilities.
/// Mirrors Python's prob_n_copies_iso.calculate_prob().
fn calculate_prob_iso(overlaps: &[Option<Vec<usize>>]) -> Vec<f64> {
    let mut all_probs: Vec<f64> = Vec::new();
    let mut prev_max: usize = 0;
    let mut linked: Vec<Vec<usize>> = Vec::new();

    for overlap in overlaps {
        match overlap {
            Some(idxs) => {
                linked.push(idxs.clone());
                let curr_max = *idxs.iter().max().unwrap();
                if prev_max > curr_max {
                    let (first, last, events) = count_events_iso(&linked);
                    all_probs.extend(apply_probs_iso(first, last, &events));
                    linked.clear();
                }
                prev_max = curr_max;
            }
            None => {
                all_probs.push(1.0);
            }
        }
    }
    all_probs
}

/// Assign probabilities to iso-enzyme fragments using the overlap graph algorithm.
/// Mirrors Python's prob_n_copies_iso.main().
pub fn apply_probabilities_iso(
    fragments: Vec<ProcessedFragment>,
    abundance: f64,
    max_len: usize,
) -> (Vec<CountFragment>, HashMap<usize, f64>) {
    if fragments.is_empty() {
        return (vec![], (0..=max_len).map(|i| (i, 0.0)).collect());
    }

    let overlaps = find_overlaps_iso(&fragments);
    let probs = calculate_prob_iso(&overlaps);

    // Zip probs with fragments; pad with 1.0 for any cluster not finalized
    // (matches Python's behavior for isolated trailing fragments).
    let mut counts: Vec<CountFragment> = fragments
        .into_iter()
        .zip(probs.into_iter().chain(std::iter::repeat(1.0_f64)))
        .map(|(f, prob)| CountFragment {
            fragment: f,
            probability: prob,
            adj_prob: prob * abundance,
            counts: 0,
        })
        .collect();

    counts.retain(|f| f.probability > 0.0);

    let mut len_freqs: HashMap<usize, f64> = (0..=max_len).map(|i| (i, 0.0)).collect();
    for f in &counts {
        if f.fragment.length <= max_len {
            *len_freqs.get_mut(&f.fragment.length).unwrap() += f.adj_prob;
        }
    }

    (counts, len_freqs)
}

/// Write counts fragments to CSV (counts_*.csv).
///
/// Mirrors Python's `df.to_csv(prob_file, index=None)` — which, because `index=None` is
/// treated as `True` by pandas, writes a leading integer index column with an empty header.
/// Column order: (idx), chrom, seq, start, end, m1, m2, internal, reverse, revc, length,
///               probability, adj_prob, counts
pub fn write_counts(fragments: &[CountFragment], path: &Path) -> Result<()> {
    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record([
        "", "chrom", "seq", "start", "end", "m1", "m2", "internal",
        "reverse", "revc", "length", "probability", "adj_prob", "counts",
    ])?;
    for (i, f) in fragments.iter().enumerate() {
        wtr.write_record([
            i.to_string(),
            f.fragment.chrom.clone(),
            f.fragment.seq.clone(),
            f.fragment.start.to_string(),
            f.fragment.end.to_string(),
            f.fragment.m1.clone(),
            f.fragment.m2.clone(),
            f.fragment.internal.to_string(),
            f.fragment.reverse.to_string(),
            f.fragment.revc.clone(),
            f.fragment.length.to_string(),
            f.probability.to_string(),
            f.adj_prob.to_string(),
            f.counts.to_string(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}
