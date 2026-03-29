use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::Path;

/// Output of size selection: per-length composition fraction and a global read-count multiplier.
/// Mirrors Python's size_selection.main() return values (fragment_comps, adjustment).
pub struct SizeSelectResult {
    /// For each fragment length 0..=max_len: the fraction of available fragments to include.
    /// Value is in [0, 1] for Gaussian mode, or 0/1 for low-pass mode.
    pub fragment_comps: HashMap<usize, f64>,
    /// Scalar: multiply fragment_comps × adj_prob × adjustment to get read counts.
    pub adjustment: f64,
}

/// Gaussian probability density function at point `x`.
///
/// IMPORTANT: mirrors Python's operator-precedence quirk:
///   `(1 / args.sd * np.sqrt(2 * np.pi)) * np.exp(...)`
/// which evaluates as `((1/sd) * sqrt(2π)) * exp(...)`,
/// NOT the standard `(1 / (sd * sqrt(2π))) * exp(...)`.
/// The scale_by computation cancels the constant, but we must match exactly.
fn gauss_pdf(x: usize, mean: f64, sd: f64) -> f64 {
    let xf = x as f64;
    (1.0 / sd * (2.0 * std::f64::consts::PI).sqrt())
        * ((-0.5) * ((xf - mean) / sd).powi(2)).exp()
}

/// Size-select fragments using a Normal distribution.
/// Mirrors Python's size_selection.main() for the Gaussian case.
///
/// Parameters:
/// - `len_dt`:  length → sum(adj_prob) across all genomes (from apply_probabilities)
/// - `n`:       total desired read count (-n flag)
/// - `mean`:    Gaussian mean (-u flag)
/// - `sd`:      Gaussian standard deviation (-sd flag)
/// - `x`:       intersection point (-x flag); when absent, defaults to `mean`
/// - `max_len`: maximum fragment length (u + 6*sd)
pub fn size_select_gaussian(
    len_dt: &HashMap<usize, f64>,
    n: usize,
    mean: usize,
    sd: usize,
    x: usize,
    max_len: usize,
) -> Result<SizeSelectResult> {
    let mean_f = mean as f64;
    let sd_f   = sd as f64;

    // x_dist = min(x, abs(max_len - sd), max_len - x)
    // Note: Python uses args.sd (not args.u) in the abs() term — this is intentional.
    let x_dist = x
        .min((max_len as isize - sd as isize).unsigned_abs())
        .min(max_len - x);

    // Average fragment density in range [x-x_dist, x+x_dist]
    let x_range: Vec<f64> = ((x - x_dist)..=(x + x_dist))
        .map(|i| *len_dt.get(&i).unwrap_or(&0.0))
        .collect();

    if x_range.is_empty() {
        anyhow::bail!("no fragment lengths in x_range around mean={mean}");
    }
    let avg_x = x_range.iter().sum::<f64>() / x_range.len() as f64;

    let scale_by = avg_x / gauss_pdf(x, mean_f, sd_f);

    // draw_dt[i] = min(gauss_pdf(i)*scale_by, len_dt[i]) * 0.65  (Sage average recovery)
    let mut draw_dt: HashMap<usize, f64> = HashMap::new();
    for i in 0..=max_len {
        let draw_counts = gauss_pdf(i, mean_f, sd_f) * scale_by;
        let data_counts = *len_dt.get(&i).unwrap_or(&0.0);
        draw_dt.insert(i, draw_counts.min(data_counts) * 0.65);
    }

    let total_draw: f64 = draw_dt.values().sum();
    if total_draw == 0.0 {
        anyhow::bail!("no fragments survive Gaussian size selection");
    }

    let adjustment = n as f64 / total_draw;

    // fragment_comps[i] = draw_dt[i] / len_dt[i], or 0 when either is 0 / NaN
    let mut fragment_comps: HashMap<usize, f64> = HashMap::new();
    for i in 0..=max_len {
        let ld = *len_dt.get(&i).unwrap_or(&0.0);
        let dd = *draw_dt.get(&i).unwrap_or(&0.0);
        let comp = if ld == 0.0 || dd == 0.0 {
            0.0
        } else {
            let v = dd / ld;
            if v.is_nan() { 0.0 } else { v }
        };
        fragment_comps.insert(i, comp);
    }

    Ok(SizeSelectResult { fragment_comps, adjustment })
}

/// Low-pass size selection: include all fragments up to max_len uniformly.
/// Mirrors Python's low-pass branch in readsynth.py (iso or -lp mode).
pub fn size_select_low_pass(
    len_dt: &HashMap<usize, f64>,
    n: usize,
    max_len: usize,
) -> SizeSelectResult {
    // fragment_comps: 1 where there are fragments, 0 otherwise (all-or-nothing)
    let fragment_comps: HashMap<usize, f64> = (0..=max_len)
        .map(|i| {
            let prob = *len_dt.get(&i).unwrap_or(&0.0);
            (i, if prob > 0.0 { 1.0 } else { 0.0 })
        })
        .collect();

    // Python: adjustment = n / sum(fragment_comps.values()) where fragment_comps still holds
    // the sum of adj_probs per length (before the 0/1 conversion). That is len_dt.
    let total: f64 = len_dt.values().sum();
    let adjustment = if total > 0.0 { n as f64 / total } else { 0.0 };

    SizeSelectResult { fragment_comps, adjustment }
}

/// Load a custom distribution JSON file and return (dist, max_len).
///
/// Call this early in the pipeline so max_len is known before the digest runs.
/// The returned HashMap maps fragment length → expected draw count.
pub fn load_custom_dist(json_path: &Path) -> Result<(HashMap<usize, f64>, usize)> {
    let content = std::fs::read_to_string(json_path)
        .with_context(|| format!("cannot open custom distribution file {}", json_path.display()))?;

    let raw: serde_json::Map<String, serde_json::Value> =
        serde_json::from_str(&content)
            .context("custom distribution JSON must be an object of length:count pairs")?;

    let mut dist: HashMap<usize, f64> = HashMap::new();
    for (k, v) in &raw {
        let length: usize = k.parse()
            .with_context(|| format!("non-integer key in custom distribution JSON: '{k}'"))?;
        let count = v.as_f64()
            .with_context(|| format!("non-numeric value for key '{k}' in custom distribution JSON"))?;
        dist.insert(length, count);
    }

    if dist.is_empty() {
        anyhow::bail!("custom distribution JSON is empty");
    }

    let max_len = *dist.keys().max().unwrap();
    Ok((dist, max_len))
}

/// Custom distribution size selection using a pre-loaded dist map.
///
/// Mirrors Python's `get_custom_dist` + `get_custom_draw_dict` path in
/// size_selection.main(). Draw counts come directly from the JSON;
/// lengths absent from the file get draw count 0.
pub fn size_select_custom(
    len_dt: &HashMap<usize, f64>,
    n: usize,
    dist: &HashMap<usize, f64>,
    max_len: usize,
) -> Result<SizeSelectResult> {
    // draw_dt[i] = dist[i] if present, else 0
    let draw_dt: HashMap<usize, f64> = (0..=max_len)
        .map(|i| (i, *dist.get(&i).unwrap_or(&0.0)))
        .collect();

    let total_draw: f64 = draw_dt.values().sum();
    if total_draw == 0.0 {
        anyhow::bail!("custom distribution has no non-zero counts");
    }

    let adjustment = n as f64 / total_draw;

    let fragment_comps: HashMap<usize, f64> = (0..=max_len)
        .map(|i| {
            let ld = *len_dt.get(&i).unwrap_or(&0.0);
            let dd = *draw_dt.get(&i).unwrap_or(&0.0);
            let comp = if ld == 0.0 || dd == 0.0 {
                0.0
            } else {
                let v = dd / ld;
                if v.is_nan() { 0.0 } else { v }
            };
            (i, comp)
        })
        .collect();

    Ok(SizeSelectResult { fragment_comps, adjustment })
}
