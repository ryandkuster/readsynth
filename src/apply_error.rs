use anyhow::{Context, Result};
use rand::Rng;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// ---------------------------------------------------------------------------
// Phred probability table
// ---------------------------------------------------------------------------

/// Build a lookup table mapping Phred index (0-41) → error probability.
///
/// Phred quality char c has index `(c as u8 - 33)`.
/// Error probability = 10^(-index/10).
///
/// Mirrors the C++ `scores_arr` initialisation in apply_error.cpp:
///   `scores_arr[i] = pow(10, -(double)i / 10)`
fn build_prob_table() -> [f64; 42] {
    let mut table = [0.0f64; 42];
    for i in 0..42 {
        table[i] = 10_f64.powf(-(i as f64) / 10.0);
    }
    table
}

// ---------------------------------------------------------------------------
// Per-base mutation helpers
// ---------------------------------------------------------------------------

/// Choose uniformly at random from the three bases that are not `base`.
///
/// Mirrors C++'s `mutate_base`: build "ACGT", remove current base,
/// pick randomly from the remaining three.
fn mutate_base(base: char, rng: &mut impl Rng) -> char {
    let others: [char; 3] = match base {
        'A' => ['C', 'G', 'T'],
        'C' => ['A', 'G', 'T'],
        'G' => ['A', 'C', 'T'],
        'T' => ['A', 'C', 'G'],
        _   => ['A', 'C', 'G'], // fallback for non-ACGT (shouldn't reach here)
    };
    others[rng.gen_range(0..3)]
}

/// Apply per-base Phred error to `seq` using the corresponding `score` string.
///
/// For each position:
///   - 'N' bases are passed through unchanged.
///   - Other bases are mutated with probability 10^(-phred/10), where phred
///     is the 0-based index of the quality character in the Phred ASCII range
///     (i.e. `q_char as u8 - 33`).
///
/// Mirrors C++'s `mutate_seq`.
fn mutate_seq(seq: &str, score: &str, prob_table: &[f64; 42], rng: &mut impl Rng) -> String {
    seq.chars()
        .zip(score.chars())
        .map(|(base, q_char)| {
            if base == 'N' {
                return base;
            }
            let idx = (q_char as u8).saturating_sub(33) as usize;
            let prob = prob_table[idx.min(41)];
            // Equivalent to C++'s: rand() < prob * (RAND_MAX + 1.0)
            if rng.gen::<f64>() < prob {
                mutate_base(base, rng)
            } else {
                base
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Read a fastq file, apply per-base Phred error simulation, and write the
/// mutated reads to a new fastq file.
///
/// Quality scores are kept unchanged — only the sequence line is mutated.
///
/// Mirrors the behaviour of `apply_error.cpp`:
///   apply_error <input_fastq> <output_fastq>
pub fn apply_error(input: &Path, output: &Path) -> Result<()> {
    let prob_table = build_prob_table();
    let mut rng = rand::thread_rng();

    let reader = BufReader::new(
        File::open(input).with_context(|| format!("cannot open {}", input.display()))?,
    );
    let mut writer = BufWriter::new(
        File::create(output).with_context(|| format!("cannot create {}", output.display()))?,
    );

    let mut lines = reader.lines();
    loop {
        let header = match lines.next() {
            Some(l) => l.context("error reading fastq header")?,
            None    => break,
        };
        let seq   = lines.next().context("truncated fastq: missing seq")?
                         .context("error reading seq")?;
        let plus  = lines.next().context("truncated fastq: missing '+' line")?
                         .context("error reading '+' line")?;
        let score = lines.next().context("truncated fastq: missing score")?
                         .context("error reading score")?;

        let mutated = mutate_seq(&seq, &score, &prob_table, &mut rng);

        writeln!(writer, "{}", header)?;
        writeln!(writer, "{}", mutated)?;
        writeln!(writer, "{}", plus)?;
        writeln!(writer, "{}", score)?;
    }

    writer.flush()?;
    Ok(())
}
