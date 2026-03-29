use anyhow::{bail, Context, Result};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use readsynth::apply_error::apply_error;
use readsynth::digest::{digest_fasta, digest_fasta_iso, write_raw_digest};
use readsynth::iupac::{
    check_for_palindrome, iupac_motifs, iupac_motifs_iso, regex_bp_len, resolve_enzyme_names,
};
use readsynth::prob_copies::{
    apply_probabilities, apply_probabilities_iso, prob_to_counts, write_counts,
};
use readsynth::size_select::{load_custom_dist, size_select_custom, size_select_gaussian, size_select_low_pass};
use readsynth::types::{CountFragment, MotifConfig, MotifConfigIso};
use readsynth::write_reads::{load_adapters, load_q_scores, write_reads, ReadConfig};

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(name = "readsynth", about = "Simulate Illumina metagenome reads from restriction enzyme digests")]
struct Cli {
    /// Path to genome abundance CSV (genome_path,abundance per line, no header)
    #[arg(short = 'g')]
    genome_file: PathBuf,

    /// Output directory (must already exist)
    #[arg(short = 'o')]
    output_dir: PathBuf,

    /// m1 motif(s): cut notation (e.g. G/AATTC) or RE name
    #[arg(long = "m1", num_args = 1..)]
    m1: Option<Vec<String>>,

    /// m2 motif(s): cut notation or RE name
    #[arg(long = "m2", num_args = 1..)]
    m2: Option<Vec<String>>,

    /// R1 read length
    #[arg(short = 'l')]
    read_len: usize,

    /// Total number of reads to simulate
    #[arg(short = 'n')]
    n_reads: usize,

    /// Mean fragment length for Gaussian size selection
    #[arg(short = 'u')]
    mean: Option<usize>,

    /// Standard deviation for Gaussian size selection
    #[arg(long = "sd")]
    sd: Option<usize>,

    /// Fragment length where the fragment distribution intersects the size distribution
    /// [default: same as -u]
    #[arg(short = 'x')]
    x: Option<usize>,

    /// Low-pass mode: max fragment length, no Gaussian size selection
    #[arg(long = "lp")]
    low_pass: Option<usize>,

    /// Per-site cut efficiency [default: 1.0]
    #[arg(short = 'c', default_value = "1.0")]
    cut_efficiency: f64,

    /// R2 read length [default: same as -l]
    #[arg(long = "l2")]
    l2: Option<usize>,

    /// Number of threads [default: 1]
    #[arg(short = 't', default_value = "1")]
    threads: usize,

    /// Iso-enzyme mode: single Type IIB enzyme name or cut notation
    /// (e.g. alfi or "NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN/")
    #[arg(long = "iso")]
    iso: Option<String>,

    /// Custom R1 adapter file (whitespace-separated: top_seq bot_seq id, one per line)
    #[arg(long = "a1")]
    a1_file: Option<PathBuf>,

    /// Custom R2 adapter file
    #[arg(long = "a2")]
    a2_file: Option<PathBuf>,

    /// R1 adapter sequencing start offset [default: adapter length]
    #[arg(long = "a1s")]
    a1s_override: Option<usize>,

    /// R2 adapter sequencing start offset [default: adapter length]
    #[arg(long = "a2s")]
    a2s_override: Option<usize>,

    /// R1 quality score file (one Phred quality string per line)
    #[arg(long = "q1")]
    q1: Option<PathBuf>,

    /// R2 quality score file
    #[arg(long = "q2")]
    q2: Option<PathBuf>,

    /// Custom fragment length distribution JSON file (length:count pairs).
    /// Overrides Gaussian size selection; max_len is derived from the JSON.
    #[arg(short = 'd')]
    custom_dist: Option<PathBuf>,

    /// Filler base used when fragment + adapter is shorter than read length [default: G]
    #[arg(short = 'e', default_value = "G")]
    filler_base: char,

    /// Test mode: skip fastq output, exit after writing count CSVs
    #[arg(long = "test")]
    test_mode: bool,
}

// ---------------------------------------------------------------------------
// Genome CSV loading
// ---------------------------------------------------------------------------

struct GenomeEntry {
    path: PathBuf,
    name: String,
    abundance: f64,
}

/// Load genome CSV: path,abundance per line (no header).
/// Paths relative to the CSV file's directory are resolved.
/// Abundances are normalised to sum to 1.0.
fn load_genomes(csv_path: &Path) -> Result<Vec<GenomeEntry>> {
    let csv_dir = csv_path.parent().unwrap_or(Path::new("."));
    let content = fs::read_to_string(csv_path)
        .with_context(|| format!("cannot open genome file {}", csv_path.display()))?;

    let mut entries: Vec<(PathBuf, String, f64)> = Vec::new();
    for (lineno, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, ',');
        let raw_path = parts.next().unwrap().trim();
        let abundance_str = parts
            .next()
            .with_context(|| format!("line {}: expected 'path,abundance'", lineno + 1))?
            .trim();
        let abundance: f64 = abundance_str
            .parse()
            .with_context(|| format!("line {}: cannot parse abundance '{}'", lineno + 1, abundance_str))?;

        let path = if Path::new(raw_path).is_absolute() {
            PathBuf::from(raw_path)
        } else {
            csv_dir.join(raw_path)
        };
        let name = path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();

        if !path.exists() {
            bail!("genome not found: {}", path.display());
        }
        entries.push((path, name, abundance));
    }

    if entries.is_empty() {
        bail!("no genomes found in {}", csv_path.display());
    }

    let total: f64 = entries.iter().map(|(_, _, a)| a).sum();
    let genomes = entries
        .into_iter()
        .map(|(path, name, a)| GenomeEntry { path, name, abundance: a / total })
        .collect();

    Ok(genomes)
}

// ---------------------------------------------------------------------------
// MotifConfig builder
// ---------------------------------------------------------------------------

fn build_motif_config(m1_raw: &[String], m2_raw: &[String]) -> Result<MotifConfig> {
    let m1_resolved = resolve_enzyme_names(m1_raw)?;
    let m2_resolved = resolve_enzyme_names(m2_raw)?;
    let m1 = check_for_palindrome(&m1_resolved);
    let m2 = check_for_palindrome(&m2_resolved);

    let motif_dt1 = iupac_motifs(&m1);
    let motif_dt2 = iupac_motifs(&m2);

    if motif_dt1.is_empty() {
        bail!("--m1 produced no valid motifs");
    }
    if motif_dt2.is_empty() {
        bail!("--m2 produced no valid motifs");
    }

    let mut motif_dt = motif_dt1.clone();
    for (k, v) in &motif_dt2 {
        motif_dt.insert(k.clone(), *v);
    }

    let motif_len: HashMap<String, usize> =
        motif_dt.keys().map(|p| (p.clone(), regex_bp_len(p))).collect();

    Ok(MotifConfig { motif_dt, motif_dt1, motif_dt2, motif_len })
}

// ---------------------------------------------------------------------------
// Combined CSV writer
// ---------------------------------------------------------------------------

/// Write combined.csv, mirroring Python's `total_freqs.to_csv(comb_file)`.
/// Columns: (index), length, sum_prob, name, counts_file, counts
fn write_combined_csv(
    path: &Path,
    genome_runs: &[(&str, &Path, &HashMap<usize, f64>)],
    fragment_comps: &HashMap<usize, f64>,
    adjustment: f64,
    max_len: usize,
) -> Result<()> {
    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record(["", "length", "sum_prob", "name", "counts_file", "counts"])?;

    let mut row_idx: usize = 0;
    for (name, counts_file, len_freqs) in genome_runs {
        for length in 0..=max_len {
            let sum_prob = *len_freqs.get(&length).unwrap_or(&0.0);
            let fc = *fragment_comps.get(&length).unwrap_or(&0.0);
            let counts = (fc * sum_prob * adjustment).round();
            wtr.write_record([
                row_idx.to_string(),
                length.to_string(),
                sum_prob.to_string(),
                name.to_string(),
                counts_file.to_string_lossy().into_owned(),
                counts.to_string(),
            ])?;
            row_idx += 1;
        }
    }
    wtr.flush()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Summary CSV writer
// ---------------------------------------------------------------------------

fn write_summary_csv(
    path: &Path,
    genome_runs: &[(String, f64, usize, usize)], // (name, abundance, reads, sites)
) -> Result<()> {
    let total_reads: usize = genome_runs.iter().map(|(_, _, r, _)| r).sum();
    let avg_depths: Vec<f64> = genome_runs
        .iter()
        .map(|(_, _, r, s)| if *s > 0 { *r as f64 / *s as f64 } else { 0.0 })
        .collect();
    let total_avg_depth: f64 = avg_depths.iter().sum();

    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record(["", "genome", "abundance", "reads", "sites",
                       "avg_depth", "depth_abundance", "read_abundance"])?;
    for (i, ((name, abundance, reads, sites), avg_depth)) in
        genome_runs.iter().zip(avg_depths.iter()).enumerate()
    {
        let depth_abundance = if total_avg_depth > 0.0 { avg_depth / total_avg_depth } else { 0.0 };
        let read_abundance  = if total_reads > 0 { *reads as f64 / total_reads as f64 } else { 0.0 };
        wtr.write_record([
            i.to_string(),
            name.clone(),
            abundance.to_string(),
            reads.to_string(),
            sites.to_string(),
            avg_depth.to_string(),
            depth_abundance.to_string(),
            read_abundance.to_string(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

const LOGO: &str = r#"
                                       _/                                  _/      _/
  _/  _/_/    _/_/      _/_/_/    _/_/_/    _/_/_/  _/    _/  _/_/_/    _/_/_/_/  _/_/_/
 _/_/      _/_/_/_/  _/    _/  _/    _/  _/_/      _/    _/  _/    _/    _/      _/    _/
_/        _/        _/    _/  _/    _/      _/_/  _/    _/  _/    _/    _/      _/    _/
_/          _/_/_/    _/_/_/    _/_/_/  _/_/_/      _/_/_/  _/    _/      _/_/  _/    _/
                                                      _/
                                                 _/_/                                      "#;

fn main() -> Result<()> {
    let cli = Cli::parse();

    // -- Print logo ----------------------------------------------------------

    // ANSI bright green
    println!("\x1b[92m{LOGO}\x1b[0m\n");

    // -- Configure thread pool -----------------------------------------------

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .context("failed to configure thread pool")?;

    // -- Validate args -------------------------------------------------------

    let m1_raw = cli.m1.as_deref().unwrap_or(&[]);
    let m2_raw = cli.m2.as_deref().unwrap_or(&[]);
    let is_iso = cli.iso.is_some();

    if !is_iso && (m1_raw.is_empty() || m2_raw.is_empty()) {
        bail!("must specify either --iso <enzyme> or both --m1 and --m2");
    }

    // Load custom distribution JSON early so max_len is known before the digest.
    // Mirrors Python's check_custom_distribution() call at startup.
    let custom_dist_data: Option<(std::collections::HashMap<usize, f64>, usize)> =
        if let Some(ref path) = cli.custom_dist {
            Some(load_custom_dist(path).context("failed to load custom distribution")?)
        } else {
            None
        };

    let max_len = if let Some((_, ml)) = &custom_dist_data {
        *ml
    } else if let Some(lp) = cli.low_pass {
        lp
    } else if let (Some(u), Some(sd)) = (cli.mean, cli.sd) {
        u + 6 * sd
    } else {
        bail!("must specify either (-u and --sd), --lp, or -d");
    };

    let x = cli.x.or(cli.mean).unwrap_or(max_len / 2);
    let l1 = cli.read_len;
    let l2 = cli.l2.unwrap_or(l1);
    let c  = cli.cut_efficiency;
    let n  = cli.n_reads;

    // -- Create output subdirectories ----------------------------------------

    let out = &cli.output_dir;
    let raw_digests_dir = out.join("raw_digests");
    let individual_dir  = out.join("individual_counts");

    fs::create_dir_all(&raw_digests_dir)
        .with_context(|| format!("cannot create {}", raw_digests_dir.display()))?;
    fs::create_dir_all(&individual_dir)
        .with_context(|| format!("cannot create {}", individual_dir.display()))?;

    // -- Build motif config --------------------------------------------------

    // Standard mode
    let motif_config = if !is_iso {
        Some(build_motif_config(m1_raw, m2_raw)?)
    } else {
        None
    };

    // Iso mode
    let motif_config_iso: Option<MotifConfigIso> = if is_iso {
        let raw = cli.iso.as_deref().unwrap();
        let motif_str = if raw.contains('/') {
            raw.to_ascii_uppercase()
        } else {
            readsynth::iupac::lookup_enzyme(raw)
                .ok_or_else(|| anyhow::anyhow!(
                    "unknown iso enzyme '{}' — use cut notation or check spelling", raw
                ))?
                .to_string()
        };
        if motif_str.chars().filter(|&c| c == '/').count() != 2 {
            bail!("--iso motif must have exactly two '/' cut markers: {}", motif_str);
        }
        let motif_dt = iupac_motifs_iso(&motif_str);
        let motif_len: HashMap<String, usize> =
            motif_dt.keys().map(|p| (p.clone(), regex_bp_len(p))).collect();
        Some(MotifConfigIso { motif_dt, motif_len })
    } else {
        None
    };

    // -- Load genomes --------------------------------------------------------

    let genomes = load_genomes(&cli.genome_file)?;
    println!("Found {} genome(s)", genomes.len());

    // -- Stage 1: digest + prob_copies per genome ----------------------------

    println!("\n1. Digesting genomes and computing fragment probabilities...");

    let pb = ProgressBar::new(genomes.len() as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "   [{bar:45.green/white}] {pos}/{len} genomes  ({elapsed})",
        )
        .unwrap()
        .progress_chars("▰▱▱"),
    );

    struct PerGenome {
        name: String,
        abundance: f64,
        count_file: PathBuf,
        len_freqs: HashMap<usize, f64>,
        frags: Vec<CountFragment>,
    }

    // Stage 1 runs in parallel across genomes — each genome is fully independent.
    let stage1: Vec<Result<Option<PerGenome>>> = genomes
        .par_iter()
        .map(|genome| -> Result<Option<PerGenome>> {
            let processed = if let Some(ref cfg) = motif_config {
                digest_fasta(&genome.path, cfg, max_len)
            } else {
                digest_fasta_iso(&genome.path, motif_config_iso.as_ref().unwrap(), max_len)
            }
            .with_context(|| format!("digest failed for {}", genome.path.display()))?;

            if processed.is_empty() {
                pb.inc(1);
                return Ok(None);
            }

            let raw_path = raw_digests_dir.join(format!("raw_digest_{}.csv", genome.name));
            write_raw_digest(&processed, &raw_path)
                .with_context(|| format!("cannot write raw digest for {}", genome.name))?;

            let (count_frags, len_freqs) = if is_iso {
                apply_probabilities_iso(processed, genome.abundance, max_len)
            } else {
                apply_probabilities(processed, c, genome.abundance, max_len)
            };

            let count_path = individual_dir.join(format!("counts_{}.csv", genome.name));
            write_counts(&count_frags, &count_path)
                .with_context(|| format!("cannot write counts for {}", genome.name))?;

            pb.inc(1);
            Ok(Some(PerGenome {
                name: genome.name.clone(),
                abundance: genome.abundance,
                count_file: count_path,
                len_freqs,
                frags: count_frags,
            }))
        })
        .collect();

    pb.finish_and_clear();

    // Propagate any per-genome errors, then discard empty genomes.
    let per_genome: Vec<PerGenome> = stage1
        .into_iter()
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();

    // Accumulate total_len_freqs sequentially after all genomes are done.
    let mut total_len_freqs: HashMap<usize, f64> = (0..=max_len).map(|i| (i, 0.0)).collect();
    for g in &per_genome {
        for (length, prob) in &g.len_freqs {
            *total_len_freqs.get_mut(length).unwrap() += prob;
        }
    }

    if per_genome.is_empty() {
        bail!("no fragments produced across all genomes — check motifs and genome files");
    }

    // -- Stage 2: size selection ---------------------------------------------

    println!("\n2. Simulating size selection...");

    let ss = if let Some((ref dist, dist_max_len)) = custom_dist_data {
        size_select_custom(&total_len_freqs, n, dist, dist_max_len)
            .context("custom distribution size selection failed")?
    } else if cli.low_pass.is_some() || is_iso {
        size_select_low_pass(&total_len_freqs, n, max_len)
    } else {
        let u  = cli.mean.unwrap();
        let sd = cli.sd.unwrap();
        size_select_gaussian(&total_len_freqs, n, u, sd, x, max_len)
            .context("Gaussian size selection failed")?
    };

    // -- Stage 3: write combined.csv -----------------------------------------

    let combined_path = out.join("combined.csv");
    let genome_run_refs: Vec<(&str, &Path, &HashMap<usize, f64>)> = per_genome
        .iter()
        .map(|g| (g.name.as_str(), g.count_file.as_path(), &g.len_freqs))
        .collect();
    write_combined_csv(&combined_path, &genome_run_refs, &ss.fragment_comps, ss.adjustment, max_len)
        .context("cannot write combined.csv")?;

    // -- Stage 4: assign read counts and rewrite count CSVs ------------------

    // -- Build ReadConfig (adapters + quality scores) -------------------------

    let mut rc = if let (Some(a1_path), Some(a2_path)) = (&cli.a1_file, &cli.a2_file) {
        // Custom adapter files supplied
        let adapters1 = load_adapters(a1_path)
            .with_context(|| format!("cannot load R1 adapters from {}", a1_path.display()))?;
        let adapters2 = load_adapters(a2_path)
            .with_context(|| format!("cannot load R2 adapters from {}", a2_path.display()))?;
        let a1s = cli.a1s_override.unwrap_or(adapters1[0].0.len());
        let a2s = cli.a2s_override.unwrap_or(adapters2[0].0.len());
        ReadConfig {
            adapters1, adapters2, a1s, a2s, l1, l2,
            filler_base: cli.filler_base, q1_scores: None, q2_scores: None,
        }
    } else if is_iso {
        // Iso-enzyme default adapters with N overhangs
        let iso_cfg = motif_config_iso.as_ref().unwrap();
        let keys: Vec<&String> = iso_cfg.motif_dt.keys().collect();
        let m1_cut1 = iso_cfg.motif_dt[keys[0]][0];
        let m2_cut1 = iso_cfg.motif_dt[keys[1]][0];
        ReadConfig::iso_adapters(l1, Some(l2), cli.filler_base, m1_cut1, m2_cut1)
    } else {
        ReadConfig::default_adapters(l1, Some(l2), cli.filler_base)
    };

    if let Some(q1_path) = &cli.q1 {
        rc.q1_scores = Some(load_q_scores(q1_path)
            .with_context(|| format!("cannot load R1 q-scores from {}", q1_path.display()))?);
    }
    if let Some(q2_path) = &cli.q2 {
        rc.q2_scores = Some(load_q_scores(q2_path)
            .with_context(|| format!("cannot load R2 q-scores from {}", q2_path.display()))?);
    }
    let mut summary_rows: Vec<(String, f64, usize, usize)> = Vec::new();

    // Stage 5 (fastq) happens here per-genome, or is skipped in test mode
    let (r1_path, r2_path) = if !cli.test_mode {
        (Some(out.join("sim_metagenome_R1.fastq")),
         Some(out.join("sim_metagenome_R2.fastq")))
    } else {
        (None, None)
    };

    let mut r1_file = r1_path
        .as_deref()
        .map(|p| fs::File::create(p).with_context(|| format!("cannot create {}", p.display())))
        .transpose()?;
    let mut r2_file = r2_path
        .as_deref()
        .map(|p| fs::File::create(p).with_context(|| format!("cannot create {}", p.display())))
        .transpose()?;

    println!("\n3. Assigning read counts and writing outputs...");

    let pb2 = ProgressBar::new(genomes.len() as u64);
    pb2.set_style(
        ProgressStyle::with_template(
            "   [{bar:45.green/white}] {pos}/{len} genomes  ({elapsed})",
        )
        .unwrap()
        .progress_chars("▰▱▱"),
    );

    // Stage 4a (parallel): assign counts and generate read bytes per genome.
    struct Stage4 {
        name:        String,
        abundance:   f64,
        count_file:  PathBuf,
        final_frags: Vec<CountFragment>,
        r1_buf:      Vec<u8>,
        r2_buf:      Vec<u8>,
    }

    let test_mode = cli.test_mode;
    let stage4: Vec<Result<Stage4>> = per_genome
        .into_par_iter()
        .map(|g| -> Result<Stage4> {
            let final_frags = prob_to_counts(g.frags, &ss.fragment_comps, ss.adjustment);
            let mut r1_buf: Vec<u8> = Vec::new();
            let mut r2_buf: Vec<u8> = Vec::new();
            if !test_mode {
                write_reads(&final_frags, &mut r1_buf, &mut r2_buf, &g.name, &rc)
                    .with_context(|| format!("write_reads failed for {}", g.name))?;
            }
            pb2.inc(1);
            Ok(Stage4 {
                name:        g.name,
                abundance:   g.abundance,
                count_file:  g.count_file,
                final_frags,
                r1_buf,
                r2_buf,
            })
        })
        .collect();

    pb2.finish_and_clear();

    // Stage 4b (sequential): flush CSVs and fastq files in genome order.
    for res in stage4.into_iter().collect::<Result<Vec<_>>>()? {
        let total_reads: u64 = res.final_frags.iter().map(|f| f.counts).sum();
        let sites = res.final_frags.len();

        write_counts(&res.final_frags, &res.count_file)
            .with_context(|| format!("cannot rewrite counts for {}", res.name))?;

        summary_rows.push((res.name, res.abundance, total_reads as usize, sites));

        if let (Some(r1), Some(r2)) = (r1_file.as_mut(), r2_file.as_mut()) {
            r1.write_all(&res.r1_buf)?;
            r2.write_all(&res.r2_buf)?;
        }
    }

    // -- Stage 5: write summary ----------------------------------------------

    let summary_path = out.join("metagenome_summary.csv");
    write_summary_csv(&summary_path, &summary_rows).context("cannot write metagenome_summary.csv")?;

    // -- Stage 5: apply_error (only when --q1 and --q2 are both supplied) ----

    let (error_r1_path, error_r2_path) =
        if cli.q1.is_some() && cli.q2.is_some() && !cli.test_mode {
            let er1 = out.join("error_sim_metagenome_R1.fastq");
            let er2 = out.join("error_sim_metagenome_R2.fastq");
            let sim_r1 = r1_path.as_deref().unwrap();
            let sim_r2 = r2_path.as_deref().unwrap();

            println!("\n4. Applying Phred error simulation...");
            apply_error(sim_r1, &er1)
                .context("apply_error failed for R1")?;
            apply_error(sim_r2, &er2)
                .context("apply_error failed for R2")?;

            (Some(er1), Some(er2))
        } else {
            (None, None)
        };

    // -- Report --------------------------------------------------------------

    let total_reads: usize = summary_rows.iter().map(|(_, _, r, _)| r).sum();
    println!("\nDone. Total reads: {}", total_reads);
    println!("  Combined:  {}", combined_path.display());
    println!("  Summary:   {}", summary_path.display());
    if let Some(p) = &r1_path {
        println!("  R1 fastq:  {}", p.display());
        println!("  R2 fastq:  {}", r2_path.as_deref().unwrap().display());
    }
    if let Some(p) = &error_r1_path {
        println!("  Error R1:  {}", p.display());
        println!("  Error R2:  {}", error_r2_path.as_deref().unwrap().display());
    }

    Ok(())
}
