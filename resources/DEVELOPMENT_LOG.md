# readsynth-rs Development Log

Rust rewrite of `readsynth.py` — a bioinformatics tool that simulates Illumina paired-end
metagenome reads from restriction enzyme (RE) digests of input genomes.

---

## Phase 1: Core digest pipeline

### `iupac.rs`
- Implemented `iupac_char_to_regex` — maps IUPAC ambiguity codes to regex character classes
- Implemented `motif_to_regex(motif: &str) -> (String, usize)` — parses cut notation (e.g. `G/AATTC`)
  into a regex pattern and cut position
- Implemented `iupac_motifs(motifs: &[String]) -> IndexMap<String, usize>` — builds the
  `motif_dt` map used throughout the pipeline
- Implemented `check_for_palindrome` — adds RC motif entry for non-palindromic enzymes
- Implemented `reverse_comp` — full IUPAC-aware reverse complement
- Implemented `regex_bp_len` — computes base-pair length of a regex pattern (handles `[ACGT]`
  bracket groups as single bases)

### `types.rs`
- Defined `RawFragment` — intermediate digest output before processing
- Defined `ProcessedFragment` — post-processing output (trimmed length, reverse flag, revc)
- Defined `CountFragment` — fragment with adj_prob and counts assigned
- Defined `MotifConfig` — holds `motif_dt`, `motif_dt1`, `motif_dt2`, `motif_len`

### `digest.rs`
- Implemented `overlapping_positions` — finds all overlapping RE site positions using
  the lookahead trick from Python (`re.finditer('(?=' + pat + ')', seq)`)
- Implemented `digest_seq` + `digest_frag` — per-chromosome fragment search, iterating
  `motif_dt` in `IndexMap` insertion order to match Python's deterministic output
- Implemented `process_fragments` — filters unviable combos, handles IUPAC→N conversion,
  reverse complement, bidirectional fragment duplication, and trimmed length calculation
- Implemented `compute_trimmed_length` — replicates Python's overhang trimming loop exactly,
  including the `back==0` edge case
- Implemented `digest_fasta` — reads gzip or plain FASTA, runs digest, returns sorted fragments
- Implemented `write_raw_digest` — writes CSV matching Python's column order

### `prob_copies.rs`
- Implemented `fragment_probability(c, internal)` — models per-site cut efficiency
- Implemented `apply_probabilities` — assigns `probability` and `adj_prob` per fragment,
  builds `len_freqs` dictionary
- Implemented `prob_to_counts` — applies `fragment_comps × adj_prob × adjustment`, rounds to integer counts
- Implemented `write_counts` — writes per-genome counts CSV

### `size_select.rs`
- Implemented `gauss_pdf` — replicates Python's non-standard operator precedence:
  `(1/sd * sqrt(2π)) * exp(...)` (not the standard `1/(sd*sqrt(2π))`; the constants cancel
  in `scale_by` but must match exactly for floating-point reproducibility)
- Implemented `size_select_gaussian` — full Gaussian size selection pipeline
  (`get_scale_gauss` → `get_draw_dict` → `fragment_comps` + `adjustment`)
- Implemented `size_select_low_pass` — flat/all-or-nothing selection for ISO and `-lp` modes

### `write_reads.rs`
- Implemented `ReadConfig` — holds adapter sequences, read lengths, filler base, and
  optional quality score pools
- Implemented `ReadConfig::default_adapters` and `ReadConfig::iso_adapters`
- Implemented `write_reads` — generates paired R1/R2 reads with adapter ligation,
  handles insert shorter than read length (adapter contamination with filler fill)

---

## Phase 2: Full pipeline wiring (`main.rs`)

- Defined `Cli` struct with all args: `-g`, `-o`, `--m1`, `--m2`, `-l`, `-n`, `-u`, `--sd`,
  `-x`, `--lp`, `-c`, `--l2`, `-t`, `--iso`, `--a1`, `--a2`, `--a1s`, `--a2s`, `--q1`, `--q2`, `--test`
- Implemented `load_genomes` — reads abundance CSV, normalises abundances to sum 1.0
- Implemented `build_motif_config` — resolves enzyme names, checks palindromes, builds `MotifConfig`
- Implemented `write_combined_csv` — mirrors Python's `total_freqs.to_csv(comb_file)`
- Implemented `write_summary_csv` — mirrors Python's `write_final_file`
- Wired Stage 1 → Stage 2 → Stage 3 → Stage 4 pipeline

---

## Phase 3: Named enzyme lookup

**Problem:** Users had to remember cut notation (e.g. `G/AATTC`); Python supported names like
`ecori`, `msei`.

**Solution:** Added `lookup_enzyme(name: &str) -> Option<&'static str>` in `iupac.rs` — a static
`match` covering 176 Type IIP and 22 Type IIB enzymes. Added `resolve_enzyme_names` to apply
this lookup (or pass cut notation through unchanged) before building `MotifConfig`.

---

## Phase 4: Rayon parallelism

**Design:** Stage 1 (digest + prob_copies per genome) runs with `par_iter()`. Stage 4 (read
generation) generates per-genome `Vec<u8>` buffers in parallel, then flushes them sequentially
to the shared fastq files.

- `total_len_freqs` accumulation is done sequentially after `par_iter` collects, avoiding
  shared mutable state
- `indicatif::ProgressBar` is `Send + Sync` — safe to call `pb.inc(1)` from rayon threads
- Default thread count is 1 (`-t 1`), configurable with `-t N`

---

## Phase 5: Progress bars

- Added `indicatif` dependency
- Stage 1: `[{bar:45.green/white}] {pos}/{len} genomes` with `▰▱▱` progress characters
- Stage 4: same style for read assignment/writing
- Both bars call `finish_and_clear()` on completion to keep terminal clean

---

## Phase 6: ASCII logo

- Added `LOGO` const as a raw string (`r#"..."#`) to avoid escaped-newline warnings
- Printed at startup in ANSI bright green: `\x1b[92m{LOGO}\x1b[0m`

---

## Phase 7: ISO enzyme mode (`--iso`)

### Parsing
- Added `iupac_motifs_iso(motif: &str) -> IndexMap<String, [usize; 2]>` — parses a Type IIB
  motif with two `/` cut markers, adds RC automatically, returns `[cut1, cut2]` per pattern
- `cut2 = slash_positions[1] - 1` (the `-1` mirrors Python's adjustment for the removed first slash)
- Added `MotifConfigIso` struct with `motif_dt: IndexMap<String, [usize; 2]>` and `motif_len`

### Digest (`digest_fasta_iso`)
- Searches genome for both forward and RC motif patterns
- Extracts fragment sequence as `seq[pos+cut1 .. pos+cut2]`
- Tracks `genome_offset` across chromosomes for genome-wide coordinates

### Probability (`apply_probabilities_iso`, `prob_n_copies_iso`)
- `find_overlaps_iso` — O(n²) overlap detection using recognition-sequence coordinates
- `calculate_prob_iso` — cluster detection via `prev_max > curr_max` boundary condition
- `count_events_iso` — builds compatibility graph, generates all valid fragment combinations
- `apply_probs_iso` — per-fragment probability = count_in_events / total_events

### ISO adapter construction
- `ReadConfig::iso_adapters(l1, l2, filler, m1_cut1, m2_cut1)` — appends `N * cut1` overhangs
  to adapter sequences, matching Python's `create_adapters_iso`

---

## Phase 8: Quality scores and custom adapters

### Quality scores (`--q1`, `--q2`)
- Added `load_q_scores(path) -> Result<Vec<String>>` — reads one Phred string per line
- In `write_reads`: for each read, randomly sample a score string, truncate both the read
  sequence and score to `min(read.len(), score.len())` — mirrors Python's `even_score_seq()`

### Custom adapters (`--a1`, `--a2`)
- Added `load_adapters(path) -> Result<Vec<Adapter>>` — parses whitespace-separated
  `top_seq bot_seq id` lines
- `--a1s`/`--a2s` override the sequencing start offset; default is adapter length

---

## Phase 9: ISO mode bug fixes (achieving exact Python parity)

Three bugs caused systematic ±1 read count differences in ISO mode:

### Bug 1 — Wrong size selection branch
**Root cause:** `main.rs` only called `size_select_low_pass` when `--lp` was explicitly set.
ISO mode with `-u`/`--sd` fell into Gaussian size selection, which gave near-zero probability
for short ALFI fragments (length 34) far from the Gaussian mean (400 bp).

**Fix:** `if cli.low_pass.is_some() || is_iso` — ISO mode always uses low-pass selection.

### Bug 2 — Wrong adjustment formula in `size_select_low_pass`
**Root cause:** Rust computed `total = sum(fragment_comps.values())` after converting to 0/1,
giving `total = count of distinct lengths with any fragments` (e.g. 1 for ALFI). Python
computes `adjustment = n / sum(adj_probs_per_length)` before the 0/1 conversion.

**Fix:** Use `len_dt.values().sum()` (sum of all adj_probs across all lengths) as the
denominator, not the count of distinct lengths.

### Bug 3 — ISO coordinates must be full recognition-sequence coordinates
**Root cause:** Rust stored `start = pos + cut1`, `end = pos + cut2` (trimmed fragment bounds).
Python stores `start = match_start`, `end = match_start + motif_len` (full recognition sequence).
For ALFI (cut1=2, motif_len=36, cut2=36): `end` happened to match but `start` differed by 2,
causing some overlaps to be missed in `find_overlaps_iso`.

**Fix:** Store `start = genome_offset + pos`, `end = genome_offset + pos + ml`.

### Bug 4 — ISO digest must sort by start, not length
**Root cause:** Rust sorted ISO fragments by `length` (all 34 for ALFI = effectively arbitrary
order). Python's `process_df_iso` sorts by `['start']`. The `calculate_prob_iso`
cluster-detection algorithm uses `prev_max > curr_max` to detect cluster boundaries, which
only works when fragments are in genomic order.

**Fix:** `fragments.sort_by_key(|f| f.start)` in `digest_fasta_iso`.

### Result
After all four fixes: all 5 genomes in the ISO test case match Python exactly — identical
fragment counts and read counts per fragment (0 differences).

---

## Validation

### Standard mode (`--m1 ecori --m2 msei`)
- 5-genome `small_test` vignette: all genomes match Python fragment counts and read counts exactly

### ISO mode (`--iso alfi`)
- 5-genome `small_test` vignette: all genomes match Python fragment counts and read counts exactly
- Confirmed that ISO fixes do not affect standard mode output
