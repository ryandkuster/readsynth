# readsynth-rs Change Manifest

This document maps every meaningful behavioral change from the Python original to the Rust
rewrite. For each area it records: the Python source of truth, the Rust equivalent, any
intentional deviations, known edge cases, and what to check if output looks wrong.

---

## 1. Fragment digest (standard mode)

| | Python | Rust |
|---|---|---|
| Source | `scripts/digest_genomes.py` | `src/digest.rs :: digest_fasta` |
| Motif iteration order | `for motif1 in args.motif_len.keys()` — dict insertion order | `IndexMap` preserves insertion order; `motif_dt` is built m1 first, m2 second |
| Overlapping site search | `re.finditer('(?=' + pat + ')', seq)` | `overlapping_positions()` advances by 1 byte after each hit |
| Genome-wide coordinates | `begin` accumulates across chromosomes | `genome_offset` accumulates across chromosomes |
| FASTA reading | gzip-transparent via `gzip_test.py` | `flate2::GzDecoder` for `.gz`, plain `BufReader` otherwise |

**Intentional deviation:** Python strips commas from chromosome names; Rust takes the first
whitespace-delimited token after `>`. If chromosome names contain commas, they will differ.

**If fragment counts are wrong:** Check `raw_digests/raw_digest_<genome>.csv`. Column order
should be `chrom, seq, start, end, m1, m2, internal, reverse, revc, length`. Verify that
`start`/`end` match the Python output and that the motif iteration order is consistent.

---

## 2. Fragment processing (trimmed length, orientation)

| | Python | Rust |
|---|---|---|
| Source | `readsynth.py :: process_df()` | `src/digest.rs :: process_fragments` + `compute_trimmed_length` |
| Bidirectional fragments | Duplicated into forward + reverse rows | Same — `is_forward && is_reverse` emits both |
| IUPAC → N | `str.replace('[RYSWKMBDHV]', 'N', regex=True)` | `to_n()` — character match on same set |
| Trimmed length calc | Temporary trim applied, length measured, seq restored | `compute_trimmed_length` — trims a clone, never modifies stored seq |
| Sort order | `df.sort_values(by=['length'])` (unstable quicksort) | `sort_by_key(|f| f.length)` (stable) |

**Known edge case:** Python's `sort_values` is unstable by default. For fragments with the
same trimmed length, row order may differ between Python and Rust. This does not affect read
counts (size selection is length-based) but will produce different row ordering in
`raw_digest_*.csv` for ties.

**If trimmed lengths are wrong:** The trimming loop in `compute_trimmed_length` must iterate
`motif_dt` in insertion order (IndexMap). The `back==0` guard (`if back > 0`) is critical —
Python's `seq[:-0]` returns an empty string, not the original; the guard replicates this.

---

## 3. Fragment digest (ISO mode)

| | Python | Rust |
|---|---|---|
| Source | `scripts/digest_genomes_iso.py` | `src/digest.rs :: digest_fasta_iso` |
| Coordinates stored | `start = match_start`, `end = match_start + motif_len` | Same — full recognition-sequence coords, NOT trimmed fragment coords |
| Fragment sequence | `seq[cut1:cut2]` within the match | `seq[pos+cut1 .. pos+cut2]` |
| Sort order | `df.sort_values(by=['start'])` | `fragments.sort_by_key(|f| f.start)` |
| RC motif handling | Both forward + RC patterns added via `iupac_motifs_iso` | Same; `is_rc = mi == 1` (second key in IndexMap) |

**Critical:** The `calculate_prob_iso` cluster-detection algorithm (`prev_max > curr_max`)
requires fragments to be in genomic order. If sort order is wrong, cluster boundaries will be
miscounted and probabilities will be incorrect.

**Critical:** The overlap-detection algorithm in `find_overlaps_iso` uses `(start, end)` from
the stored `ProcessedFragment`. These must be full recognition-sequence coordinates (not cut
coordinates) to match Python. Storing trimmed coordinates was a prior bug — see Phase 9 in
`DEVELOPMENT_LOG.md`.

**If ISO fragment counts are wrong:** Verify that `start` in `raw_digest_*.csv` equals the
regex match start (not `match_start + cut1`), and that rows are sorted by `start`.

---

## 4. Per-fragment probability (standard mode)

| | Python | Rust |
|---|---|---|
| Source | `scripts/prob_n_copies.py` | `src/prob_copies.rs :: apply_probabilities` |
| Probability formula | `(c ** (internal + 1)) * ((1 - c) ** (1 - 1))` simplified | `fragment_probability(c, internal) = c.powi(internal + 1)` |
| adj_prob | `probability * args.comp` | `prob * abundance` |
| len_freqs | `df[df['length'] == i]['adj_prob'].sum()` for each `i` | Same accumulation, 0.0 for absent lengths |

**If probabilities are wrong:** Check that `adj_prob` in `individual_counts/counts_*.csv`
matches Python. A mismatch here flows into size selection and final read counts.

---

## 5. Per-fragment probability (ISO mode)

| | Python | Rust |
|---|---|---|
| Source | `scripts/prob_n_copies_iso.py` | `src/prob_copies.rs :: apply_probabilities_iso` |
| Overlap detection | `find_overlaps()` — O(n²), uses `(start, end)` coordinates | `find_overlaps_iso()` — same logic |
| Cluster detection | `calculate_prob()` — `prev_max > curr_max` boundary | `calculate_prob_iso()` — same logic |
| Event counting | `count_events()` — deduplicates via `set(tuple(sorted(...)))` | `count_events_iso()` — deduplicates via `HashSet<Vec<usize>>` after sorting |
| Probability assignment | `count_in_events / total_events` | `apply_probs_iso()` — same |

**If ISO probabilities differ:** The most likely causes are (1) fragment order (must be by
genomic start), (2) coordinate system (must be full recognition-sequence coords). Both were
bugs that have been fixed — see Phase 9 in `DEVELOPMENT_LOG.md`.

---

## 6. Size selection — Gaussian

| | Python | Rust |
|---|---|---|
| Source | `scripts/size_selection.py :: get_scale_gauss + get_draw_dict` | `src/size_select.rs :: size_select_gaussian` |
| `gauss_pdf` formula | `(1 / sd * sqrt(2π)) * exp(...)` — non-standard precedence | Replicated exactly: `(1.0 / sd * sqrt(2π)) * exp(...)` |
| `x_dist` | `min(x, abs(max - sd), max - x)` | Same; uses `unsigned_abs()` for the middle term |
| Sage recovery factor | `min(draw, data) * 0.65` | Same |
| `fragment_comps` | `draw_dt[i] / len_dt[i]`, 0 if either is 0 or NaN | Same, with `is_nan()` guard |

**If Gaussian size selection produces wrong counts:** The `gauss_pdf` formula is intentionally
non-standard. Do not "fix" it to `1/(sd*sqrt(2π))` — the constants cancel in `scale_by` but
must match Python exactly for floating-point reproducibility.

---

## 7. Size selection — Low-pass / ISO

| | Python | Rust |
|---|---|---|
| Source | `readsynth.py` lines 899–906 (ISO branch) | `src/size_select.rs :: size_select_low_pass` |
| `fragment_comps` | `1 if sum_adj_prob > 0 else 0` per length | Same |
| `adjustment` | `n / sum(adj_probs_per_length)` — sum BEFORE 0/1 conversion | `n / len_dt.values().sum()` — same |
| Triggered by | `args.iso or args.lp` | `cli.low_pass.is_some() \|\| is_iso` |

**Critical:** `adjustment` must use the sum of raw adj_probs from `len_dt`, not the count of
lengths that have fragments. Using the wrong denominator was a prior bug — see Phase 9.

**Critical:** ISO mode must always use low-pass even when `-u`/`--sd` are supplied (those args
only set `max_len`, not the selection distribution). This was a prior bug — see Phase 9.

---

## 8. Size selection — Custom distribution (`-d`)

| | Python | Rust |
|---|---|---|
| Source | `readsynth.py :: check_custom_distribution` + `size_selection.py :: get_custom_dist` | `src/size_select.rs :: load_custom_dist + size_select_custom` |
| `max_len` source | `max(json_keys)` | Same — `load_custom_dist` returns `(dist, max_len)` |
| Draw counts | JSON value for each length, 0 if absent | Same |
| `adjustment` | `n / sum(draw_dt.values())` | Same |
| `fragment_comps` | `draw_dt[i] / len_dt[i]`, 0 if either 0/NaN | Same |

**Note:** When `-d` is provided, `-u`/`--sd` are not required. `max_len` comes entirely from
the JSON. The JSON is loaded once at startup (before the digest) so `max_len` is available
for the digest's search window.

---

## 9. Read writing

| | Python | Rust |
|---|---|---|
| Source | `scripts/write_reads.py` | `src/write_reads.rs :: write_reads` |
| R1 ligated string | `a1[0] + seq + a2[1]` | `format!("{}{}{}", a1.0, frag.seq, a2.1)` |
| R2 ligated string | `a2[0] + revc + a1[1]` | `format!("{}{}{}", a2.0, frag.revc, a1.1)` |
| Read slice | `ligated[a?s : a?s + l?]` | `slice_pad(ligated, a?s, a?s+l?, l?, filler)` |
| Filler padding | `str.ljust(l, args.e)` | `while out.len() < len { out.push(filler) }` |
| Quality (no q-scores) | `'I' * l` | `"I".repeat(l)` |
| Quality (with q-scores) | `even_score_seq`: truncate both seq and score to `min(len(seq), len(score))` | Same — `let len = r1_read.len().min(s.len())` |
| Header format | `@{r_no}:{idx}:{full}:{length}:{a1_id}:{a2_id}:{genome_name} 1` | Same |

**Default adapter offsets:** `a1s = 62`, `a2s = 58` (length of the default Illumina adapter
sequences). Custom adapter files default to `a1s = len(adapters[0].top)`.

**If reads look wrong:** Check that `seq` and `revc` in the counts CSV are the full
(untrimmed) sequences. Trimmed length is only used for size selection; the actual bases written
are the full motif-to-motif sequence.

---

## 10. Error simulation (`apply_error`)

| | Python | Rust |
|---|---|---|
| Source | `src/apply_error.cpp` | `src/apply_error.rs :: apply_error` |
| Trigger | Separate manual step after fastq generation | Automatic Stage 4 when both `--q1` and `--q2` are provided |
| Output files | `error_sim_metagenome_R1/R2.fastq` (replaces originals after removing simulated) | Written alongside originals: `error_sim_metagenome_R1/R2.fastq` |
| Probability formula | `pow(10, -idx/10.0)` where `idx = phred_char - '!'` | `10f64.powf(-(idx as f64) / 10.0)` — same |
| Mutation | `rand() < prob * (RAND_MAX+1.0)` | `rng.gen::<f64>() < prob` — statistically equivalent |
| N bases | Passed through unchanged | Same |
| Mutation target | Uniform over the 3 non-self bases | Same |

**Intentional deviation:** Python removes the clean `sim_metagenome_*.fastq` after running
`apply_error`, leaving only the error-mutated files. Rust keeps both — the clean fastq and
the error-mutated fastq are written side-by-side. This makes it easier to compare the two.

**If mutation rate is wrong:** Verify the Phred character mapping. ASCII `'!'` = Phred 0 =
100% error, `'5'` = Phred 20 = 1% error, `'I'` = Phred 40 = 0.01% error. The index is
`q_char as u8 - 33`.

---

## 11. Parameters: Python vs. Rust flag names

| Python flag | Rust flag | Notes |
|-------------|-----------|-------|
| `-m1` | `--m1` | Now long-form |
| `-m2` | `--m2` | Now long-form |
| `-iso` | `--iso` | Now long-form |
| `-sd` | `--sd` | Now long-form |
| `-lp` | `--lp` | Now long-form |
| `-a1` | `--a1` | Now long-form |
| `-a2` | `--a2` | Now long-form |
| `-a1s` | `--a1s` | Now long-form |
| `-a2s` | `--a2s` | Now long-form |
| `-q1` | `--q1` | Now long-form |
| `-q2` | `--q2` | Now long-form |
| `-l2` | `--l2` | Now long-form |
| `-test` | `--test` | Now long-form |
| `-g` | `-g` | Unchanged |
| `-o` | `-o` | Unchanged |
| `-l`, `-l1` | `-l` | `-l1` alias removed |
| `-n` | `-n` | Unchanged |
| `-u` | `-u` | Unchanged |
| `-x` | `-x` | Unchanged |
| `-d` | `-d` | Unchanged |
| `-c` | `-c` | Unchanged |
| `-e` | `-e` | Unchanged |
| _(none)_ | `-t` | New: thread count, default 1 |

---

## 12. Output files

| File | Present when | Notes |
|------|-------------|-------|
| `raw_digests/raw_digest_<genome>.csv` | Always | All possible fragments before probability weighting |
| `individual_counts/counts_<genome>.csv` | Always | Per-fragment adj_prob and final read counts |
| `combined.csv` | Always | All genomes' len_freqs with size-selection counts |
| `metagenome_summary.csv` | Always | Per-genome read totals, site counts, abundances |
| `sim_metagenome_R1.fastq` | Unless `--test` | Clean simulated reads (quality = `I` if no `--q1`) |
| `sim_metagenome_R2.fastq` | Unless `--test` | Clean simulated reads |
| `error_sim_metagenome_R1.fastq` | `--q1` + `--q2` + not `--test` | Phred-mutated reads |
| `error_sim_metagenome_R2.fastq` | `--q1` + `--q2` + not `--test` | Phred-mutated reads |

**Note:** Python also produced per-genome histogram PNGs (`individual_histograms/`) and
combined histogram PDFs (`fragment_distributions.pdf`, `read_distributions.pdf`). These are
not produced by the Rust rewrite.

---

## 13. Validated test cases

| Test | Parameters | Result |
|------|-----------|--------|
| Standard 5-genome | `--m1 ecori --m2 msei -n 1000000 -u 400 --sd 100 -l 150` | Fragment counts and read counts match Python exactly across all 5 genomes |
| ISO 5-genome | `--iso alfi -n 1000000 -u 400 --sd 100 -l 150` | Fragment counts and read counts match Python exactly across all 5 genomes |
| ISO does not affect standard | Run standard after ISO fixes applied | Confirmed — 0 differences from pre-fix standard output |
| Q-score attachment | `--q1 miseq_151bp_R1 --q2 miseq_151bp_R2` | 0 score/seq length mismatches, 0 invalid Phred chars across 999,397 reads |
| Error simulation (Q20) | `--q1 20_phred_score --q2 20_phred_score` | Observed mutation rate 1.00% across 149,851,474 bases (expected 1.00%) |
| Filler base (`-e`) | `-e N` vs default `-e G`, `-l 300` | 2,894 reads differ between runs, all correctly using specified filler |
| Custom distribution (`-d`) | JSON with 10:100:1000 count ratios | Read counts reflect input ratios exactly |
| Custom adapters (`--a1`/`--a2`) | 2-entry A1 file, 1-entry A2 file | Correct adapter IDs in all read headers, correct read length |
| Vignette: `small_test` | Both standard and ISO | Stored in `vignettes/small_test/` |
