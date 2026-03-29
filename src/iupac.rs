use indexmap::IndexMap;

// ---------------------------------------------------------------------------
// Named enzyme lookup
// ---------------------------------------------------------------------------

/// Look up a restriction enzyme by name (case-insensitive) and return its
/// cut-notation string (e.g. "ecori" → "G/AATTC").
///
/// Covers 176 Type IIP enzymes (from resources/type_iip_enzymes.pickle) and
/// 22 Type IIB enzymes (from resources/type_iib_enzymes.pickle).
pub fn lookup_enzyme(name: &str) -> Option<&'static str> {
    match name.to_ascii_lowercase().as_str() {
        // Type IIP
        "aatii"     => Some("GACGT/C"),
        "abasi"     => Some("CNNNNNNNNNNN/NNNNNNNNNG"),
        "acc65i"    => Some("G/GTACC"),
        "acci"      => Some("GT/MKAC"),
        "acli"      => Some("AA/CGTT"),
        "afei"      => Some("AGC/GCT"),
        "aflii"     => Some("C/TTAAG"),
        "afliii"    => Some("A/CRYGT"),
        "agei"      => Some("A/CCGGT"),
        "ahdi"      => Some("GACNNN/NNGTC"),
        "alei"      => Some("CACNN/NNGTG"),
        "alui"      => Some("AG/CT"),
        "alwni"     => Some("CAGNNN/CTG"),
        "apai"      => Some("GGGCC/C"),
        "apali"     => Some("G/TGCAC"),
        "apeki"     => Some("G/CWGC"),
        "apoi"      => Some("R/AATTY"),
        "asci"      => Some("GG/CGCGCC"),
        "asei"      => Some("AT/TAAT"),
        "asisi"     => Some("GCGAT/CGC"),
        "avai"      => Some("C/YCGRG"),
        "avaii"     => Some("G/GWCC"),
        "avrii"     => Some("C/CTAGG"),
        "baegi"     => Some("GKGCM/C"),
        "bamhi"     => Some("G/GATCC"),
        "bani"      => Some("G/GYRCC"),
        "banii"     => Some("GRGCY/C"),
        "bcli"      => Some("T/GATCA"),
        "bfai"      => Some("C/TAG"),
        "bgli"      => Some("GCCNNNN/NGGC"),
        "bglii"     => Some("A/GATCT"),
        "blpi"      => Some("GC/TNAGC"),
        "bmti"      => Some("GCTAG/C"),
        "bsaai"     => Some("YAC/GTR"),
        "bsabi"     => Some("GATNN/NNATC"),
        "bsahi"     => Some("GR/CGYC"),
        "bsaji"     => Some("C/CNNGG"),
        "bsawi"     => Some("W/CCGGW"),
        "bsiei"     => Some("CGRY/CG"),
        "bsihkai"   => Some("GWGCW/C"),
        "bsiwi"     => Some("C/GTACG"),
        "bsli"      => Some("CCNNNNN/NNGG"),
        "bsmbi"     => Some("CGTCTC"),
        "bsp1286i"  => Some("GDGCH/C"),
        "bspdi"     => Some("AT/CGAT"),
        "bspei"     => Some("T/CCGGA"),
        "bsphi"     => Some("T/CATGA"),
        "bsrfi"     => Some("R/CCGGY"),
        "bsrgi"     => Some("T/GTACA"),
        "bsshii"    => Some("G/CGCGC"),
        "bstapi"    => Some("GCANNNN/NTGC"),
        "bstbi"     => Some("TT/CGAA"),
        "bsteii"    => Some("G/GTNACC"),
        "bstni"     => Some("CC/WGG"),
        "bstui"     => Some("CG/CG"),
        "bstxi"     => Some("CCANNNNN/NTGG"),
        "bstyi"     => Some("R/GATCY"),
        "bstz17i"   => Some("GTATAC"),
        "bsu36i"    => Some("CC/TNAGG"),
        "btgi"      => Some("C/CRYGG"),
        "cac8i"     => Some("GCN/NGC"),
        "cviaii"    => Some("C/ATG"),
        "cviki-1"   => Some("RG/CY"),
        "cviqi"     => Some("G/TAC"),
        "ddei"      => Some("C/TNAG"),
        "dpni"      => Some("GA/TC"),
        "dpnii"     => Some("/GATC"),
        "drai"      => Some("TTT/AAA"),
        "draiii"    => Some("CACNNN/GTG"),
        "drdi"      => Some("GACNNNN/NNGTC"),
        "eaei"      => Some("Y/GGCCR"),
        "eagi"      => Some("C/GGCCG"),
        "eco53ki"   => Some("GAG/CTC"),
        "econi"     => Some("CCTNN/NNNAGG"),
        "ecoo109i"  => Some("RG/GNCCY"),
        "ecori"     => Some("G/AATTC"),
        "ecorv"     => Some("GAT/ATC"),
        "fati"      => Some("/CATG"),
        "fnu4hi"    => Some("GC/NGC"),
        "fsei"      => Some("GGCCGG/CC"),
        "fspi"      => Some("TGC/GCA"),
        "haeii"     => Some("RGCGC/Y"),
        "haeiii"    => Some("GG/CC"),
        "hhai"      => Some("GCG/C"),
        "hincii"    => Some("GTY/RAC"),
        "hindiii"   => Some("A/AGCTT"),
        "hinfi"     => Some("G/ANTC"),
        "hinp1i"    => Some("G/CGC"),
        "hpai"      => Some("GTT/AAC"),
        "hpaii"     => Some("C/CGG"),
        "hpy166ii"  => Some("GTN/NAC"),
        "hpy188i"   => Some("TCN/GA"),
        "hpy188iii" => Some("TC/NNGA"),
        "hpy99i"    => Some("CGWCG/"),
        "hpych4iii" => Some("ACN/GT"),
        "hpych4iv"  => Some("A/CGT"),
        "hpych4v"   => Some("TG/CA"),
        "kasi"      => Some("G/GCGCC"),
        "kpni"      => Some("GGTAC/C"),
        "mfei"      => Some("C/AATTG"),
        "mluci"     => Some("/AATT"),
        "mlui"      => Some("A/CGCGT"),
        "msci"      => Some("TGG/CCA"),
        "msei"      => Some("T/TAA"),
        "msli"      => Some("CAYNN/NNRTG"),
        "mspa1i"    => Some("CMG/CKG"),
        "mwoi"      => Some("GCNNNNN/NNGC"),
        "naei"      => Some("GCC/GGC"),
        "nari"      => Some("GG/CGCC"),
        "nb.bbvci"  => Some("CCTCAGC"),
        "nb.bsmi"   => Some("GAATGC"),
        "nb.bsrdi"  => Some("GCAATG"),
        "nb.bsssi"  => Some("CACGAG"),
        "nb.btsi"   => Some("GCAGTG"),
        "ncii"      => Some("CC/SGG"),
        "ncoi"      => Some("C/CATGG"),
        "ndei"      => Some("CA/TATG"),
        "ngomiv"    => Some("G/CCGGC"),
        "nhei"      => Some("G/CTAGC"),
        "nlaiii"    => Some("CATG/"),
        "nlaiv"     => Some("GGN/NCC"),
        "noti"      => Some("GC/GGCCGC"),
        "nrui"      => Some("TCG/CGA"),
        "nsii"      => Some("ATGCA/T"),
        "nspi"      => Some("RCATG/Y"),
        "paci"      => Some("TTAAT/TAA"),
        "paer7i"    => Some("C/TCGAG"),
        "pcii"      => Some("A/CATGT"),
        "pflfi"     => Some("GACN/NNGTC"),
        "pflmi"     => Some("CCANNNN/NTGG"),
        "pluti"     => Some("GGCGC/C"),
        "pmei"      => Some("GTTT/AAAC"),
        "pmli"      => Some("CAC/GTG"),
        "ppumi"     => Some("RG/GWCCY"),
        "pshai"     => Some("GACNN/NNGTC"),
        "psii"      => Some("TTA/TAA"),
        "pspgi"     => Some("/CCWGG"),
        "pspomi"    => Some("G/GGCCC"),
        "pspxi"     => Some("VC/TCGAGB"),
        "psti"      => Some("CTGCA/G"),
        "pvui"      => Some("CGAT/CG"),
        "pvuii"     => Some("CAG/CTG"),
        "rsai"      => Some("GT/AC"),
        "rsrii"     => Some("CG/GWCCG"),
        "saci"      => Some("GAGCT/C"),
        "sacii"     => Some("CCGC/GG"),
        "sali"      => Some("G/TCGAC"),
        "sau96i"    => Some("G/GNCC"),
        "sbfi"      => Some("CCTGCA/GG"),
        "scai"      => Some("AGT/ACT"),
        "scrfi"     => Some("CC/NGG"),
        "sexai"     => Some("A/CCWGGT"),
        "sfci"      => Some("C/TRYAG"),
        "sfii"      => Some("GGCCNNNN/NGGCC"),
        "sfoi"      => Some("GGC/GCC"),
        "sgrai"     => Some("CR/CCGGYG"),
        "smai"      => Some("CCC/GGG"),
        "smli"      => Some("C/TYRAG"),
        "snabi"     => Some("TAC/GTA"),
        "spei"      => Some("A/CTAGT"),
        "sphi"      => Some("GCATG/C"),
        "srfi"      => Some("GCCC/GGGC"),
        "sspi"      => Some("AAT/ATT"),
        "stui"      => Some("AGG/CCT"),
        "styd4i"    => Some("/CCNGG"),
        "styi"      => Some("C/CWWGG"),
        "swai"      => Some("ATTT/AAAT"),
        "taqi"      => Some("T/CGA"),
        "tfii"      => Some("G/AWTC"),
        "tsp45i"    => Some("/GTSAC"),
        "tspmi"     => Some("C/CCGGG"),
        "tspri"     => Some("NNCASTGNN/"),
        "xbai"      => Some("T/CTAGA"),
        "xcmi"      => Some("CCANNNNN/NNNNTGG"),
        "xmni"      => Some("GAANN/NNTTC"),
        "zrai"      => Some("GAC/GTC"),
        // Type IIB
        "ajui"  => Some("NNNNN/NNNNNNNGAANNNNNNNTTGGNNNNNNNNNNN/"),
        "alfi"  => Some("NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN/"),
        "aloi"  => Some("NNNNN/NNNNNNNGAACNNNNNNTCCNNNNNNNNNNNN/"),
        "arsi"  => Some("NNNNN/NNNNNNNNGACNNNNNNTTYGNNNNNNNNNNN/"),
        "baei"  => Some("NNNNN/NNNNNNNNNNACNNNNGTAYCNNNNNNNNNNNN/"),
        "bbri"  => Some("GAAGAC"),
        "bcodi" => Some("NNNNN/NNNNNNNNGAANNNNNNNTTCNNNNNNNNNNN/"),
        "bfii"  => Some("ACTGGG"),
        "bmri"  => Some("ACTGGG"),
        "bpli"  => Some("NNNNN/NNNNNNNNGAGNNNNNTCTNNNNNNNNNN/"),
        "bsai"  => Some("GGTCTC"),
        "bsbi"  => Some("GAAGAC"),
        "bsei"  => Some("NNNNN/NNNNNNNNGAANNNNNNNTTCNNNNNNNNNNN/"),
        "bsgii" => Some("GAGCAG"),  // approximate
        "bsmi"  => Some("GAATGC"),
        "bsmfi" => Some("GGGAC"),
        "bspqi" => Some("GCTCTTC"),
        "btsi"  => Some("GCAGTG"),
        "eari"  => Some("CTCTTC"),
        "esp3i" => Some("CGTCTC"),
        "sapi"  => Some("GCTCTTC"),
        "sfi"   => Some("GGCCNNNNN/NGGCC"),
        _ => None,
    }
}

/// Resolve a list of motif strings, replacing named enzyme identifiers with
/// their cut-notation equivalents (e.g. "ecori" → "G/AATTC").
/// Strings that already contain '/' are passed through unchanged.
/// Returns an error if a name is not recognised.
pub fn resolve_enzyme_names(motifs: &[String]) -> anyhow::Result<Vec<String>> {
    motifs
        .iter()
        .map(|m| {
            if m.contains('/') {
                Ok(m.to_ascii_uppercase())
            } else {
                lookup_enzyme(m)
                    .map(|s| s.to_string())
                    .ok_or_else(|| anyhow::anyhow!("unknown enzyme name '{}' — use cut notation (e.g. G/AATTC) or check spelling", m))
            }
        })
        .collect()
}

/// Complement a single IUPAC character (including '/' for cut-site notation).
fn complement(c: char) -> char {
    match c {
        'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
        'R' => 'Y', 'Y' => 'R', 'S' => 'S', 'W' => 'W',
        'K' => 'M', 'M' => 'K',
        'B' => 'V', 'V' => 'B', 'D' => 'H', 'H' => 'D',
        'N' => 'N', '/' => '/',
        _ => panic!("Unknown IUPAC character: {c}"),
    }
}

/// Reverse complement of an IUPAC sequence.
/// Handles the '/' cut-site marker and all IUPAC ambiguity codes.
/// Mirrors Python's reverse_comp().
pub fn reverse_comp(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

/// Convert a single IUPAC character to its regex equivalent string.
fn iupac_char_to_regex(c: char) -> &'static str {
    match c {
        'A' => "A", 'C' => "C", 'G' => "G", 'T' => "T",
        'R' => "[AG]", 'Y' => "[CT]", 'S' => "[GC]", 'W' => "[AT]",
        'K' => "[GT]", 'M' => "[AC]",
        'B' => "[CGT]", 'D' => "[AGT]", 'H' => "[ACT]", 'V' => "[ACG]",
        'N' => "[ACGT]",
        _ => panic!("Unknown IUPAC character for regex: {c}"),
    }
}

/// Convert a motif string like "GCG/C" into (regex_pattern, cut_position).
///
/// The '/' marks the cut site and is stripped from the regex.
/// Example: "GCG/C" → ("GCGC", 3)
///          "T/TAA" → ("TTAA", 1)
pub fn motif_to_regex(motif: &str) -> (String, usize) {
    let cut_pos = motif
        .chars()
        .position(|c| c == '/')
        .expect("motif must contain a '/' cut-site marker");

    let regex: String = motif
        .chars()
        .filter(|&c| c != '/')
        .map(|c| iupac_char_to_regex(c.to_ascii_uppercase()))
        .collect::<Vec<_>>()
        .join("");

    (regex, cut_pos)
}

/// Return the actual bp length represented by a compiled regex pattern.
/// Bracket groups like `[ACGT]` count as a single bp.
/// Mirrors Python's get_motif_regex_len().
pub fn regex_bp_len(pattern: &str) -> usize {
    let mut len = 0usize;
    let mut in_bracket = false;
    for c in pattern.chars() {
        match c {
            '[' => in_bracket = true,
            ']' => {
                in_bracket = false;
                len += 1;
            }
            _ if !in_bracket => len += 1,
            _ => {}
        }
    }
    len
}

/// For non-palindromic motifs, append the reverse-complement motif (with same cut position).
/// Mirrors Python's check_for_palindrome().
pub fn check_for_palindrome(motifs: &[String]) -> Vec<String> {
    let mut result = motifs.to_vec();
    let mut extras: Vec<String> = Vec::new();

    for motif in motifs {
        let cut_pos = motif
            .chars()
            .position(|c| c == '/')
            .expect("motif must contain '/'");
        let orig: String = motif.chars().filter(|&c| c != '/').collect();
        let revc = reverse_comp(&orig);
        if orig != revc {
            let mut revc_with_cut = revc;
            revc_with_cut.insert(cut_pos, '/');
            extras.push(revc_with_cut);
        }
    }

    result.extend(extras);
    result
}

/// For iso-enzyme mode: parse a single motif with exactly two '/' cut markers,
/// automatically add its reverse complement, and return an IndexMap of
/// regex_pattern → [cut_pos1, cut_pos2].
///
/// Mirrors Python's iupac_motifs_iso().
pub fn iupac_motifs_iso(motif: &str) -> IndexMap<String, [usize; 2]> {
    let mut all_motifs = vec![motif.to_string(), reverse_comp(motif)];
    let mut map = IndexMap::new();
    for m in &all_motifs {
        let slash_positions: Vec<usize> = m
            .char_indices()
            .filter(|(_, c)| *c == '/')
            .map(|(i, _)| i)
            .collect();
        assert_eq!(slash_positions.len(), 2, "iso motif must have exactly two '/' markers: {m}");
        let cut1 = slash_positions[0];
        let cut2 = slash_positions[1] - 1; // Python's -1 adjustment for the removed first '/'
        let regex: String = m
            .chars()
            .filter(|&c| c != '/')
            .map(|c| iupac_char_to_regex(c.to_ascii_uppercase()))
            .collect::<Vec<_>>()
            .join("");
        map.insert(regex, [cut1, cut2]);
    }
    map
}

/// Build an IndexMap of regex_pattern → cut_position for a list of motifs.
/// Preserves insertion order (IndexMap, not HashMap) because the order in which
/// motifs appear in motif_dt determines the trimming sequence in process_fragments.
/// Mirrors Python's iupac_motifs().
pub fn iupac_motifs(motifs: &[String]) -> IndexMap<String, usize> {
    let mut map = IndexMap::new();
    for motif in motifs {
        let (regex, cut_pos) = motif_to_regex(motif);
        map.insert(regex, cut_pos);
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_comp_basic() {
        assert_eq!(reverse_comp("ACGTN"), "NACGT");
    }

    #[test]
    fn test_reverse_comp_iupac() {
        assert_eq!(reverse_comp("RVDAB"), "VTHBY");
    }

    #[test]
    fn test_reverse_comp_complex() {
        assert_eq!(
            reverse_comp("ACYGRYSCCTATTACNNCTNVBBVAYSR"),
            "YSRTBVVBNAGNNGTAATAGGSRYCRGT"
        );
    }

    #[test]
    fn test_motif_to_regex_agei() {
        let (regex, cut) = motif_to_regex("A/CCGGT");
        assert_eq!(regex, "ACCGGT");
        assert_eq!(cut, 1);
    }

    #[test]
    fn test_motif_to_regex_iupac() {
        let (regex, cut) = motif_to_regex("C/YCGRG");
        assert_eq!(regex, "C[CT]CG[AG]G");
        assert_eq!(cut, 1);
    }

    #[test]
    fn test_regex_bp_len() {
        assert_eq!(regex_bp_len("C[CT]CG[AG]G"), 6);
        assert_eq!(regex_bp_len("GCGC"), 4);
        assert_eq!(regex_bp_len("TTAA"), 4);
    }

    #[test]
    fn test_iupac_motifs_accggt() {
        let motifs = vec!["A/CCGGT".to_string()];
        let dt = iupac_motifs(&motifs);
        assert_eq!(dt.get("ACCGGT"), Some(&1));
    }

    #[test]
    fn test_iupac_motifs_iupac() {
        let motifs = vec!["C/YCGRG".to_string()];
        let dt = iupac_motifs(&motifs);
        assert_eq!(dt.get("C[CT]CG[AG]G"), Some(&1));
    }

    #[test]
    fn test_palindrome_gcgc_is_palindromic() {
        // GCGC is palindromic — no extra motif should be added
        let motifs = vec!["GCG/C".to_string()];
        let result = check_for_palindrome(&motifs);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_palindrome_ttaa_is_palindromic() {
        let motifs = vec!["T/TAA".to_string()];
        let result = check_for_palindrome(&motifs);
        assert_eq!(result.len(), 1);
    }
}
