# pCRE v2 — permuted CRE library design

Generates the permuted CRE (pCRE) oligo library used in the paper. For each
input mouse CRE tile, `shuffleCREFasta_chromTiles.py` emits a large
combinatorial set of synthetic sequences that probe the sufficiency,
arrangement, and context of the CRE's predicted TFBS motifs.

This directory is self-contained: every input the script needs ships under
`data/`, and all paths are resolved relative to the script's location, so the
final pCRE TSV can be regenerated bit-for-bit by cloning the repo, installing
the dependencies, and running `python shuffleCREFasta_chromTiles.py`.

## Contents

| File | Purpose |
| --- | --- |
| `shuffleCREFasta_chromTiles.py` | Main driver; produces `JB_tiles_300bp_subtile_shuffled_v4.tsv` in the current working directory. |
| `SequencePremutationTools.py` | Helper class: dinucleotide-preserving shuffling (Altschul & Erickson Eulerian-path algorithm), interval-aware TFBS extraction, TFBS implantation into backgrounds, random-segment shuffling / thripsis. |
| `data/DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa` | Source CRE tiles (mouse 300 bp max-activity subtiles from the DMS screen). |
| `data/TFBS_ProBound_Gata6_Sox17_Foxa2_Klf4_ParEndo_maxTileDMS_20240409.txt` | ProBound TFBS predictions (Gata6 / Sox17 / Foxa2 / Klf4) in ParEndo on the max-activity DMS tiles. |
| `data/TFBS_flanks_minus1_ProBound_..._20240409.txt` | Same TFBS table with each hit trimmed by 1 bp on each flank. |
| `data/TFBS_flanks_minus2_ProBound_..._20240409.txt` | Same TFBS table with each hit trimmed by 2 bp on each flank. |
| `data/max_endo_chromBPnet_320bp_subtile_shuffled_v2.tsv` | Background sequences carried over from the pCRE v1 library. |
| `data/endo_bg_2k.fa` | 2 000-sequence pool of endoderm background candidates; 100 new sequences are sampled from this on top of the v1 set. |

## Requirements

- Python ≥ 3.8
- `pandas`
- `numpy`
- `biopython` (`Bio.SeqIO`, `Bio.Seq`)

Install with e.g. `conda create -n pcre python=3.10 pandas numpy biopython`.

## Usage

From this directory:

```bash
python shuffleCREFasta_chromTiles.py
```

Output is written to the current working directory as
`JB_tiles_300bp_subtile_shuffled_v4.tsv` (three tab-separated columns: header,
sequence, `tile<index>`).

Reproducibility notes:

- The script sets `random.seed(42)` at the top level and also re-seeds
  `numpy.random` and `random` inside `shuffle_RandomSegments` /
  `shuffle_TFBS_implantation` on a per-sequence basis, so the generated TSV
  is deterministic given the inputs shipped under `data/`.
- Biopython's `SeqIO.parse` iterates FASTA records in file order, so the
  order of CREs in the output follows the order in
  `DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa`.

## Library design — what each class of sequence is

For every CRE tile, the script emits:

1. **Originals** (1 per CRE)
   - `<cre>|original_seq` — the unmodified 300 bp tile.
2. **Dinucleotide shuffles — whole CRE** (200 per CRE)
   - `<cre>|dishuffle_all|shuffle_<i>` — full sequence shuffled while
     preserving dinucleotide frequencies. Destroys everything.
3. **TFBS-only shuffles** (200 per CRE)
   - `<cre>|dishuffle_TFBS_full|shuffle_<i>` — dinucleotide-shuffle only
     within each merged TFBS interval; flanking sequence preserved.
4. **CRE "thripsis" — shatter and reassemble** (200 × 3 break-counts per CRE)
   - `<cre>|CRE_thripsis_{5,10,20}breaks_TFBS_full|ss_<i>` — sequence is cut
     at *N* random sites outside TFBS, segments are shuffled, TFBS remain
     intact. Tests the importance of motif *order* and *spacing*.
5. **Reconstitutions** — TFBS implanted into backgrounds at original spacing
   - `<cre>|reconstitution_TFBS_full|bg_<i>` — into a v1 background.
   - `<cre>|reconstitution_TFBS_full|shuffled_bg_<i>` — into a
     dinucleotide-shuffled version of that same background.
   - Tests *sufficiency* of the TFBS motif set + spacing.
6. **Random depositions** — TFBS implanted at random non-overlapping positions
   - `<cre>|random_deposition_all_bkg_TFBS_full|bg_<i>|fixed_deposition_{1,2}`
   - `<cre>|random_deposition_all_bkg_TFBS_full|shuffled_bg_<i>|fixed_deposition_{1,2}`
   - 2 fixed seeds per background — null for "does spacing matter?".
7. **High-replicate random depositions in 6 curated backgrounds**
   - `<cre>|random_deposition_fixed_bkg_TFBS_full|<type>_bg_<i>|deposition_<seed>`
   - `<type>` is `permissive` / `intermediate` / `recalcitrant` — pairs of
     backgrounds chosen from pCRE v1 based on how readily they accommodate
     implanted activity. 100 seeds per background per CRE.
8. **Flanking reconstitutions** (uses the `-1` / `-2` TFBS flank tables)
   - `<cre>|reconstitution_TFBS_minus1_flank|{bg,shuffled_bg}_<i>`
   - `<cre>|reconstitution_TFBS_minus2_flank|{bg,shuffled_bg}_<i>`
   - Tests how much of the immediate TFBS-flanking sequence is required.
9. **Raw backgrounds** (once per library)
   - `bg_<i>|background_seq` and `shuffled_bg_<i>|shuffled_background_seq`

After every implantation, the script re-checks that every TFBS substring is
still present in the resulting sequence and prints an error line if one was
clobbered (sanity check for overlapping TFBS coordinates).

## How the shuffles work (implementation notes)

- Dinucleotide shuffle: exact Altschul–Erickson algorithm (Eulerian-path
  walk on the dinucleotide graph) implemented in
  `SequencePremutationTools.py`. Preserves mononucleotide *and*
  dinucleotide counts, which controls for CpG/GC/dinucleotide composition
  when comparing against originals.
- `SequencePremutationTools` merges overlapping TFBS intervals first, then
  fills the gaps between them, so "shuffle flanks only" and "shuffle TFBS
  only" operate on a clean partition of the sequence.
- TFBS implantation takes an arbitrary background, truncates to the CRE's
  length, and writes each TFBS substring into the background at the
  TFBS's *original* coordinates (reconstitution) or at fresh
  non-overlapping random coordinates (deposition).

## Downstream

The output TSV is the master list of designed sequences. Those were then
scored for predicted marginal footprint with chromBPNet and filtered /
assembled into the Twist oligo order. The scoring scripts and chromBPNet
dependencies are documented in `../README.md` (top-level `oligo_design/`)
and `../../chrombpnet_modeling/`.
