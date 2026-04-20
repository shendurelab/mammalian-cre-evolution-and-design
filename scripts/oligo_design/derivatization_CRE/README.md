# Derivatized CRE library design

`shuffleCREFasta_chromTiles.py` generates the derivatized CRE oligo library from the paper: per input mouse CRE tile, a large combinatorial set of synthetic sequences that probe the sufficiency, arrangement, and context of the CRE's predicted TFBS motifs.

Self-contained — every input is under `data/`. Run `python shuffleCREFasta_chromTiles.py` to regenerate the library (output: `tiles_300bp_subtile_shuffled.tsv`, columns: header, sequence, `tile<idx>`).

## Inputs (`data/`)

| File | What |
|---|---|
| `DMS_300bp_tiles_..._subtile.fa` | Source CRE tiles (mouse 300 bp max-activity subtiles). |
| `TFBS_ProBound_Gata6_Sox17_Foxa2_Klf4_ParEndo_maxTileDMS.txt` | ProBound TFBS predictions in ParEndo on the max-activity tiles. |
| `TFBS_flanks_minus{1,2}_ProBound_..._maxTileDMS.txt` | Same table with TFBS hits trimmed by 1 / 2 bp each flank. |
| `max_endo_chromBPnet_320bp_subtile_shuffled_v2.tsv` | Background sequences from the pCRE v1 library. |
| `endo_bg_2k.fa` | 2 000-sequence pool; 100 new backgrounds are sampled on top of v1. |

`SequencePremutationTools.py` implements the dinucleotide-preserving shuffle (Altschul–Erickson Eulerian walk) and the TFBS-aware interval primitives used below.

## Requirements

Python ≥ 3.8 with `pandas`, `numpy`, `biopython`.

## Library classes (per CRE)

1. **Original** — 1 per CRE.
2. **Dinucleotide shuffles, whole CRE** (200×).
3. **TFBS-only shuffles** (200×) — shuffle inside merged TFBS intervals only.
4. **CRE thripsis** (200× at 5 / 10 / 20 break counts) — cut outside TFBS, shuffle segments, keep TFBS intact. Tests motif order / spacing.
5. **Reconstitutions** — implant TFBS at their *original* spacing into each background (both real and dinuc-shuffled versions). Tests sufficiency.
6. **Random depositions** — implant TFBS at random non-overlapping positions (2 fixed seeds per background). Null for spacing.
7. **High-replicate depositions in 6 curated backgrounds** (permissive / intermediate / recalcitrant pairs from v1; 100 seeds each).
8. **Flanking reconstitutions** (±1 / ±2 bp flanks around TFBS).
9. **Raw backgrounds** (included once).

Each implant is sanity-checked for TFBS substring preservation. Deterministic given the shipped inputs — `random.seed(42)` at load and per-call re-seeding inside the shuffler/implanter; FASTA iteration order drives output order.
