# MPRA Barcode Pipeline

End-to-end pipeline for MPRA barcode dictionary construction and quantification. Configured for the **PYS2 oCREv2 / CRE trajectory MPRA**:

- 7 MPRA libraries (`oCREv2_DNA_1-3`, `oCREv2_RNA_1-3`, `oCREv2_plasmid`)
- 7 Bowtie2 indexes (`all_oCRE_cactus`, `compact`, `CRE_trajectory`, `evo_model`, `evo_phyloP`, `oCREs_Gata4_Human`, `pCREs_v2`)
- 6 BC dictionaries merged in Step 5 (`oCRE`, `CRE_traj`, `full_CRE`, `Engreitz_control`, `neg_control`, `pos_control`)
- Sample naming: `{lib}_{material}_{rep}_S{n}` (e.g., `oCREv2_DNA_1_S4`)

## Steps

| | Script | What |
|---|---|---|
| 1 | `01_build_index.sh` | Build Bowtie2 index(es) from CRE FASTA |
| 2 | `02_build_dictionary.sh` | Align subassembly PE reads, extract BC, pileup BC↔CRE |
| 2b | `02b_build_dictionary_compact.sh` | Length-stratified BC↔CRE for the compact library (SE + PEAR) |
| 3 | `03_filter_dictionary.R` | Filter dictionary by min reads, drop polyG |
| 4 | `04_quantify_barcodes.sh` | PEAR-correct BC, link UMI, pileup per mBC |
| 5 | `05_merge_and_count.R` | Merge the 6 dictionaries, join with per-sample counts |
| 5b | `05b_add_dups_in_oCREs.py` | Copy counts to duplicate CRE IDs |
| 6 | `06_analyze_mpra.R` | QC, controls, replicate correlation, oCRE activity |

SGE array wrappers: `sge_01_build_indexes.sh`, `sge_02b_build_dictionary_compact.sh`, `sge_04_quantify_barcodes.sh`.

## Data

Raw **subassembly** FASTQs (Step 2 input) are in GEO **[GSE328309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328309)**. Place them under `DICT_FASTQ_DIR` from `config.sh`. MPRA expression FASTQs (Step 4) are not part of this submission.

## Quick start

All paths are in `config.sh`. Typical run:

```bash
# 1) indexes (SGE) and dictionary
qsub sge_01_build_indexes.sh
bash 02_build_dictionary.sh
Rscript 03_filter_dictionary.R output/dict_pileups/*.txt.gz 10 output/BC_dictionary.txt

# 2) quantify + merge + dedup
qsub sge_04_quantify_barcodes.sh
Rscript 05_merge_and_count.R output/bc_quant output/PYS2
python 05b_add_dups_in_oCREs.py output/PYS2_count_table_wo_zeroes.txt.gz \
       data/dictionaries/all_oCRE_cactusv2_reoriented_marked_dups.txt \
       output/PYS2

# 3) analysis
Rscript 06_analyze_mpra.R output/PYS2_count_table_wo_zeroes_w_dups.txt.gz plots/
```

## FASTQ format

| Run type | R1 | R2 | R3 |
|---|---|---|---|
| Subassembly (Step 2) | CRE fwd | 16 bp barcode | CRE rev |
| MPRA expression (Step 4) | 15 bp mBC fwd | 15 bp mBC rev | 10 bp UMI |

Files: `{sample}_R{1,2,3}_001.fastq.gz`.

## Install

```bash
conda env create -f environment.yml && conda activate mpra-pipeline
# or: conda install -c bioconda bowtie2 samtools pear
#     pip install -r requirements.txt
#     Rscript install_r_packages.R
```

Needs `bowtie2 ≥ 2.5`, `samtools ≥ 1.14`, `pear ≥ 0.9.11`, Python ≥ 3.6 (`pysam`, `biopython`, `pandas`), R ≥ 4.3 (`tidyverse`, `dplyr`, `ggplot2`, `cowplot`, `reshape2`, `ggrastr`, `DescTools`, `castor`, `ape`, `phytools`).

## Compact library (Step 2b)

Compacted CREs from greedy 1 bp-deletion (see `../scripts/oligo_design/compact_optimization/`) end up at varying lengths 40–300 bp. A single pooled Bowtie2 index lets short compacted CREs align to prefixes of longer ones, so `02b_build_dictionary_compact.sh` builds one index **per exact length** and aligns length-matched reads only. For each length *L* ∈ 40..300:

1. `seqtk subseq` the CRE FASTA and paired PEAR-assembled CRE-read + BC FASTQs to length *L*.
2. `bowtie2-build` on just those CREs; `bowtie2 --end-to-end` on the reads.
3. Extract 16 bp BC (`scripts/subassembly_1BC_SE_CRE.py`) and pileup (`scripts/pileup_BC_pair_CRE_SE.R`).

The 261 per-length pileups are concatenated and filtered in Step 3 exactly like the PE output.

Config (`config.sh`): `COMPACT_REF_FASTA`, `COMPACT_PEAR_FASTQ`, `COMPACT_BC_FASTQ`, `COMPACT_MIN_MAPQ` (default 2), optional `COMPACT_LENGTHS`.

If you only have raw PE CRE reads, run PEAR + trim_galore first:

```bash
pear -f *_R1_*.fastq.gz -r *_R3_*.fastq.gz -o compact_pear -j 8
trim_galore compact_pear.assembled.fastq -o trim_galore_fastqs/ --length 40
```

## Outputs

- `output/dict_pileups/*_full_pileup_*.txt.gz` — raw BC-CRE pileups (Step 2)
- `output/*_BC_dictionary.txt` — filtered BC↔CRE map (Step 3)
- `output/bc_quant/*_BC_quant_*.txt.gz` — per-sample mBC quantification (Step 4)
- `output/*_count_table.txt.gz`, `output/*_count_table_wo_zeroes_w_dups.txt.gz` — count matrices (Steps 5 / 5b)
- `plots/*.pdf` and `plots/*_mean_act.txt` — QC + activity figures/tables (Step 6)
