# Mammalian CRE Evolution and Design

> **Note:** This repository is under active construction and will be updated with additional code, data, and documentation as the manuscripts are finalized.

Code and data to reproduce the analyses and figures from our two companion studies on the **evolution** and **design** of mammalian cis-regulatory elements (CREs), using MPRA, ATAC-seq, and deep-learning models of chromatin accessibility across mammals.

## Papers

1. Retracing and rewriting the evolutionary trajectories of mammalian developmental enhancers ([bioRxiv](https://www.biorxiv.org/content/10.64898/2026.04.20.719714v1))
2. Multi-scale dissection, compaction and derivatization of mammalian developmental enhancers ([bioRxiv](https://www.biorxiv.org/content/10.64898/2026.04.20.719625v1))

## Repository layout

```
mammalian-cre-evolution-and-design/
├── data/                       # MPRA count tables, ATAC/RNA-seq, phyloP, MSAs, predictions
├── extdata/                    # Pre-computed R binary data for figure scripts
├── imgs/                       # Phylogenetic silhouette images (from https://www.phylopic.org/)
├── scripts/
│   ├── make_figure1.R … make_figure5.R          # Main figures (paper 1)
│   ├── make_supp_figure*.R, make_ext_data_fig1.R # Supplementary figures (paper 1)
│   ├── utils/                  # Shared plotting, analysis, tree, MSA helpers
│   ├── cactus_liftover/        # Ortholog extraction from 241-way mammalian cactus
│   ├── oligo_design/           # Synthetic CRE library design
│   │   ├── evolvability_optimization/   # DMS + pairwise ancestral SNP rescue
│   │   ├── compact_optimization/        # Iterative 1bp deletion
│   │   └── derivatization_CRE/          # TFBS shuffling / reconstitution / deposition
│   └── chrombpnet_modeling/    # Cell-type chromBPNet model training from scATAC
└── mpra-barcode-quant-pipeline/      # End-to-end MPRA BC dictionary + quantification pipeline
```

## What's in each directory

- **`scripts/make_figure*.R`** – R scripts that regenerate the main and supplementary figure panels **for paper 1** (the evolution paper), from files in `data/` and `extdata/`. Shared helpers live in `scripts/utils/`. Figure scripts for paper 2 will be added as the manuscript is finalized.
- **`mpra-barcode-quant-pipeline/`** – Configured for the PYS2 oCREv2 MPRA. Takes raw FASTQs → BC-CRE dictionary → BC quantification → per-CRE activity. See [`mpra-barcode-quant-pipeline/README.md`](mpra-barcode-quant-pipeline/README.md).
- **`scripts/oligo_design/`** – Generates the synthetic CRE sequences tested in MPRA, using a trained chromBPNet model. See [`scripts/oligo_design/README.md`](scripts/oligo_design/README.md).
- **`scripts/chrombpnet_modeling/`** – Trains cell-type-specific chromBPNet models from mEB scATAC-seq used by the oligo-design scripts. See [`scripts/chrombpnet_modeling/README.md`](scripts/chrombpnet_modeling/README.md).
- **`scripts/cactus_liftover/`** – Pulls orthologous CRE sequences across 241 mammals from UCSC cactus alignments.

## Reproducing results

```bash
# 1. Regenerate a figure (e.g. Fig 4)
Rscript scripts/make_figure4.R

# 2. Re-run MPRA processing from FASTQs
cd mpra-barcode-quant-pipeline
bash 01_build_index.sh && bash 02_build_dictionary.sh && Rscript 03_filter_dictionary.R …
bash 04_quantify_barcodes.sh && Rscript 05_merge_and_count.R … && Rscript 06_analyze_mpra.R …

# 3. Retrain chromBPNet models (GPU required)
cd scripts/chrombpnet_modeling
./sm_script.sh                    # preprocessing
qsub train_gpu.sh                 # bias models
qsub train_factorized_gpu.sh      # factorized chromBPNet models

# 4. Regenerate oligo designs (requires trained chromBPNet)
python scripts/oligo_design/evolvability_optimization/dms_marginal_footprinting.py
python scripts/oligo_design/compact_optimization/compact_CREs.py
python scripts/oligo_design/derivatization_CRE/shuffleCREFasta_chromTiles.py
```

## Dependencies

- **R** (≥4.3) with tidyverse, ggplot2, cowplot, ggtree, Biostrings, GenomicRanges, castor, phytools, ape, DescTools, ggrastr
- **Python** (≥3.7) with tensorflow, chrombpnet, pyfaidx, pyBigWig, biopython, pandas, numpy, python-Levenshtein
- **CLI**: bowtie2, samtools, bedtools, pear, macs2, kentUtils, snakemake, mafft

See subdirectory READMEs for component-specific dependencies.

## Data

Large data files (MPRA count tables, TFBS predictions, cactus alignments, ATAC/RNA-seq signals) are in `data/`. See individual figure scripts for the specific files each one consumes.

## License

Released under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license. PhyloPic silhouettes under `imgs/` retain their original individual licenses (see [`imgs/README.md`](imgs/README.md)).
