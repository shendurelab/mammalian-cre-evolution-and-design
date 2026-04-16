# Cactus liftover of mouse CREs to mammalian species

Standalone reproduction of the `generate_lift_regions` step from the paper's
pipeline. Lifts a mouse (`Mus_musculus`) BED of regions to any target species
in the Zoonomia 241-mammalian Cactus alignment, then stitches fragmented
lifted intervals back into a single interval per source region.

This is the bash port of the `generate_lift_regions` Snakemake rule, packaged
so the per-species Cactus orthologs can be regenerated without Snakemake.

## Contents

| File | Purpose |
| --- | --- |
| `generate_lift_regions.sh` | Driver script: runs `halLiftover` + stitching for one or more target species. |
| `stitchHalFrags_v2.R` | Post-processing in R: merges lifted fragments per source region and filters by `minMatch` / `maxFrac`. |
| `active_CRE_for_tiling_20230404.bed` | Default source mm10 regions used in the paper. Override with `-b` for other inputs. |

## Requirements

- **halLiftover** from the [HAL tools](https://github.com/ComparativeGenomicsToolkit/hal)
  (easiest via the `cactus` conda env).
- **Rscript** with: `GenomicFeatures`, `GenomicRanges`, `rtracklayer`,
  `stringr`, `dplyr`.
- The **241-mammalian Cactus HAL** (`241-mammalian-2020v2.hal`). This file is
  ~1 TB and is **not** shipped in this repo. Download from the Zoonomia /
  UCSC Cactus data page:
  <https://cglgenomics.ucsc.edu/data/cactus/>

## Usage

Single species:

```bash
./generate_lift_regions.sh \
    -H /path/to/241-mammalian-2020v2.hal \
    -s Homo_sapiens \
    -o lifted_hits/
```

Multiple species (quoted, space-separated):

```bash
./generate_lift_regions.sh \
    -H /path/to/241-mammalian-2020v2.hal \
    -s "Homo_sapiens Rattus_norvegicus Canis_lupus_familiaris" \
    -o lifted_hits/
```

Many species from a file (one name per line; e.g. the
`241-mammalian-2020v2.genomes` list that ships with the HAL):

```bash
./generate_lift_regions.sh \
    -H /path/to/241-mammalian-2020v2.hal \
    -S species_list.txt \
    -o lifted_hits/
```

Lift a different source BED:

```bash
./generate_lift_regions.sh \
    -H /path/to/241-mammalian-2020v2.hal \
    -b my_regions.bed \
    -s Homo_sapiens \
    -o lifted_hits/
```

## Options

| Flag | Description | Default |
| --- | --- | --- |
| `-H` | Cactus `.hal` file | *required* |
| `-o` | Output directory | *required* |
| `-s` | Target species, or quoted space-separated list | — |
| `-S` | File with one target species per line | — |
| `-b` | Source BED in the query genome | `active_CRE_for_tiling_20230404.bed` |
| `-q` | Source (query) species in the HAL | `Mus_musculus` |
| `-m` | `minMatch`: min fraction of original interval covered | `0.5` |
| `-f` | `maxFrac`: max fraction of original interval covered | `1.5` |
| `-l` | Path to `halLiftover` binary | `halLiftover` on `$PATH` |
| `-r` | Path to `stitchHalFrags_v2.R` | alongside this script |
| `-h` | Show help | — |

At least one of `-s` or `-S` must be supplied. If the source (query) species
appears in the target list it is skipped.

## Output

For each target species `<SP>` the script writes:

```
<out_dir>/<SP>_active_CRE_lift.bed
```

Columns are `chrom start end name` (one row per successfully lifted source
region, after the `minMatch` / `maxFrac` filters in
`stitchHalFrags_v2.R`). Regions that drop below `minMatch` coverage or
balloon above `maxFrac` of their original width are filtered out.

## How it works

For each target species the script runs:

1. `halLiftover --noDupes --keepExtra <HAL> <query_sp> <src.bed> <target_sp> <out.bed>`
   to project each source interval into the target genome. Without
   post-processing, a single source region typically lands as many small
   fragments.
2. `Rscript stitchHalFrags_v2.R <src.bed> <out.bed> <out.bed> <minMatch> <maxFrac>`
   to group fragments by the source `name` column, keep only sources whose
   summed lifted width is ≥ `minMatch * source_width`, stitch the fragments
   into `[min_start, max_end]` per (chrom, name), and discard stitched
   intervals wider than `maxFrac * source_width`.

The result is one interval per source region per target species, suitable
for downstream tiling / scoring.

## Notes

- `halLiftover` is single-threaded and I/O-heavy on the HAL. For the full
  241-species run the original pipeline parallelized across species on a
  cluster; the equivalent here is to launch one `generate_lift_regions.sh -s
  <SP> ...` job per species.
- `minMatch` and `maxFrac` match the values used in the paper (0.5 and 1.5).
  Adjust with `-m` / `-f` if you want a stricter or looser ortholog
  definition.
