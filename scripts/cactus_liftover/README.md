# Cactus liftover of mouse CREs to mammalian species

Standalone bash port of the `generate_lift_regions` Snakemake rule from the paper. Lifts a mouse (`Mus_musculus`) BED of regions to any target species in the Zoonomia 241-mammalian Cactus alignment, then stitches fragmented lifted intervals back into one interval per source region.

## Files

| File | Purpose |
|---|---|
| `generate_lift_regions.sh` | Driver: `halLiftover` + stitching for one or more target species. |
| `stitchHalFrags_v2.R` | Merges lifted fragments per source region; filters by `minMatch` / `maxFrac`. |
| `active_CRE_for_tiling.bed` | Default source mm10 regions (override with `-b`). |

## Requirements

- **halLiftover** from [HAL tools](https://github.com/ComparativeGenomicsToolkit/hal) (easiest via the `cactus` conda env).
- **Rscript** with `GenomicFeatures`, `GenomicRanges`, `rtracklayer`, `stringr`, `dplyr`.
- The 241-mammalian Cactus HAL (`241-mammalian-2020v2.hal`, ~1 TB, **not shipped**) — <https://cglgenomics.ucsc.edu/data/cactus/>.

## Usage

```bash
# one species
./generate_lift_regions.sh -H 241-mammalian-2020v2.hal -s Homo_sapiens -o lifted/

# many (space-separated)
./generate_lift_regions.sh -H 241-mammalian-2020v2.hal \
    -s "Homo_sapiens Rattus_norvegicus Canis_lupus_familiaris" -o lifted/

# many (one per line in a file)
./generate_lift_regions.sh -H 241-mammalian-2020v2.hal -S species_list.txt -o lifted/

# custom source BED
./generate_lift_regions.sh -H 241-mammalian-2020v2.hal -s Homo_sapiens -b my_regions.bed -o lifted/
```

## Options

| Flag | Default | |
|---|---|---|
| `-H` | — | Cactus `.hal` *(required)* |
| `-o` | — | Output dir *(required)* |
| `-s` / `-S` | — | Target species (quoted list) or file with one per line *(need one)* |
| `-b` | `active_CRE_for_tiling.bed` | Source BED |
| `-q` | `Mus_musculus` | Query species in the HAL |
| `-m` | `0.5` | `minMatch` (min covered fraction of source) |
| `-f` | `1.5` | `maxFrac` (max span vs source width) |
| `-l` | `halLiftover` on PATH | halLiftover binary |
| `-r` | alongside script | `stitchHalFrags_v2.R` |

If the query species appears in the target list it is skipped. `-m`/`-f` match the paper values (0.5 / 1.5).

## Output

`<out>/<species>_active_CRE_lift.bed` — one row per successfully lifted source region, columns `chrom start end name`. Regions below `minMatch` coverage or above `maxFrac * source_width` are filtered out by the stitching step.

## Notes

For the full 241-species run, parallelize across species — one job per `-s <SP>` (halLiftover is single-threaded and HAL-I/O heavy).
