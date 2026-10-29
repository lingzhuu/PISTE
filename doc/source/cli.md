# Command line

SpaPath can be run from a Linux shell without writing R code. One command
reads the data, normalizes it and runs every gene set of a
GMT file. It reads the files and then runs the same workflow as
[`run_spapath()`](api.md#run_spapath) in R.

## Run with Docker

The Docker image contains R, INLA, SpaPath and the `spapath` command. Mount
the folder with your data at `/work`:

```bash
docker build -t spapath .        # once, from the SpaPath repository
docker run --rm -v "$PWD":/work spapath \
  spapath run --input sample.h5ad --region-key layer \
              --gmt pathways.gmt --out results
```

## Run without Docker

After [installing the R package](installation.md), put the command on your
`PATH` once:

```bash
ln -s "$(Rscript -e 'cat(system.file("bin/spapath", package = "PathwayLGM"))')" \
      ~/.local/bin/spapath
spapath --version
```

## Input formats

| Input | Formats |
|-------|---------|
| Expression (`--input`) | AnnData `.h5ad` (`X` or `--layer`), 10x `.h5`, 10x matrix directory (`matrix.mtx`, `features.tsv`, `barcodes.tsv`), `.csv`/`.tsv` with spots as rows (or `--genes-in-rows`), `.rds` matrix |
| Region labels | An `obs` column of the `.h5ad` (`--region-key`), or `--regions`: a table with spot ID and region |
| Section labels (optional) | An `obs` column of the `.h5ad` (`--section-key`), or `--sections`: a table with spot ID and section |
| Gene sets | `--gmt`: GMT file |
| Cell-type proportions (adjusted mode) | `--proportions`: table of spots x cell types, e.g. RCTD weights |
| Reference (adjusted mode) | `--reference`: table of cell types x genes on a linear scale, or single-cell raw counts (`.h5ad`, `.h5` or matrix directory) with `--celltype-key` |

Spots are matched by their IDs across all inputs; spots without a region
label are dropped. Spot coordinates are not needed.

## Observed mode

The input is raw counts, normalized with SCTransform (sctransform v2, with
Pearson residuals clipped as in Seurat). Use `--normalized` if the input is
already normalized.

```bash
spapath run --input 151673.h5ad --region-key layer_guess_reordered \
            --target WM --gmt DLPFC_brain_and_reference.gmt \
            --out results_observed
```

`--target WM` compares white matter with all other spots. Without it, every
region is compared with the others and with each other.

## Adjacent sections

Several sections can be given in one run with `--section-key` or
`--sections`. Each section is normalized separately and fitted on its own.
For adjacent sections of the same tissue block with the same regions,
`--joint-sections` fits them jointly: the regional effects are shared, each
section has its own deviation, and π is computed from the shared part (see
[`observed_fit()`](api.md#observed_fit)).

```bash
spapath run --input 151673_151674.h5ad --region-key layer_guess_reordered \
            --section-key sample_id --joint-sections \
            --target WM --gmt DLPFC_brain_and_reference.gmt \
            --out results_joint
```

## Composition-adjusted mode

The input is raw counts, together with cell-type proportions of each spot and
a reference.

```bash
spapath run --mode adjusted \
            --input 151673.h5ad --region-key layer_guess_reordered \
            --target WM --gmt DLPFC_brain_and_reference.gmt \
            --proportions rctd_proportions.csv \
            --reference dlpfc_snrna.h5ad --celltype-key cellType_broad_hc \
            --out results_adjusted
```

The reference is pooled per cell type and scaled to counts per 10,000. Gene
dispersions are estimated once on all genes with DESeq2.

## Output

| File | Content |
|------|---------|
| `overall_scores.csv` | π score, evidence (`evidence_higher` = P(π > 0.5), `evidence_lower` = P(π < 0.5)) and call of each region against the rest, for every gene set; `p_value` with `--p-value` |
| `pairwise_scores.csv` | The same for every pair of regions |
| `gene_sets.csv` | Number of genes used and status (`done`, `skipped`, `failed`) of each gene set |
| `<gene set>/` | All files written by [`screen_write()`](api.md#screen_write) for that gene set (`<section>/<gene set>/` when sections are fitted one by one) |

## Options

```text
--input INPUT                  Expression: .h5ad, 10x .h5, 10x matrix directory, .csv/.tsv (spots x genes) or .rds. Raw counts unless --normalized.
--layer LAYER                  .h5ad only: read layers/<LAYER> instead of X.
--genes-in-rows                .csv/.tsv only: rows are genes and columns are spots.
--regions REGIONS              Table with spot ID and region label.
--region-key REGION_KEY        .h5ad only: obs column with the region labels.
--sections SECTIONS            Optional table with spot ID and section label; sections are normalised separately and fitted on their own.
--section-key SECTION_KEY      .h5ad only: obs column with the section labels.
--joint-sections               Observed mode: fit the sections jointly, sharing the regional effects (adjacent sections of the same tissue).
--gmt GMT                      Gene sets in GMT format.
--out OUT                      Output directory.
--mode MODE                    observed (Gaussian model of normalised expression) or adjusted (NB model of raw counts adjusted for cell-type composition). [default: observed]
--normalized                   Observed mode: the input is already normalised; skip SCTransform.
--proportions PROPORTIONS      Adjusted mode: table of cell-type proportions, spots x cell types, e.g. RCTD weights.
--reference REFERENCE          Adjusted mode: reference expression, either a table of cell types x genes on a linear scale (e.g. counts per 10,000) or a single-cell .h5ad/.h5/matrix directory of raw counts with --celltype-key.
--celltype-key CELLTYPE_KEY    Adjusted mode: obs column (in a .h5ad reference) or table with cell ID and cell type.
--target TARGET                Compare this region with all other spots, which are pooled into 'rest'.
--baseline BASELINE            Baseline region of the model (default: the first label found).
--min-genes MIN_GENES          Skip gene sets with fewer usable genes than this (at least 2). [default: 5]
--threads THREADS              Adjusted mode: INLA threads, e.g. 4:1 (default: 1:1).
--alpha ALPHA                  A region is called higher or lower when the posterior probability of that direction is at least 1 - ALPHA/2. [default: 0.05]
--p-value                      Also report p_value = 2 min(P(pi > 0.5), P(pi < 0.5)), the doubled posterior probability of the less likely direction (not a frequentist p-value).
--no-pi                        Skip pi scores and evidence; report regional means only.
--help                         Show this help.
```
