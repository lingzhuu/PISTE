# Command line

SpaPath can be run from a Linux shell without writing R code. One command
reads the data, normalizes it, builds the mesh and runs every gene set of a
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
| Coordinates | `obsm['spatial']` of the `.h5ad`, or `--coords`: a table with spot ID, `x` and `y`, or Space Ranger `tissue_positions.csv` / `tissue_positions_list.csv` |
| Region labels | An `obs` column of the `.h5ad` (`--region-key`), or `--regions`: a table with spot ID and region |
| Gene sets | `--gmt`: GMT file |
| Cell-type proportions (adjusted mode) | `--proportions`: table of spots x cell types, e.g. RCTD weights |
| Reference (adjusted mode) | `--reference`: table of cell types x genes on a linear scale, or single-cell raw counts (`.h5ad`, `.h5` or matrix directory) with `--celltype-key` |

Spots are matched by their IDs across all inputs; spots without coordinates
or a region label are dropped.

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
| `overall_scores.csv` | π score, evidence, δ evidence and call of each region against the rest, for every gene set |
| `pairwise_scores.csv` | π score, evidence, δ evidence and call for every pair of regions, for every gene set |
| `gene_sets.csv` | Number of genes used and status (`done`, `skipped`, `failed`) of each gene set |
| `<gene set>/` | All files written by [`screen_write()`](api.md#screen_write) for that gene set |

## Options

```text
--input INPUT                  Expression: .h5ad, 10x .h5, 10x matrix directory, .csv/.tsv (spots x genes) or .rds. Raw counts unless --normalized.
--layer LAYER                  .h5ad only: read layers/<LAYER> instead of X.
--genes-in-rows                .csv/.tsv only: rows are genes and columns are spots.
--coords COORDS                Spot coordinates: table with spot ID and x, y, or Space Ranger tissue_positions.csv.
--spatial-key SPATIAL_KEY      .h5ad only: obsm key of the coordinates. [default: spatial]
--regions REGIONS              Table with spot ID and region label.
--region-key REGION_KEY        .h5ad only: obs column with the region labels.
--gmt GMT                      Gene sets in GMT format.
--out OUT                      Output directory.
--mode MODE                    observed or adjusted. [default: observed]
--normalized                   Observed mode: the input is already normalised; skip SCTransform.
--proportions PROPORTIONS      Adjusted mode: table of cell-type proportions, spots x cell types.
--reference REFERENCE          Adjusted mode: table of cell types x genes, or single-cell raw counts with --celltype-key.
--celltype-key CELLTYPE_KEY    Adjusted mode: obs column of a .h5ad reference, or table with cell ID and cell type.
--target TARGET                Compare this region with all other spots, pooled into 'rest'.
--baseline BASELINE            Baseline region of the model (default: the first label found).
--min-genes MIN_GENES          Skip gene sets with fewer usable genes than this (at least 2). [default: 5]
--h H                          Spatial length scale, in coordinates scaled to the unit square. [default: 0.3]
--max-edge MAX_EDGE            Largest mesh triangle edge, one value or inner,outer (default: h/3,h).
--cutoff CUTOFF                Minimum distance between mesh nodes (default: h/9).
--mesh-offset MESH_OFFSET      Inner and outer mesh extension around the spots (default: h/2,1.5h).
--threads THREADS              INLA threads, e.g. 4:1 (default: 1:1).
--delta DELTA                  Minimal effect: a region is called higher when P(pi > 0.5 + DELTA) > LEVEL and lower when P(pi < 0.5 - DELTA) > LEVEL. [default: 0.05]
--level LEVEL                  Posterior probability required for a call (see --delta). [default: 0.9]
--no-pi                        Skip pi scores and evidence; report regional means only.
```
