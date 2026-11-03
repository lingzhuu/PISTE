# piste compare

Compare regions on every gene set of a GMT file, in observed or composition-adjusted mode. Runs [`compare()`](../api.md#compare).

## Usage

```text
piste compare --input FILE --region LABELS --gmt FILE --out DIR [options]
```

## Arguments

`--input` INPUT
:   Expression: .h5ad, 10x .h5, 10x matrix directory, .csv/.tsv (spots x genes) or .rds. Raw counts unless --normalized.

`--region` REGION
:   Region labels: an obs column of the .h5ad, or a table with spot ID and region label.

`--section` SECTION
:   Optional section labels: an obs column of the .h5ad, or a table with spot ID and section label; sections are normalised separately and fitted on their own.

`--gmt` GMT
:   Gene sets in GMT format.

`--out` OUT
:   Output directory.

`--mode` MODE
:   observed (Gaussian model of normalised expression) or adjusted (NB model of raw counts adjusted for cell-type composition). Default: `observed`.

`--normalized`
:   Observed mode: the input is already normalised; skip SCTransform.

`--proportions` PROPORTIONS
:   Adjusted mode: table of cell-type proportions, spots x cell types, e.g. RCTD weights.

`--reference` REFERENCE
:   Adjusted mode: reference expression, either a table of cell types x genes on a linear scale (e.g. counts per 10,000) or a single-cell .h5ad/.h5/matrix directory of raw counts with --celltype.

`--celltype` CELLTYPE
:   Adjusted mode: obs column (in a .h5ad reference) or table with cell ID and cell type.

`--target` TARGET
:   Compare this region with all other spots, which are pooled into 'rest'.

`--min-genes` MIN_GENES
:   Skip gene sets with fewer usable genes than this (at least 2). Default: `5`.

`--threads` THREADS
:   Adjusted mode: INLA threads, e.g. 4:1 (default: 1:1).

`--alpha` ALPHA
:   A region is called higher or lower when the posterior probability of that direction is at least 1 - ALPHA/2. Default: `0.05`.
