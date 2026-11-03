# piste compare -joint_sections

Compare regions on every gene set of a GMT file, fitting adjacent sections of the same tissue jointly (observed mode). Runs [`compare(joint_sections = TRUE)`](../api.md#compare).

## Usage

```text
piste compare -joint_sections --input FILE --region LABELS --section LABELS --gmt FILE --out DIR [options]
```

## Arguments

`--input` INPUT
:   Expression: .h5ad, 10x .h5, 10x matrix directory, .csv/.tsv (spots x genes) or .rds. Raw counts unless --normalized.

`--layer` LAYER
:   .h5ad only: read layers/&lt;LAYER&gt; instead of X.

`--genes-in-rows`
:   .csv/.tsv only: rows are genes and columns are spots.

`--region` REGION
:   Region labels: an obs column of the .h5ad, or a table with spot ID and region label.

`--section` SECTION
:   Section labels: an obs column of the .h5ad, or a table with spot ID and section label. The sections are fitted jointly, sharing the regional effects.

`--gmt` GMT
:   Gene sets in GMT format.

`--out` OUT
:   Output directory.

`--normalized`
:   The input is already normalised; skip SCTransform.

`--target` TARGET
:   Compare this region with all other spots, which are pooled into 'rest'.

`--baseline` BASELINE
:   Baseline region of the model (default: the first label found).

`--min-genes` MIN_GENES
:   Skip gene sets with fewer usable genes than this (at least 2). Default: `5`.

`--alpha` ALPHA
:   A region is called higher or lower when the posterior probability of that direction is at least 1 - ALPHA/2. Default: `0.05`.

`--p-value`
:   Also report p_value = 2 min(P(pi > 0.5), P(pi < 0.5)), the doubled posterior probability of the less likely direction (not a frequentist p-value).

`--no-pi`
:   Skip pi scores and evidence; report regional means only.
