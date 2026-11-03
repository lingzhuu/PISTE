# piste map

Map the spatial activity of every gene set of a GMT file, without region labels. Runs [`map_activity()`](../api.md#map_activity).

## Usage

```text
piste map --input FILE --coords FILE --gmt FILE --out DIR [options]
```

## Arguments

`--input` INPUT
:   Expression: .h5ad, 10x .h5, 10x matrix directory, .csv/.tsv (spots x genes) or .rds. Raw counts unless --normalized.

`--layer` LAYER
:   .h5ad only: read layers/&lt;LAYER&gt; instead of X.

`--genes-in-rows`
:   .csv/.tsv only: rows are genes and columns are spots.

`--coords` COORDS
:   Spot coordinates: table with spot ID, x and y, or a Space Ranger tissue_positions file.

`--gmt` GMT
:   Gene sets in GMT format.

`--out` OUT
:   Output directory.

`--normalized`
:   The input is already normalised; skip SCTransform.

`--no-gene-field`
:   Leave out the spatial field of each gene: about 50 times faster but less accurate, as patterns of single genes can appear as pathway regions.

`--d` D
:   Minimum effect of a region, in units of the noise sd. Default: `0.1`.

`--min-genes` MIN_GENES
:   Skip gene sets with fewer usable genes than this (at least 2). Default: `5`.

`--threads` THREADS
:   INLA threads. Default: `8:1`.

`--alpha` ALPHA
:   A region of higher or lower activity has joint posterior probability at least 1 - ALPHA. Default: `0.05`.
