# Command line

`piste` runs the analyses from the shell. Every command reads an expression file
and a GMT file of gene sets, and writes the results to the output folder.
`piste <command> --help` lists the options.

These options are the same for every command:

| Option | What it takes | Default |
|--------|---------------|---------|
| `--input` | Expression: .h5ad, 10x .h5, 10x matrix folder, .csv/.tsv (spots × genes) or .rds. Normalized expression; raw counts in the adjusted mode. | |
| `--gene_set` | Gene sets, as a GMT file. | |
| `--out` | Output folder. | |
| `--min-genes` | Skip gene sets with fewer usable genes than this, at least 2. | `5` |

## piste compare

Region comparison, in the observed or the adjusted mode. Runs
[`compare()`](api.md#compare).

```text
piste compare --input FILE --region LABELS \
              --gene_set FILE --out DIR [options]
```

| Option | What it takes | Default |
|--------|---------------|---------|
| `--region` | Region labels: an obs column of the .h5ad, or a table of spot ID and region. | |
| `--target` | Compare this region with all other spots, pooled as `rest`. | |
| `--section` | Section labels, given like `--region`. Each section is fitted on its own. | |
| `--mode` | `observed` (Gaussian model of normalized expression) or `adjusted` (negative binomial model of raw counts, adjusted for cell-type composition). | `observed` |
| `--proportions` | Adjusted mode. Cell-type proportions, spots × cell types, for example RCTD weights. | |
| `--reference` | Adjusted mode. Reference expression: cell types × genes on a linear scale (for example counts per 10,000), or single-cell raw counts (.h5ad, .h5 or matrix folder) with `--celltype`. | |
| `--celltype` | Adjusted mode. Cell type of each cell: an obs column of the .h5ad reference, or a table of cell ID and cell type. | |
| `--threads` | Adjusted mode. INLA threads, for example `4:1`. | `1:1` |
| `--alpha` | A region is called higher or lower when the posterior probability of that direction is at least 1 − α/2. | `0.05` |

## piste compare -joint_sections

Region comparison with adjacent sections of the same tissue fitted jointly, in
the observed mode. Runs [`compare(joint_sections = TRUE)`](api.md#compare).

```text
piste compare -joint_sections --input FILE --region LABELS \
              --section LABELS --gene_set FILE --out DIR [options]
```

`--region`, `--target` and `--alpha` work as in `piste compare`.

| Option | What it takes | Default |
|--------|---------------|---------|
| `--section` | Section labels, given like `--region`. The sections are fitted together and share the regional effects. | |

## piste map

Spatial map of every gene set, without region labels. Runs
[`map_activity()`](api.md#map_activity).

```text
piste map --input FILE --coords FILE \
          --gene_set FILE --out DIR [options]
```

| Option | What it takes | Default |
|--------|---------------|---------|
| `--coords` | Spot coordinates: a table of spot ID, x and y, or a Space Ranger tissue_positions file. | |
| `--no-gene-field` | Leave out the spatial field of each gene: about 50 times faster but less accurate, since patterns of single genes can then appear as regions of the gene set. | |
| `--threads` | INLA threads, outer:inner. | `2:4` |
| `--alpha` | A region of higher or lower activity has joint posterior probability at least 1 − α. | `0.05` |
