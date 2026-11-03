# Command line

## Synopsis

```text
piste [--help] [--version] {compare, compare -joint_sections, map} ...
```

## Subcommands

[`compare`](compare.md)
:   Compare regions on every gene set of a GMT file.

[`compare -joint_sections`](compare_joint_sections.md)
:   Compare regions, fitting adjacent sections of the same tissue jointly.

[`map`](map.md)
:   Map the spatial activity of every gene set, without region labels.

## Input files

| Argument | Formats |
|----------|---------|
| `--input` | AnnData `.h5ad`, 10x `.h5`, 10x matrix directory, `.csv`/`.tsv`, `.rds` matrix |
| `--region`, `--section` | An `obs` column of the `.h5ad` input, or a table with spot ID and label |
| `--coords` | A table with spot ID, `x` and `y`, or a Space Ranger `tissue_positions` file |
| `--gmt` | GMT file |
| `--proportions` | Table of spots x cell types, e.g. RCTD weights |
| `--reference` | Table of cell types x genes, or single-cell raw counts with `--celltype-key` |

Spots are matched by ID across the files. The output folder holds the files
that the R function writes with `out`; see [`compare()`](../api.md#compare)
and [`map_activity()`](../api.md#map_activity).
