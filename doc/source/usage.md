# Example usage

## Compare regions

```bash
piste compare --input <normalized.h5ad> --region <region_column> --target <region> \
              --gene_set <gene_sets.gmt> --out <output_dir>
```

```r
library(PISTE)
res <- compare(<expression>, <region>, "<gene_sets.gmt>", out = "<output_dir>")
```

## Adjust for cell-type composition

```bash
piste compare --mode adjusted --input <counts.h5ad> \
              --region <region_column> --target <region> \
              --proportions <proportions.csv> \
              --reference <reference.h5ad> --celltype <celltype_column> \
              --gene_set <gene_sets.gmt> --out <output_dir>
```

```r
res <- compare(<counts>, <region>, "<gene_sets.gmt>", mode = "adjusted",
               proportions = <proportions>, reference = <reference>,
               out = "<output_dir>")
```

## Fit adjacent sections jointly

```bash
piste compare -joint_sections --input <normalized.h5ad> --region <region_column> \
                              --section <section_column> --target <region> \
                              --gene_set <gene_sets.gmt> --out <output_dir>
```

```r
res <- compare(<expression>, <region>, "<gene_sets.gmt>",
               section = <section>, joint_sections = TRUE, out = "<output_dir>")
```

## Map spatial activity

```bash
piste map --input <normalized.h5ad> --coords <tissue_positions.csv> \
          --gene_set <gene_sets.gmt> --out <output_dir>
```

```r
res <- map_activity(<expression>, <coords>, "<gene_sets.gmt>", out = "<output_dir>")
```
