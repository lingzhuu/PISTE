# Example usage

The examples use the DLPFC section 151673, with the cortical layers in
`layer_guess_reordered` and the gene sets in `DLPFC_brain_and_reference.gmt`.
All arguments are listed in the [command line](cli/index.md) and
[R API](api.md) reference.

## Compare regions

White matter against the other layers:

```bash
piste compare --input 151673.h5ad --region layer_guess_reordered --target WM \
              --gmt DLPFC_brain_and_reference.gmt --out results
```

```r
library(PISTE)
res <- compare(counts, region, "DLPFC_brain_and_reference.gmt", out = "results")
```

## Adjust for cell-type composition

```bash
piste compare --mode adjusted --input 151673.h5ad \
              --region layer_guess_reordered --target WM \
              --proportions rctd_proportions.csv \
              --reference dlpfc_snrna.h5ad --celltype-key cellType_broad_hc \
              --gmt DLPFC_brain_and_reference.gmt --out results_adjusted
```

```r
res <- compare(counts, region, "DLPFC_brain_and_reference.gmt", mode = "adjusted",
               proportions = proportions, reference = reference,
               out = "results_adjusted")
```

## Fit adjacent sections jointly

```bash
piste compare -joint_sections --input 151673_151674.h5ad --region layer_guess_reordered \
                              --section sample_id --target WM \
                              --gmt DLPFC_brain_and_reference.gmt --out results_joint
```

```r
res <- compare(counts, region, "DLPFC_brain_and_reference.gmt",
               section = section, joint_sections = TRUE, out = "results_joint")
```

## Map spatial activity

```bash
piste map --input 151673.h5ad --coords tissue_positions.csv \
          --gmt DLPFC_brain_and_reference.gmt --out results_map
```

```r
res <- map_activity(counts, coords, "DLPFC_brain_and_reference.gmt", out = "results_map")
```
