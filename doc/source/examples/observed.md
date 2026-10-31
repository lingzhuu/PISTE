# Observed mode on DLPFC

<!-- Generated from doc/examples/observed.Rmd by doc/examples/render.R; edit the .Rmd. -->

This example tests pathways in white matter (WM) against the cortical layers
of the human dorsolateral prefrontal cortex (DLPFC) Visium section 151673
from spatialLIBD (Maynard et al., 2021). The observed mode fits a Gaussian
model to normalized expression and reports the total regional difference in
pathway activity, including differences that come from cell composition.

## Requirements

This example also uses:

| Package | Source | Used for |
|---------|--------|----------|
| spatialLIBD | Bioconductor | Provides the DLPFC Visium data and its layer annotation, and plots spots on the H&E image |


``` r
install.packages("BiocManager")
BiocManager::install("spatialLIBD")
```


## 1. Load the data


``` r
library(spatialLIBD)
library(SpatialExperiment)
library(PathwayLGM)

spe <- fetch_data(type = "spe")
spe <- spe[, spe$sample_id == "151673"]

# Keep spots with a layer annotation (Layer1 to Layer6, WM).
spe <- spe[, !is.na(spe$layer_guess_reordered)]

counts <- t(counts(spe))                       # spots x genes, raw counts
colnames(counts) <- make.unique(rowData(spe)$gene_name)
region <- as.character(spe$layer_guess_reordered)
table(region)
```

```
#> region
#> Layer1 Layer2 Layer3 Layer4 Layer5 Layer6     WM 
#>    273    253    989    218    673    692    513
```

`counts` and `region` follow the same spot order.

The H&E image of the section and its manual layer annotation, Layer1 to
Layer6 and white matter:


``` r
plot(imgRaster(spe, sample_id = "151673", image_id = "lowres"))
```

![plot of chunk he-image](figures/observed-he-image-1.png)


``` r
vis_clus(spe, sampleid = "151673", clustervar = "layer_guess_reordered",
         colors = libd_layer_colors)
```

![plot of chunk layers](figures/observed-layers-1.png)

## 2. Gene sets

20 MSigDB gene sets are stored as a GMT file:


``` r
gmt <- "benchmark/analysis/DLPFC/data/DLPFC_brain_and_reference.gmt"
```

## 3. Run SpaPath


``` r
res <- compare_regions(counts, region, gmt, baseline = "WM",
                   out = "results_observed")
```

`compare_regions()` normalizes the counts with SCTransform and fits every gene
set, with WM as the baseline region. Results are written to
`results_observed`. See [API](../api.md#compare_regions) for all options.

## 4. Read the results


``` r
wm <- res$overall[res$overall$region == "WM", ]
wm <- wm[order(wm$evidence, decreasing = TRUE), ]
rownames(wm) <- NULL
wm[, c("gene_set", "genes", "score", "evidence")]
```

```
#>                                                         gene_set genes
#> 1                                         REACTOME_MRNA_SPLICING   210
#> 2                             GOBP_MONOCARBOXYLIC_ACID_TRANSPORT   126
#> 3                               REACTOME_SPHINGOLIPID_METABOLISM    96
#> 4                                             GOCC_MYELIN_SHEATH    42
#> 5               GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM    23
#> 6                            REACTOME_NUCLEOTIDE_EXCISION_REPAIR   107
#> 7              GOBP_VERY_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS    32
#> 8                                        GOBP_MYELIN_MAINTENANCE    16
#> 9                              REACTOME_CHOLESTEROL_BIOSYNTHESIS    26
#> 10                                         REACTOME_MRNA_CAPPING    29
#> 11                       REACTOME_TRNA_PROCESSING_IN_THE_NUCLEUS    57
#> 12                                      REACTOME_PROTEIN_FOLDING    88
#> 13                          GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC    59
#> 14                              GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS   102
#> 15                      GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC   100
#> 16             REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE    24
#> 17 REACTOME_ACTIVATION_OF_NMDA_RECEPTORS_AND_POSTSYNAPTIC_EVENTS    84
#> 18                        REACTOME_TRAFFICKING_OF_AMPA_RECEPTORS    31
#> 19                 GOMF_VOLTAGE_GATED_POTASSIUM_CHANNEL_ACTIVITY    86
#> 20                                GOBP_OXIDATIVE_PHOSPHORYLATION   135
#>        score     evidence
#> 1  0.6973029 1.0000000000
#> 2  0.6160519 0.9992882060
#> 3  0.6215357 0.9987039153
#> 4  0.6868095 0.9983546889
#> 5  0.7101165 0.9955955086
#> 6  0.5676972 0.9722121547
#> 7  0.6292330 0.9696651274
#> 8  0.6594246 0.9580399610
#> 9  0.6274162 0.9524832882
#> 10 0.5795515 0.8521399667
#> 11 0.5591023 0.8510963319
#> 12 0.5025770 0.5249499053
#> 13 0.3373526 0.0001317621
#> 14 0.2380503 0.0000000000
#> 15 0.2670079 0.0000000000
#> 16 0.1769775 0.0000000000
#> 17 0.3238279 0.0000000000
#> 18 0.2031705 0.0000000000
#> 19 0.3206181 0.0000000000
#> 20 0.2656734 0.0000000000
```

The `overall` table compares each region with the spot-weighted average of all
other regions. For the WM row, `score` is the posterior mean of π, how often a
pathway gene's expression state in WM exceeds its state in the cortical
layers, and `evidence` is the posterior probability that π > 0.5. Evidence
near 1 indicates higher activity in WM and evidence near 0 higher activity in
the cortical layers. `res$pairwise` gives the same quantities for every pair
of regions. See [API](../api.md#compare_regions) for all outputs.
