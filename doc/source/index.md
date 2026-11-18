# Overview

PISTE (**P**athway **I**nference from **S**patial **T**ranscriptomic
**E**xpression) finds the gene sets that are more or less active in a tissue
region than in the rest, and maps where in a section a gene set is more or
less active. It models all genes of a gene set together.

![Overview of PISTE](figures/method_overview.png){ .pd-figure }

<p class="pd-caption">Spatial data go into an expression model; a region comparison or a spatial map comes out.</p>

## Two analyses

- **Region comparison** needs region labels. For each gene set and region it
  reports a score π and the evidence that activity is higher or lower than in
  the other regions. Adjacent sections can be fitted jointly, and the
  adjusted mode removes differences that come from cell-type composition.
- **Spatial map** needs spot coordinates and no labels. It returns the regions
  of a section where a gene set is more or less active.

PISTE takes normalized expression, for example from SCTransform. Only the
adjusted mode takes raw counts.

## Three ways to run it

- **R package**: [`compare()`](api.md#compare) and
  [`map_activity()`](api.md#map_activity).
- **Command line**: [`piste compare` and `piste map`](cli.md).
- **[Web app](https://lingzhuu.github.io/PISTE/)**: region comparison of one
  section in the browser, with nothing to install.

## Start here

1. [Installation](installation.md)
2. [Quick start](usage.md): one command for each analysis.
3. [Region comparison](examples/observed.md): a full tutorial on a brain section.

A manuscript describing PISTE is in preparation.
