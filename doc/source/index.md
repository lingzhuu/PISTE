# Welcome to PISTE's documentation

## Introduction

PISTE (**P**athway **I**nference from **S**patial **T**ranscriptomic
**E**xpression) tests whether a gene pathway is more active in one tissue
region than in another, and maps where in the tissue a pathway is more or
less active. It models all genes of a gene set jointly and is distributed as
an R package with a command-line tool.

## How to cite

A manuscript describing PISTE is in preparation.

## Key features

- **Region comparison**: a directional score π between 0 and 1 and its
  posterior evidence, for every pair of regions and for each region against
  the rest of the tissue.
- **Composition adjustment**: the regional difference that cell-type
  composition does not explain, from raw counts, cell-type proportions and a
  single-cell reference.
- **Adjacent sections**: sections of the same tissue fitted jointly, so that
  the score reflects the difference reproduced across sections.
- **Spatial map**: regions of higher and lower pathway activity from a spatial
  field shared by the genes of the set, without region labels.

## Overview of PISTE

![Overview of PISTE](figures/method_overview.png)
