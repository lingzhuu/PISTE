# SpaPath

SpaPath tests whether a gene pathway is more active in one tissue region than
in another in spatial transcriptomics data. It models all genes of a pathway
jointly, accounts for spatial dependence between neighbouring spots with
gene-specific spatial random fields, and fits the model with INLA. For every
pair of regions, and for every region against the rest of the tissue, it
reports a directional score and a posterior evidence.

SpaPath is distributed as the R package PathwayLGM.

- [Installation](installation.md)
- [Tutorials](tutorials.md)
- [Example usages](examples.md)
- [About](about.md)
