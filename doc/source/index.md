# PISTE

PISTE tests whether a gene pathway is more active in one tissue region than
in another in spatial transcriptomics data. It models all genes of a pathway
jointly and, for every pair of regions and for every region against the rest
of the tissue, reports a directional score and its posterior evidence.
Adjacent sections of the same tissue can be fitted jointly, so that the score
reflects the regional difference reproduced across sections.

PISTE is distributed as an R package.

![Overview of PISTE](figures/method_overview.png)

- [Installation](installation.md)
- [Command line](cli.md)
- [API](api.md)
- Examples: [observed mode](examples/observed.md) and
  [composition-adjusted mode](examples/composition.md)
- [About](about.md)
