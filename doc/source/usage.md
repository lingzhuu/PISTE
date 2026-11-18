# Quick start

Each analysis runs from the command line or from R. Words in angle brackets
are your own files and names. The [tutorials](examples/observed.md) run them on
real data.

## Region comparison

Compares a region with the rest of the section, on every gene set of a GMT
file.

=== "Command line"

    ```bash
    piste compare --input <normalized.h5ad> \
                  --region <region_column> --target <region> \
                  --gene_set <gene_sets.gmt> --out <output_folder>
    ```

=== "R"

    ```r
    library(PISTE)
    res <- compare(<expression>, <region>, "<gene_sets.gmt>",
                   out = "<output_folder>")
    ```

## Adjusted mode

Region comparison on raw counts, without the differences that come from
cell-type composition. It needs cell-type proportions and a reference.

=== "Command line"

    ```bash
    piste compare --mode adjusted --input <counts.h5ad> \
                  --region <region_column> --target <region> \
                  --proportions <proportions.csv> \
                  --reference <reference.h5ad> --celltype <celltype_column> \
                  --gene_set <gene_sets.gmt> --out <output_folder>
    ```

=== "R"

    ```r
    res <- compare(<counts>, <region>, "<gene_sets.gmt>", mode = "adjusted",
                   proportions = <proportions>, reference = <reference>,
                   out = "<output_folder>")
    ```

## Joint sections

Region comparison with adjacent sections of the same tissue fitted together.

=== "Command line"

    ```bash
    piste compare -joint_sections --input <normalized.h5ad> \
                  --region <region_column> --target <region> \
                  --section <section_column> \
                  --gene_set <gene_sets.gmt> --out <output_folder>
    ```

=== "R"

    ```r
    res <- compare(<expression>, <region>, "<gene_sets.gmt>",
                   section = <section>, joint_sections = TRUE,
                   out = "<output_folder>")
    ```

## Spatial map

Maps where each gene set is more or less active, without region labels.

=== "Command line"

    ```bash
    piste map --input <normalized.h5ad> \
              --coords <tissue_positions.csv> \
              --gene_set <gene_sets.gmt> --out <output_folder>
    ```

=== "R"

    ```r
    res <- map_activity(<expression>, <coords>, "<gene_sets.gmt>",
                        out = "<output_folder>")
    ```
