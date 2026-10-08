# Knit the example R Markdown files into the documentation pages.
#
# Run from the repository root:
#   Rscript doc/examples/render.R                # both examples
#   Rscript doc/examples/render.R observed       # one example
#
# Each <page>.Rmd is knitted to doc/source/examples/<page>.md, with figures in
# doc/source/examples/figures/. Code runs in a work folder outside the
# repository (SPAPATH_DOC_WORK, default ~/spapath_doc_work), which also holds
# the knitr cache and the result folders written by the examples. The DLPFC
# GMT file is copied there under the same relative path, since the examples
# read it from benchmark/analysis/DLPFC/data/.

library(knitr)
knit_hooks$set(crop = hook_pdfcrop)   # crop = TRUE trims the white margins of a figure

pages <- commandArgs(trailingOnly = TRUE)
if (!length(pages)) pages <- c("observed", "composition")

source_dir <- normalizePath("doc/examples")
page_dir   <- normalizePath("doc/source/examples")
work_dir   <- path.expand(Sys.getenv("SPAPATH_DOC_WORK", "~/spapath_doc_work"))
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
work_dir   <- normalizePath(work_dir)
gmt_dir    <- "benchmark/analysis/DLPFC/data"
dir.create(file.path(work_dir, gmt_dir), showWarnings = FALSE, recursive = TRUE)
invisible(file.copy(file.path(gmt_dir, "DLPFC_brain_and_reference.gmt"),
                     file.path(work_dir, gmt_dir), overwrite = TRUE))

for (page in pages) {
  opts_knit$set(base.dir = page_dir, root.dir = work_dir)
  opts_chunk$set(fig.path = paste0("figures/", page, "-"),
                 cache.path = file.path(work_dir, "cache", page, ""),
                 comment = "#>", fig.width = 6, fig.height = 6, dpi = 100,
                 error = FALSE)   # stop at the first error rather than knit it into the page
  knit(file.path(source_dir, paste0(page, ".Rmd")),
       file.path(page_dir, paste0(page, ".md")),
       envir = new.env())
}
