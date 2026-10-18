# DLPFC 12 sections: WM (target) versus all other labelled spots (rest),
# SpaPath / PathwayLGM 0.7.0 (lingzhuu/pathwayLGM main, src/spapath),
# composition-adjusted mode (negative binomial on raw counts, measurement
# error kappa2 in pi = harmonic mean of the gene kappa2_g), 30 gene sets of
# benchmark/analysis/DLPFC/data/DLPFC12_gene_sets.gmt (the sets of Figure X).
# Inputs of each section:
#   * raw UMI counts of <section>_filtered_feature_bc_matrix.h5; gene-set
#     symbols matched through genes.csv (first Ensembl ID of a repeated symbol)
#   * spots with an author layer label that RCTD fitted
#   * RCTD full-mode proportions of the 7 LIBD types
#     (RCTD_LIBD/<section>/<section>_spot_celltype_proportions.csv) and the
#     reference profiles RCTD used (mean nucleus UMI fractions,
#     DLPFC_151673_RCTD_LIBD/reference/reference_profiles.rds)
#   * library size = total UMI of the spot over all genes
#   * DESeq2 dispersion scales estimated on all genes of the section
# Other settings as the observed run
# (SpaPath_WM_brain30_observed_0.6.1.1/scripts/run_wm_brain30.R):
#   * pixel coordinates shifted to 0 and divided by the larger range
#   * mesh max_edge .25, cutoff .08, offset (.1, .2); h = .3; intercept only
#   * baseline 'rest', pi = TRUE, delta .05 and level .95 (package defaults)
#   * seed 42 before each fit; INLA threads 2:1 (lab PC)
# Resumable: each gene set is written to <job>.partial and renamed to <job>
# only after it finished; the dispersion scales of a section are saved once
# and reused. screen_summary.rds keeps `joint` (CCD weights, regional means
# and kappa2 per configuration), so pi and evidence can be recomputed without
# refitting.
# The three myelin gene sets run first in all 12 sections, then the rest.
# Usage: Rscript run_wm_adjusted.R [SECTION ...]   (default: all 12)
Sys.setenv(OMP_NUM_THREADS = "2", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
Sys.setlocale("LC_CTYPE", "English_United States.utf8")
st2 <- "C:/Users/lxs888/Documents/ChatGPT/ST 2"
.libPaths(file.path(st2, "tools", c("R-library-spapath-0.7.0", "R-library-rctd",
                                    "R-library-gsva", "R-library",
                                    "R-library-spapath-0.6.1", "R-library-sct-v2")))
suppressPackageStartupMessages({library(Matrix); library(PathwayLGM)})
version <- "0.7.0"
stopifnot(as.character(packageVersion("PathwayLGM")) == version,
          packageVersion("INLA") == package_version("26.08.07"),
          startsWith(find.package("PathwayLGM"),
                     file.path(st2, "tools/R-library-spapath-0.7.0")))

data_root <- file.path(st2, "data/DLPFC_12_sections_20261007")
out_root <- file.path(data_root, "SpaPath_WM_DLPFC30_adjusted_0.7.0")
gmt_file <- file.path(out_root, "DLPFC12_gene_sets.gmt")
reference_file <- file.path(st2, "data/DLPFC_151673_RCTD_LIBD/reference/reference_profiles.rds")
types <- c("Astro", "EndoMural", "Excit", "Inhib", "Micro", "Oligo", "OPC")
h <- .3
cfg <- list(threads = "2:1", pi = TRUE, baseline_region = "rest",
            pi_delta = .05, pi_level = .95)
sections <- c("151507", "151508", "151509", "151510", "151669", "151670",
              "151671", "151672", "151673", "151674", "151675", "151676")
args <- commandArgs(TRUE)
if (length(args)) {stopifnot(all(args %in% sections)); sections <- args}

read_gmt <- function(path) {
  lines <- strsplit(readLines(path), "\t", fixed = TRUE)
  setNames(lapply(lines, function(x) unique(x[-(1:2)])), vapply(lines, `[`, "", 1))
}
pathways <- read_gmt(gmt_file)
stopifnot(length(pathways) == 30,
          unname(tools::md5sum(gmt_file)) == "0a6c5ecd7ae54ad58d5da140fc0c667e")
profiles <- readRDS(reference_file)$profiles
stopifnot(all(types %in% colnames(profiles)))

# The run settings are fixed when the folder is first used; a later launch
# with different settings stops instead of mixing results.
manifest <- c(paste("PathwayLGM", version, "main; INLA", packageVersion("INLA")),
              paste("gmt md5", unname(tools::md5sum(gmt_file))),
              paste("reference md5", unname(tools::md5sum(reference_file))),
              "mode adjusted; raw UMI; RCTD full proportions (7 LIBD types); library size = total UMI; DESeq2 dispersion scales on all genes",
              "labelled spots fitted by RCTD; WM = target, rest = other layers",
              "mesh max_edge .25 cutoff .08 offset (.1, .2); h .3; X intercept; seed 42 per fit",
              paste("cfg:", paste(names(cfg), unlist(cfg), sep = "=", collapse = "; ")))
manifest_file <- file.path(out_root, "run_manifest.txt")
if (file.exists(manifest_file)) {
  if (!identical(readLines(manifest_file), manifest)) stop("run_manifest.txt differs from this script's settings; use a new output folder.")
} else writeLines(manifest, manifest_file)

log_line <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ...)
  cat(line, "\n"); flush.console()
  cat(line, "\n", file = file.path(out_root, "progress.log"), append = TRUE)
}

read_counts <- function(file) {
  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(h5$close_all())
  counts <- sparseMatrix(i = h5[["matrix/indices"]][] + 1, p = h5[["matrix/indptr"]][],
                         x = as.numeric(h5[["matrix/data"]][]), dims = h5[["matrix/shape"]][])
  dimnames(counts) <- list(h5[["matrix/features/id"]][], h5[["matrix/barcodes"]][])
  counts
}

jobs <- sprintf("P%02d", seq_along(pathways))
first <- jobs[names(pathways) %in% c("GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM",
                                     "GOBP_MYELIN_MAINTENANCE", "GOCC_MYELIN_SHEATH")]
stopifnot(length(first) == 3)
for (pass_jobs in list(first, jobs)) for (section in sections) {
  out <- file.path(out_root, section)
  if (all(dir.exists(file.path(out, pass_jobs)))) next
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  sec_dir <- file.path(data_root, section)

  counts <- read_counts(file.path(sec_dir, paste0(section, "_filtered_feature_bc_matrix.h5")))
  stopifnot(!anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)),
            all(counts@x == round(counts@x)))
  genes_csv <- read.csv(file.path(sec_dir, "genes.csv"), stringsAsFactors = FALSE)
  stopifnot(setequal(genes_csv$gene_id, rownames(counts)))
  genes_csv <- genes_csv[!duplicated(genes_csv$symbol), ]
  symbol_id <- setNames(genes_csv$gene_id, genes_csv$symbol)

  rctd <- read.csv(file.path(data_root, "RCTD_LIBD", section,
                             paste0(section, "_spot_celltype_proportions.csv")),
                   stringsAsFactors = FALSE)
  spots <- rctd[!is.na(rctd$layer) & nzchar(rctd$layer) & rctd$fit_status == "fitted", ]
  region <- ifelse(spots$layer == "WM", "target", "rest")
  stopifnot(sum(region == "target") >= 20, sum(region == "rest") >= 20,
            all(spots$spot %in% colnames(counts)))
  proportions <- as.matrix(spots[, types])
  rownames(proportions) <- spots$spot

  Y_all <- t(counts[, spots$spot])                  # spots x genes, sparse
  rm(counts); gc(FALSE)
  library_size <- rowSums(Y_all)
  dispersion_file <- file.path(out, "dispersion_scale.csv")
  if (!file.exists(dispersion_file)) {
    log_line(section, "estimating dispersion scales on", ncol(Y_all), "genes")
    dispersion <- estimate_dispersion_scale(Y_all, library_size)
    write.csv(dispersion, paste0(dispersion_file, ".partial"), row.names = FALSE)
    stopifnot(file.rename(paste0(dispersion_file, ".partial"), dispersion_file))
  }
  dispersion <- read.csv(dispersion_file, stringsAsFactors = FALSE)
  scale <- setNames(dispersion$dispersion_scale, dispersion$gene)

  # Gene-set genes by symbol: present in the counts and the reference, with a
  # dispersion scale and a positive predicted expression in every spot.
  panel <- unique(unlist(pathways, use.names = FALSE))
  ids <- symbol_id[intersect(panel, names(symbol_id))]
  ids <- ids[ids %in% rownames(profiles) & !is.na(scale[ids])]
  log_m <- composition_offset(proportions, t(profiles[ids, types, drop = FALSE]))
  colnames(log_m) <- names(ids)
  ids <- ids[colSums(!is.finite(log_m)) == 0]
  log_m <- log_m[, names(ids), drop = FALSE]
  Y <- as.matrix(Y_all[, ids, drop = FALSE])
  colnames(Y) <- names(ids)
  storage.mode(Y) <- "integer"
  gene_scale <- setNames(scale[ids], names(ids))
  sets <- lapply(pathways, intersect, names(ids))
  rm(Y_all); gc(FALSE)
  stopifnot(all(lengths(sets) >= 2))

  xy <- as.matrix(spots[, c("pixel_col", "pixel_row")])
  xy <- sweep(xy, 2, apply(xy, 2, min))
  xy <- xy / max(xy)
  dimnames(xy) <- list(spots$spot, c("x", "y"))
  mesh <- build_mesh(xy, max_edge = .25, cutoff = .08, offset = c(.1, .2))
  X <- matrix(1, nrow(xy), 1, dimnames = list(rownames(xy), "intercept"))

  write.csv(data.frame(job = jobs, pathway = names(sets),
                       panel_genes = lengths(pathways), matched_genes = lengths(sets)),
            file.path(out, "gene_counts.csv"), row.names = FALSE)
  write.csv(data.frame(spot = spots$spot, layer = spots$layer, region = region,
                       library_size = library_size, proportions, check.names = FALSE),
            file.path(out, "spots_used.csv"), row.names = FALSE)
  writeLines(c(paste("Section", section, "labelled RCTD-fitted spots", nrow(spots),
                     "WM", sum(region == "target"), "rest", sum(region == "rest")),
               paste("Rest layers:", paste(sort(unique(spots$layer[region == "rest"])), collapse = "; ")),
               paste("Usable gene-set genes:", length(ids), "of", length(panel)),
               paste("Mesh nodes:", mesh$n)),
             file.path(out, "parameters.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

  status_file <- file.path(out, "status.csv")
  for (i in seq_along(sets)) {
    job <- jobs[i]; p <- names(sets)[i]; genes <- sets[[p]]
    final <- file.path(out, job)
    if (!job %in% pass_jobs || dir.exists(final)) next
    partial <- paste0(final, ".partial")
    for (attempt in 1:2) {
      unlink(partial, recursive = TRUE)
      start <- proc.time()[["elapsed"]]
      error <- ""; ok <- FALSE
      log_line("Starting", section, job, p, "genes", length(genes), "attempt", attempt)
      tryCatch({
        set.seed(42)
        fit <- screen_joint_fit(Y[, genes, drop = FALSE], xy, region, mesh, h,
                                mode = "adjusted", X = X, library_size = library_size,
                                dispersion_scale = gene_scale[genes],
                                composition_offset = log_m[, genes, drop = FALSE],
                                cfg = cfg)
        summ <- screen_summarize(fit)
        stopifnot(nrow(summ$pi$overall) == 2,
                  setequal(summ$pi$overall$region, c("target", "rest")),
                  identical(summ$metadata$baseline_region, "rest"),
                  all(is.finite(summ$pi$overall$score)),
                  all(is.finite(summ$pi$overall$evidence)),
                  summ$metadata$model_version == version,
                  length(summ$joint$weights) >= 2)
        screen_write(summ, partial)
        stopifnot(any(grepl("converged: TRUE", readLines(file.path(partial, "complete.txt")))),
                  file.rename(partial, final))
        ok <- TRUE
      }, error = function(e) error <<- conditionMessage(e))
      secs <- proc.time()[["elapsed"]] - start
      tgt <- if (ok) summ$pi$overall[summ$pi$overall$region == "target", ] else NULL
      row <- data.frame(section = section, job = job, pathway = p,
        status = if (ok) "ok" else "error", attempt = attempt, n_genes = length(genes),
        target_spots = sum(region == "target"), rest_spots = sum(region == "rest"),
        pi_target = if (ok) tgt$score else NA, evidence_target = if (ok) tgt$evidence else NA,
        fit_seconds = secs, finished = format(Sys.time()), error = error)
      write.table(row, status_file, sep = ",", row.names = FALSE,
                  col.names = !file.exists(status_file), append = TRUE)
      log_line("Finished", section, job, if (ok) "ok" else paste("error:", error),
               "seconds", round(secs))
      for (nm in c("fit", "summ")) if (exists(nm, inherits = FALSE)) rm(list = nm)
      gc(FALSE)
      if (ok) break
    }
  }
  rm(Y, log_m); gc(FALSE)
  done <- sum(dir.exists(file.path(out, jobs)))
  log_line(section, "has", done, "of", length(jobs), "gene sets")
}
log_line("Run ended")
