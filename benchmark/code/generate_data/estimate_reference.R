# Estimate NB reference parameters from real DLPFC 151673.
#
#   Y_gi ~ NB(mu_gi, theta_g),  mu_gi = c_i * m_g
#   Var(Y_gi) = mu_gi + mu_gi^2 / theta_g
#
# DESeq2 (intercept-only) gives spot depth factors c_i and per-gene
# dispersion = 1/theta_g; m_g is the depth-corrected mean expression.
#
# Input : ../../reference_151673/{genes.csv, spots.csv, counts.mtx}
# Output: ../../reference_151673/reference_parameters.csv  (per-gene m, theta)
#         ../../reference_151673/reference_depth.csv       (per-spot c_i)

library(Matrix)
library(DESeq2)

base <- "../../reference_151673"

genes <- read.csv(file.path(base, "genes.csv"))$gene
spots <- read.csv(file.path(base, "spots.csv"))
y <- as(readMM(file.path(base, "counts.mtx")), "dgCMatrix")
stopifnot(nrow(y) == length(genes), ncol(y) == nrow(spots))
rownames(y) <- genes
colnames(y) <- spots$spot

# Keep genes with enough counts for stable estimates.
keep <- rowSums(y) >= 10 & rowSums(y > 0) >= 3
mat <- as.matrix(y[keep, ])
storage.mode(mat) <- "integer"
rm(y)

dds <- DESeqDataSetFromMatrix(mat, data.frame(row.names = colnames(mat)), ~1)
rm(mat)
dds <- estimateSizeFactors(dds, type = "poscounts")
sizeFactors(dds) <- sizeFactors(dds) / exp(mean(log(sizeFactors(dds))))
dds <- estimateDispersions(dds, quiet = TRUE)

mat  <- counts(dds)
norm <- sweep(mat, 2, sizeFactors(dds), "/")   # depth-corrected counts
logy <- log1p(norm)
out <- data.frame(
  gene     = rownames(dds),
  m        = rowMeans(norm),
  theta    = 1 / dispersions(dds),
  a        = rowMeans(logy),
  sigma2   = apply(logy, 1, var),
  total    = rowSums(mat),
  detected = rowSums(mat > 0))

ok <- is.finite(out$m) & out$m > 0 &
      is.finite(out$theta) & out$theta > 0 &
      is.finite(out$a) & is.finite(out$sigma2) & out$sigma2 > 0
out <- out[ok, ]

write.csv(out, file.path(base, "reference_parameters.csv"), row.names = FALSE)
write.csv(data.frame(spot = colnames(dds), c_i = sizeFactors(dds)),
          file.path(base, "reference_depth.csv"), row.names = FALSE)
