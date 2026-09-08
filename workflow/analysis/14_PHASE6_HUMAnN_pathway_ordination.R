ph6 <- "/path/to/PFAS_mSystems_revision/05_humann_functional_profiles"
amat <- file.path(ph6, "04_analysis_matrices")
ord <- file.path(ph6, "06_ordination")
dir.create(ord, recursive=TRUE, showWarnings=FALSE)

infile <- file.path(amat, "Ronneby_HUMAnN_pathabundance_cpm_unstratified_no_unmapped_prevalence10.tsv")

x <- read.table(infile, header=TRUE, sep="\t", check.names=FALSE, quote="", comment.char="")
features <- x[[1]]
mat <- as.matrix(x[, -1, drop=FALSE])
rownames(mat) <- features
storage.mode(mat) <- "numeric"

# Remove zero-variance rows.
vars <- apply(log1p(mat), 1, var)
mat <- mat[is.finite(vars) & vars > 0, , drop=FALSE]
features <- rownames(mat)

# Top-variable pathways.
vars <- apply(log1p(mat), 1, var)
top <- data.frame(
  Pathway=names(sort(vars, decreasing=TRUE)),
  variance_log1p=as.numeric(sort(vars, decreasing=TRUE)),
  row.names=NULL
)
write.table(top, file.path(ord, "Ronneby_HUMAnN_top_variable_pathways.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# PCA on log1p CPM, samples as rows.
sample_mat <- t(log1p(mat))
pca <- prcomp(sample_mat, center=TRUE, scale.=TRUE)
ve <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- data.frame(
  sample_id=rownames(pca$x),
  PC1=pca$x[, 1],
  PC2=pca$x[, 2],
  PC3=pca$x[, 3],
  row.names=NULL
)
write.table(scores, file.path(ord, "Ronneby_HUMAnN_pathway_PCA_sample_scores.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

loads <- data.frame(
  Pathway=rownames(pca$rotation),
  PC1=pca$rotation[, 1],
  PC2=pca$rotation[, 2],
  PC3=pca$rotation[, 3],
  abs_PC1=abs(pca$rotation[, 1]),
  abs_PC2=abs(pca$rotation[, 2]),
  row.names=NULL
)
loads <- loads[order(-pmax(loads$abs_PC1, loads$abs_PC2)), ]
write.table(loads, file.path(ord, "Ronneby_HUMAnN_pathway_PCA_feature_loadings.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Bray-Curtis PCoA from CPM matrix.
raw <- t(mat)
n <- nrow(raw)
d <- matrix(0, n, n)
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    den <- sum(raw[i, ] + raw[j, ])
    if (den > 0) {
      d[i, j] <- sum(abs(raw[i, ] - raw[j, ])) / den
    } else {
      d[i, j] <- 0
    }
  }
}
rownames(d) <- rownames(raw)
colnames(d) <- rownames(raw)

pcoa <- cmdscale(as.dist(d), eig=TRUE, k=3)
eig <- pcoa$eig
pos_eig <- eig[eig > 0]
pcoa_ve <- pos_eig / sum(pos_eig)

pcoa_scores <- data.frame(
  sample_id=rownames(raw),
  PCoA1=pcoa$points[, 1],
  PCoA2=pcoa$points[, 2],
  PCoA3=pcoa$points[, 3],
  row.names=NULL
)
write.table(pcoa_scores, file.path(ord, "Ronneby_HUMAnN_pathway_Bray_PCoA_sample_scores.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

summary <- data.frame(
  metric=c(
    "input_file",
    "n_samples",
    "n_pathways_after_prevalence_filter",
    "PCA_PC1_variance",
    "PCA_PC2_variance",
    "PCA_PC3_variance",
    "Bray_PCoA1_positive_eigen_fraction",
    "Bray_PCoA2_positive_eigen_fraction",
    "Bray_PCoA3_positive_eigen_fraction"
  ),
  value=c(
    infile,
    ncol(mat),
    nrow(mat),
    ve[1],
    ve[2],
    ve[3],
    pcoa_ve[1],
    pcoa_ve[2],
    pcoa_ve[3]
  )
)
write.table(summary, file.path(ord, "Ronneby_HUMAnN_pathway_ordination_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

pdf(file.path(ord, "Figure_S11_Ronneby_HUMAnN_pathway_PCA.pdf"), width=6, height=5)
plot(scores$PC1, scores$PC2,
     xlab=paste0("PC1 (", round(100*ve[1], 1), "%)"),
     ylab=paste0("PC2 (", round(100*ve[2], 1), "%)"),
     main="Ronneby HUMAnN pathway PCA")
text(scores$PC1, scores$PC2, labels=scores$sample_id, cex=0.5, pos=3)
dev.off()

pdf(file.path(ord, "Figure_S11_Ronneby_HUMAnN_pathway_Bray_PCoA.pdf"), width=6, height=5)
plot(pcoa_scores$PCoA1, pcoa_scores$PCoA2,
     xlab=paste0("PCoA1 (", round(100*pcoa_ve[1], 1), "% positive eig.)"),
     ylab=paste0("PCoA2 (", round(100*pcoa_ve[2], 1), "% positive eig.)"),
     main="Ronneby HUMAnN pathway Bray-Curtis PCoA")
text(pcoa_scores$PCoA1, pcoa_scores$PCoA2, labels=pcoa_scores$sample_id, cex=0.5, pos=3)
dev.off()

cat("DONE\n")
cat(file.path(ord, "Ronneby_HUMAnN_pathway_ordination_summary.tsv"), "\n")
