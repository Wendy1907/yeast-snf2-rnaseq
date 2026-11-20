#!/usr/bin/env Rscript
# Differential expression of SNF2 knockout vs wild type with DESeq2.
# Usage: deseq2.R <gene_counts.tsv> <samplesheet.csv>
# Outputs (in the working directory): tables, a text summary and five PNG figures.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})
options(bitmapType = "cairo")
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: deseq2.R <gene_counts.tsv> <samplesheet.csv>")
counts_file <- args[1]
sheet_file  <- args[2]

# ---- Settings ---------------------------------------------------------------
control   <- "WT"       # reference level: log2 fold changes are SNF2 relative to WT
treatment <- "SNF2"
alpha     <- 0.05       # adjusted p-value threshold
lfc_min   <- 1          # |log2 fold change| threshold used for "strong" changes
knockout_gene <- "YOR290C"   # SNF2 itself: expected to be lost in the knockout
cols <- c(WT = "#3B6EA8", SNF2 = "#D95F02")

# ---- Read inputs ------------------------------------------------------------
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
sheet  <- read.csv(sheet_file, stringsAsFactors = FALSE)
stopifnot(all(sheet$condition %in% c(control, treatment)))
coldata <- data.frame(
  condition = factor(sheet$condition, levels = c(control, treatment)),
  row.names = sheet$sample
)
stopifnot(all(colnames(counts) %in% rownames(coldata)))
coldata <- coldata[colnames(counts), , drop = FALSE]

# ---- Differential expression ------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts), colData = coldata, design = ~ condition
)
# keep genes with at least 10 reads in at least 4 samples (the smallest group size)
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", treatment, control), alpha = alpha)
shrunk <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm", res = res)

out <- data.frame(
  gene_id           = rownames(res),
  baseMean          = signif(res$baseMean, 6),
  log2FoldChange    = signif(res$log2FoldChange, 6),
  log2FoldChange_shrunk = signif(shrunk$log2FoldChange, 6),
  pvalue            = signif(res$pvalue, 6),
  padj              = signif(res$padj, 6)
)
out <- out[order(out$padj, na.last = TRUE), ]
write.csv(out, "deseq2_results.csv", row.names = FALSE)
write.csv(data.frame(gene_id = rownames(dds),
                     signif(counts(dds, normalized = TRUE), 6), check.names = FALSE),
          "normalized_counts.csv", row.names = FALSE)

# ---- Text summary -----------------------------------------------------------
sig    <- !is.na(out$padj) & out$padj < alpha
strong <- sig & abs(out$log2FoldChange_shrunk) >= lfc_min
ko     <- out[out$gene_id == knockout_gene, ]
summary_lines <- c(
  sprintf("Comparison: %s vs %s (positive log2 fold change = higher in %s)", treatment, control, treatment),
  sprintf("Samples: %d %s, %d %s", sum(coldata$condition == control), control,
          sum(coldata$condition == treatment), treatment),
  sprintf("Genes tested after filtering: %d (of %d in the count table)", nrow(dds), nrow(counts)),
  sprintf("Significant (padj < %.2f): %d  [up: %d, down: %d]", alpha, sum(sig),
          sum(sig & out$log2FoldChange_shrunk > 0), sum(sig & out$log2FoldChange_shrunk < 0)),
  sprintf("Significant and |shrunken log2FC| >= %g: %d  [up: %d, down: %d]", lfc_min, sum(strong),
          sum(strong & out$log2FoldChange_shrunk > 0), sum(strong & out$log2FoldChange_shrunk < 0)),
  if (nrow(ko) == 1) sprintf("Knockout check %s: log2FC = %.2f, padj = %s", knockout_gene,
                             ko$log2FoldChange, format(ko$padj, digits = 3))
  else sprintf("Knockout check %s: gene not present after filtering", knockout_gene)
)
writeLines(summary_lines, "deseq2_summary.txt")

# ---- Figures ----------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
ann <- data.frame(condition = coldata$condition, row.names = rownames(coldata))
ann_cols <- list(condition = cols)

# 1. PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p1 <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 4) + geom_text_repel(show.legend = FALSE, size = 3.5) +
  scale_colour_manual(values = cols) +
  labs(title = "PCA of variance-stabilised counts",
       x = sprintf("PC1: %d%% variance", pv[1]), y = sprintf("PC2: %d%% variance", pv[2])) +
  theme_bw(base_size = 12)
ggsave("pca.png", p1, width = 6.5, height = 5, dpi = 150)

# 2. Sample-to-sample distances
dmat <- as.matrix(dist(t(assay(vsd))))
pheatmap(dmat, annotation_col = ann, annotation_row = ann, annotation_colors = ann_cols,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distances", filename = "sample_distance_heatmap.png",
         width = 7, height = 6)

# 3. MA plot
plot_df <- out[!is.na(out$padj), ]
plot_df$status <- ifelse(plot_df$padj < alpha, "padj < 0.05", "not significant")
p3 <- ggplot(plot_df, aes(baseMean, log2FoldChange_shrunk, colour = status)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_x_log10() + scale_colour_manual(values = c("padj < 0.05" = "#D95F02", "not significant" = "grey60")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "MA plot (shrunken log2 fold changes)", x = "Mean of normalised counts",
       y = sprintf("log2 fold change (%s vs %s)", treatment, control)) +
  theme_bw(base_size = 12)
ggsave("ma_plot.png", p3, width = 7, height = 5, dpi = 150)

# 4. Volcano plot
plot_df$neglog10 <- -log10(plot_df$padj)
plot_df$label <- ""
top <- head(plot_df[order(plot_df$padj), "gene_id"], 10)
plot_df$label[plot_df$gene_id %in% c(top, knockout_gene)] <-
  plot_df$gene_id[plot_df$gene_id %in% c(top, knockout_gene)]
plot_df$group <- ifelse(plot_df$padj < alpha & abs(plot_df$log2FoldChange_shrunk) >= lfc_min,
                        "padj < 0.05 and |log2FC| >= 1", "other")
p4 <- ggplot(plot_df, aes(log2FoldChange_shrunk, neglog10, colour = group, label = label)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_text_repel(colour = "black", size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c("padj < 0.05 and |log2FC| >= 1" = "#D95F02", "other" = "grey60")) +
  geom_vline(xintercept = c(-lfc_min, lfc_min), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  labs(title = "Volcano plot", x = sprintf("log2 fold change (%s vs %s)", treatment, control),
       y = "-log10 adjusted p-value") +
  theme_bw(base_size = 12)
ggsave("volcano.png", p4, width = 7, height = 5.5, dpi = 150)

# 5. Heatmap of the most significant genes
top_genes <- head(out$gene_id[!is.na(out$padj)], 30)
pheatmap(assay(vsd)[top_genes, ], scale = "row", annotation_col = ann, annotation_colors = ann_cols,
         main = "Top 30 genes by adjusted p-value (row-scaled)", filename = "top_genes_heatmap.png",
         width = 7, height = 8)

writeLines(capture.output(sessionInfo()), "deseq2_sessioninfo.txt")
