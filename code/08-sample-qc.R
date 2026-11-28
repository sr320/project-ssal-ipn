#!/usr/bin/env Rscript
# Sample-level QC on windowed CpG methylation. Called by code/08-sample-qc.sh:
#   Rscript code/08-sample-qc.R samples.tsv tiles_dir out_dir tile_bp
# Base R only, so it does not depend on which packages the container has.
#
# Writes to out_dir:
#   filter-counts.tsv        windows passing each minimum-read filter in all libraries
#   sample-summary.tsv       per library: windows covered, reads, methylation, PCs
#   chromosome-methylation.tsv   per library and chromosome
#   correlation.tsv          pairwise Pearson correlation of window methylation
#   pca-variance.tsv         variance explained per component
#   sample-qc.pdf            all plots; also one PNG per plot when the build allows

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: 08-sample-qc.R samples.tsv tiles_dir out_dir tile_bp")
samples_file <- args[1]; tiles_dir <- args[2]; out_dir <- args[3]
tile_bp <- as.integer(args[4])
min_reads <- as.integer(Sys.getenv("MIN_READS", "10"))   # filter used for correlation and PCA
min_cpgs  <- as.integer(Sys.getenv("MIN_CPGS", "3"))     # covered CpGs per window, per library
max_pca   <- as.integer(Sys.getenv("MAX_PCA_WINDOWS", "500000"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

ss <- read.delim(samples_file, stringsAsFactors = FALSE)
libs <- ss$sample
grp <- setNames(ss$group, libs)

# ---- read windows ---------------------------------------------------------
read_tiles <- function(s) {
  d <- read.delim(gzfile(file.path(tiles_dir, paste0(s, ".tsv.gz"))), header = FALSE,
                  col.names = c("chr", "start", "meth", "total", "cpgs"),
                  colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
  d$key <- paste(d$chr, format(d$start, scientific = FALSE, trim = TRUE), sep = ":")
  message(s, ": ", nrow(d), " windows")
  d
}
tl <- lapply(libs, read_tiles); names(tl) <- libs

keys <- unique(unlist(lapply(tl, `[[`, "key"), use.names = FALSE))
n <- length(keys)
M <- matrix(0, n, length(libs), dimnames = list(NULL, libs)); Tt <- M; N <- M
for (s in libs) {
  i <- match(tl[[s]]$key, keys)
  M[i, s] <- tl[[s]]$meth; Tt[i, s] <- tl[[s]]$total; N[i, s] <- tl[[s]]$cpgs
}
chr <- sub(":[0-9]+$", "", keys)
rm(tl); invisible(gc())
message(n, " windows of ", tile_bp, " bp with reads in at least one library")

# ---- how many windows survive each filter ---------------------------------
min_total <- apply(Tt, 1, min); min_n <- apply(N, 1, min)
grp_total <- sapply(unique(grp), function(g) rowSums(Tt[, grp == g, drop = FALSE]))
filt <- do.call(rbind, lapply(c(1, 3, 5, 10, 15, 20, 30), function(k) data.frame(
  min_reads_per_library = k,
  windows_all_libraries = sum(min_total >= k),
  windows_all_libraries_min_cpgs = sum(min_total >= k & min_n >= min_cpgs),
  pct_of_windows = round(100 * sum(min_total >= k) / n, 1))))
write.table(filt, file.path(out_dir, "filter-counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
grpfilt <- do.call(rbind, lapply(c(10, 20, 30, 50), function(k) data.frame(
  min_reads_per_group = k, windows_all_groups = sum(apply(grp_total, 1, min) >= k))))
write.table(grpfilt, file.path(out_dir, "filter-counts-group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

keep <- which(min_total >= min_reads & min_n >= min_cpgs)
message(length(keep), " windows with >= ", min_reads, " reads and >= ", min_cpgs, " CpGs in every library")
if (length(keep) < 1000) stop("too few windows pass the filter; lower MIN_READS")
B <- M[keep, , drop = FALSE] / Tt[keep, , drop = FALSE]

# ---- per-library summary ---------------------------------------------------
summ <- data.frame(
  sample = libs, group = grp[libs],
  windows_with_reads = colSums(Tt > 0),
  windows_ge_min_reads = colSums(Tt >= min_reads),
  reads_in_windows = colSums(Tt),
  mean_reads_per_window = round(colSums(Tt) / n, 2),
  methylation_all = round(100 * colSums(M) / colSums(Tt), 2),          # read-weighted, all windows
  methylation_shared = round(100 * colMeans(B), 2),                    # mean of shared windows
  median_window_shared = round(100 * apply(B, 2, median), 2),
  pct_windows_below_20 = round(100 * colMeans(B < 0.2), 2),            # unmethylated windows
  pct_windows_above_80 = round(100 * colMeans(B > 0.8), 2))

# ---- per-chromosome methylation (named chromosomes only) -------------------
is_chr <- grepl("^NC_", chr)
chrs <- sort(unique(chr[is_chr]))
if (length(chrs)) {
  cm <- rowsum(M[is_chr, , drop = FALSE], chr[is_chr]) / rowsum(Tt[is_chr, , drop = FALSE], chr[is_chr])
  cr <- rowsum(Tt[is_chr, , drop = FALSE], chr[is_chr])
  chr_tab <- data.frame(chromosome = rownames(cm), round(100 * cm, 2), check.names = FALSE)
  write.table(chr_tab, file.path(out_dir, "chromosome-methylation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # share of each library's reads per chromosome: flags sex or copy-number differences
  share <- sweep(cr, 2, colSums(cr), "/")
  write.table(data.frame(chromosome = rownames(share), round(100 * share, 3), check.names = FALSE),
              file.path(out_dir, "chromosome-read-share.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- correlation, clustering, PCA -----------------------------------------
cc <- cor(B)
write.table(data.frame(sample = libs, round(cc, 4), check.names = FALSE),
            file.path(out_dir, "correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summ$mean_cor_with_others <- round((rowSums(cc) - 1) / (length(libs) - 1), 4)

pi <- if (length(keep) > max_pca) sort(sample(length(keep), max_pca)) else seq_along(keep)
v <- apply(B[pi, , drop = FALSE], 1, var)
pc <- prcomp(t(B[pi[v > 0], , drop = FALSE]), center = TRUE, scale. = FALSE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
write.table(data.frame(component = paste0("PC", seq_along(ve)), pct_variance = round(100 * ve, 2)),
            file.path(out_dir, "pca-variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
npc <- min(4, ncol(pc$x))
summ <- cbind(summ, round(pc$x[libs, seq_len(npc), drop = FALSE], 3))
write.table(summ, file.path(out_dir, "sample-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Do the main components follow group, or depth / overall methylation (technical)?
assoc <- do.call(rbind, lapply(seq_len(npc), function(k) {
  x <- pc$x[libs, k]
  data.frame(component = paste0("PC", k), pct_variance = round(100 * ve[k], 2),
             group_r2 = round(summary(lm(x ~ factor(grp[libs])))$r.squared, 3),
             group_p = signif(anova(lm(x ~ factor(grp[libs])))[1, "Pr(>F)"], 3),
             cor_depth = round(cor(x, summ$reads_in_windows), 3),
             cor_methylation = round(cor(x, summ$methylation_all), 3))
}))
write.table(assoc, file.path(out_dir, "pca-associations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- plots -----------------------------------------------------------------
gcol <- setNames(c("#3b6fb6", "#d1803a", "#4a9d6f", "#8a5aa8", "#b84a4a", "#666666")[seq_along(unique(grp))], unique(grp))
cols <- gcol[grp[libs]]
plots <- list(
  "methylation-by-library" = function() {
    par(mar = c(6, 4.5, 3, 1))
    barplot(summ$methylation_all, names.arg = libs, las = 2, col = cols, border = NA,
            ylim = range(summ$methylation_all) + c(-1.5, 1), xpd = FALSE,
            ylab = "CpG methylation, all reads (%)", main = "Overall CpG methylation per library")
    legend("topright", names(gcol), fill = gcol, border = NA, bty = "n", horiz = TRUE)
  },
  "window-methylation-density" = function() {
    par(mar = c(4.5, 4.5, 3, 1))
    ds <- lapply(libs, function(s) density(100 * B[, s], from = 0, to = 100, n = 256))
    plot(NA, xlim = c(0, 100), ylim = c(0, max(sapply(ds, function(d) max(d$y)))),
         xlab = paste0("Window methylation (%), ", tile_bp, " bp windows"), ylab = "Density",
         main = paste0("Distribution over ", format(length(keep), big.mark = ","), " shared windows"))
    for (i in seq_along(libs)) lines(ds[[i]], col = cols[i])
    legend("topleft", names(gcol), col = gcol, lwd = 2, bty = "n")
  },
  "correlation-heatmap" = function() {
    o <- hclust(as.dist(1 - cc), "average")$order
    z <- cc[o, o]; par(mar = c(6, 6, 3, 1))
    image(seq_along(libs), seq_along(libs), z, axes = FALSE, xlab = "", ylab = "",
          col = hcl.colors(50, "YlGnBu", rev = TRUE), zlim = range(z[z < 1]),
          main = sprintf("Pairwise correlation of window methylation (%.3f to %.3f)", min(z), max(z[z < 1])))
    axis(1, seq_along(libs), libs[o], las = 2, cex.axis = 0.8)
    axis(2, seq_along(libs), libs[o], las = 2, cex.axis = 0.8)
    text(rep(seq_along(libs), each = length(libs)), rep(seq_along(libs), length(libs)),
         sub("^0", "", sprintf("%.2f", t(z))), cex = 0.5)
  },
  "clustering" = function() {
    par(mar = c(2, 4.5, 3, 1))
    h <- hclust(as.dist(1 - cc), "average")
    plot(h, labels = libs, main = "Clustering on 1 - correlation (average linkage)", xlab = "", sub = "", ylab = "1 - r")
  },
  "pca" = function() {
    par(mar = c(4.5, 4.5, 3, 1))
    plot(pc$x[, 1], pc$x[, 2], col = cols, pch = 19, cex = 1.6,
         xlab = sprintf("PC1 (%.1f%%)", 100 * ve[1]), ylab = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
         main = paste0("PCA of window methylation, ", format(sum(v > 0), big.mark = ","), " windows"))
    text(pc$x[, 1], pc$x[, 2], libs, pos = 3, cex = 0.7, xpd = NA)
    legend("topright", names(gcol), col = gcol, pch = 19, bty = "n")
  },
  "pca-pc3-pc4" = function() {
    if (ncol(pc$x) < 4) return(invisible())
    par(mar = c(4.5, 4.5, 3, 1))
    plot(pc$x[, 3], pc$x[, 4], col = cols, pch = 19, cex = 1.6,
         xlab = sprintf("PC3 (%.1f%%)", 100 * ve[3]), ylab = sprintf("PC4 (%.1f%%)", 100 * ve[4]), main = "PC3 and PC4")
    text(pc$x[, 3], pc$x[, 4], libs, pos = 3, cex = 0.7, xpd = NA)
  },
  "chromosome-methylation" = function() {
    if (!length(chrs)) return(invisible())
    par(mar = c(7, 4.5, 3, 1))
    d <- 100 * sweep(cm, 1, rowMeans(cm))   # each library's difference from the chromosome mean
    matplot(seq_along(chrs), d, type = "l", lty = 1, col = cols, axes = FALSE,
            xlab = "", ylab = "Difference from chromosome mean (points)",
            main = "Per-chromosome methylation, each library against the mean of all")
    axis(1, seq_along(chrs), chrs, las = 2, cex.axis = 0.6); axis(2); box(); abline(h = 0, lty = 3)
    lo <- which.min(colMeans(d)); text(length(chrs), d[length(chrs), lo], libs[lo], pos = 2, cex = 0.7)
    legend("topleft", names(gcol), col = gcol, lwd = 2, bty = "n", horiz = TRUE)
  })

pdf(file.path(out_dir, "sample-qc.pdf"), width = 9, height = 7)
for (f in plots) f()
invisible(dev.off())
if (isTRUE(capabilities("cairo")) || isTRUE(capabilities("png"))) {
  for (nm in names(plots)) {
    ok <- tryCatch({
      png(file.path(out_dir, paste0(nm, ".png")), width = 1350, height = 1050, res = 150,
          type = if (isTRUE(capabilities("cairo"))) "cairo" else getOption("bitmapType"))
      plots[[nm]](); invisible(dev.off()); TRUE
    }, error = function(e) { message("no PNG for ", nm, ": ", conditionMessage(e)); FALSE })
  }
}

cat("\nWindows passing each filter in all libraries\n"); print(filt, row.names = FALSE)
cat("\nPer-library summary\n"); print(summ, row.names = FALSE)
cat("\nComponents against group, depth and overall methylation\n"); print(assoc, row.names = FALSE)
