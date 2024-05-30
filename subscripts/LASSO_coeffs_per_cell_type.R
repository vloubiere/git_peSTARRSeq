setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data ----
models <- list("S2 - ecd."= readRDS("db/lasso_models/full_dataset_large_WT_lasso.rds"),
               "S2 + ecd."= readRDS("db/lasso_models/full_dataset_ECD_lasso.rds"),
               "OSC cells"= readRDS("db/lasso_models/full_dataset_OSC_lasso.rds"))
dat <- lapply(models, function(x) as.data.table(as.matrix(x$log2FoldChange$beta), keep.rownames= T))
dat <- rbindlist(dat, idcol = "cdition")
dat[, cdition:= factor(cdition, 
                       c("S2 - ecd.",
                         "S2 + ecd.",
                         "OSC cells"))]

# Select top motifs ----
dat[, rank:= rank(-abs(s0)), cdition]
sel <- unique(dat[rank<=5, rn])
dat <- dat[rn %in% sel]

# Dcast Left and right motif ----
mat <- dcast(dat, rn~cdition, value.var = "s0")
mat <- as.matrix(mat, 1)

# Retrieve motif names ----
enr <- readRDS("db/motif_counts/lib8_motifs_enrichments.rds")
idx <- match(rownames(mat), enr$motif_ID)
rownames(mat) <- enr[idx, name]
pwms <- enr[idx, pwm]

# Plot ----
pdf("pdf/draft/LASSO_coeff_per_cell_type.pdf", 3.9, 3.25)
vl_par(mai= c(.9, 2, .9, 1.5),
       cex.main= 7/12,
       xpd= NA,
       bty= "o",
       lwd= .5)
hm <- vl_heatmap(mat,
                 cluster.cols = F,
                 breaks = c(-.5, 0, .5),
                 legend.title = "LASSO coef.",
                 tilt.colnames = T, 
                 legend.cex = 6/12,
                 box.lwd = .5)
vl_seqlogo(pwms,
           x= -5.5,
           y= hm$rows$y,
           pos = 2,
           add= T,
           cex.width = .5,
           min_content = .1)
dev.off()