setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import motif counts ----
dat <- readRDS("db/motif_counts/twist008_motif_counts_selected.rds")
dat <- dat[grepl("^dev", ID)]
mat <- as.matrix(dat, 1)

# Cluster enhancers based on their motifs ----
set.seed(1)
cl <- vl_heatmap(log2(mat+1),
                 kmeans.k = 4)

# Motif enrichment compared to negative control sequences ----
mot <- split(dat, cl$rows$cl)
mot[["ctl"]] <- readRDS("db/motif_counts/twist008_motif_counts_selected_negative_controls.rds")
mot <- lapply(mot, as.matrix, 1)
enr <- vl_motif_cl_enrich(mot, control.cl = "ctl")

# Select top motifs ----
setorderv(enr, "padj")
sel <- enr[variable %in% enr[, .SD[1:3], cl]$variable]
res <- dcast(sel,
             name~cl,
             value.var = "log2OR",
             fun.aggregate = mean)
mat <- as.matrix(res, 1)
colnames(mat) <- paste0("Cluster ", colnames(mat))

# Import full dataset ----
dat <- readRDS("db/linear_models/FC_DSCP_large_WT_lm_predictions.rds")
dat[cl$rows, clL:= i.cl, on= "L==name"]
dat[cl$rows, clR:= i.cl, on= "R==name"]
dat <- dat[!is.na(clL) & !is.na(clR)]
dat[, class:= ifelse(clL==clR, "Homotypic", "Heterotypic")]

# Plot ----
pdf("pdf/draft/reviewer_figure_residuals_homotypic_hetero_pairs.pdf", 7.5, 4)
vl_par(mai= c(.9, 1.5, .9, 1.5),
       mfrow= c(1,2))
vl_heatmap(mat,
           breaks = c(-3, 0, 3),
           cluster.cols = F,
           legend.title = "Odd ratio (log2)",
           tilt.colnames = T,
           legend.cex = 7/12)
vl_par(mai= c(1.4, .9, 1.4, 1.2))
plot(density(dat[class=="Homotypic", residuals]),
     xlim= c(-4, 4),
     xlab= "Residuals",
     main= NA)
lines(density(dat[class=="Heterotypic", residuals]),
      col= "red")
vl_legend(legend= c("Homotypic", "Heterotypic"),
          col= c("black", "red"),
          lwd= 1)
dev.off()