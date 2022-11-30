setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)
require(parallel)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/linear_models/FC_dev_pairs_vllib002_with_predictions.rds")

# Select variables of interest
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
lib <- lib[ID_vl %in% dat[, c(L, R)]]
sel <- c("flyfactorsurvey__CG16778_SANGER_5_FBgn0003715", "cisbp__M2014")
counts <- vl_motif_counts(lib$enh_sequence, sel= sel, genome= "dm3", collapse_overlapping = TRUE)
counts <- cbind(lib, counts)
mat <- dcast(counts, flyfactorsurvey__CG16778_SANGER_5_FBgn0003715~cisbp__M2014, fill= NA)
mat <- mat[nrow(mat):1]
vl_heatmap(as.matrix(mat, 1), 
           cluster_rows= F, 
           cluster_cols= F, 
           col= c("cornflowerblue", "yellow"),
           legend_title = "Motif counts", 
           display_numbers = T)

# plotting object
dat[counts, c("twist_L", "Trl_L"):= .(i.flyfactorsurvey__CG16778_SANGER_5_FBgn0003715, i.cisbp__M2014), on= "L==ID_vl"]
dat[counts, c("twist_R", "Trl_R"):= .(i.flyfactorsurvey__CG16778_SANGER_5_FBgn0003715, i.cisbp__M2014), on= "R==ID_vl"]

enr <- dat[twist_L>=2 & Trl_L>=2 & twist_R>=2 & Trl_R>=2]
ctls <- dat[twist_L+Trl_L+twist_R+Trl_R==0]
ctls <- enr[, {
  idx <- seq(.N)
  res <- .SD[, {
    ctls[, dist:= sqrt((cL-indL)^2+(cR-indR)^2)]
    ctls <- ctls[order(dist)]
    .c <- ctls[1]
    ctls <- ctls[-1]
    .c
  }, .(cL= indL, cR= indR, idx)]
}]
res <- list(ctls$indL,
            enr$indL,
            ctls$indR,
            enr$indR,
            ctls$log2FoldChange,
            enr$log2FoldChange)

xpos <- seq(1, 5, length.out= 3)
vl_boxplot(res,
           at= rep(xpos, each= 2)+c(-0.35, 0.35),
           tilt.names= T, 
           compute_pval= list(c(1,2), c(3,4), c(5,6)),
           main= paste0("Trl x Twist (", 
                        length(unique(enr$L)), " x 5'; ",
                        length(unique(enr$R)), " x 3'; ",
                        nrow(enr), " pairs)"),
           xaxt= "n",
           col= c("lightgrey", "rosybrown1"),
           ylab= "Activity (log2)")
axis(1, at= xpos, labels = c("5'", "3'", "Pair"))
legend(par("usr")[2],
       par("usr")[4],
       fill= c("lightgrey", "rosybrown1"),
       legend= c("Controls (no motifs)", 
                 paste0("Motif pairs (>= ", 2, ")")),
       bty= "n")
