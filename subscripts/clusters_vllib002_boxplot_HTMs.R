setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Get ChIP signal
vars <- c("ATAC", "H3K27Ac", "H3K4me1", "H3K4me3")
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
feat <- feat[, c("ID", vars), with= F]
cl <- readRDS("Rdata/clustering_lm_residuals_vllib002.rds")
rows <- feat[cl$rows, on= "ID==name"]
rows <- melt(rows, id.vars = c("ID", "cl"), measure.vars = vars)
rows[, cl:= factor(cl, levels= c("No syn.", "Weak syn.", "Medium syn.", "High syn."))]
cols <- feat[cl$cols, on= "ID==name"]
cols <- melt(cols, id.vars = c("ID", "cl"), measure.vars = vars)

Cc <- adjustcolor(c("grey90", "grey60", "limegreen", "tomato"), 0.6)
.c <- substitute(vl_boxplot(value~cl, 
                            tilt.names= T,
                            col= Cc[.GRP], 
                            ylab= ifelse(.GRP==1, "Enrichment", NA),
                            main= variable, 
                            # compute_pval= list(c(1,2), c(2,3), c(3,4)))
                            compute_pval= list(c(1,2), c(1,3), c(1,4))))

pdf("pdf/draft/clusters_vllib002_boxplot_HTMs.pdf",
    width= 3.8,
    height= 2.5)
par(mfrow=c(1,4),
    mgp= c(2, 0.5, 0),
    oma= c(0,3,0,0),
    mar= c(4.75,1,1.5,0.5),
    tcl= -0.2,
    las= 1,
    xpd= NA)
rows[, {eval(.c, .SD); print("")}, variable]
cols[, {eval(.c, .SD); print("")}, variable]
dev.off()