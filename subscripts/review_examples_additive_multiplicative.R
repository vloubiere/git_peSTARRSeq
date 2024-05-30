setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data ----
dev <- readRDS("db/linear_models/FC_DSCP_focused_lm_predictions.rds")[grepl("^dev", L) & grepl("^dev", R)]
hk <- readRDS("db/linear_models/FC_RpS12_focused_lm_predictions.rds")[grepl("^hk", L) & grepl("^hk", R)]
model <- readRDS("db/linear_models/lm_DSCP_large_WT.rds")
hk[, `Multiplicative model`:= predict(model, .SD)]
dev[, `Multiplicative model`:= predict(model, .SD)]

# Select individual enhancers with similar act within and between screens ----
value <- log2(16)
sel <- substitute(between(indL, value-.25, value+.25) & between(indR, value-.25, value+.25))
dev <- dev[eval(sel)]
hk <- hk[eval(sel)]
dev <- dev[order(abs(`Additive model`-value))]
hk <- hk[order(abs(`Additive model`-value))]

# Melt and plot ----
sel <- list(hk= hk[1:20],
            dev= dev[1:20])
sel <- rbindlist(sel,
                 idcol = "CP")
sel[, CP:= factor(CP, c("hk", "dev"))]
sel[, col:= c("tomato", "limegreen")[.GRP], CP]
setnames(sel,
         c("indL", "indR", "Additive model", "Multiplicative model", "log2FoldChange"),
         c("5' enh.", "3' enh.", "Pred. add.", "Pred. fit. mult.", "Combined"))
pl <- melt(sel,
           id.vars = c("CP", "L", "R", "col"),
           measure.vars = c("5' enh.", "3' enh.", "Pred. add.", "Pred. fit. mult.", "Combined"))
pl[, value:= 2^value]
pl <- pl[, .(mean= mean(value), sd= sd(value), value= .(value), col= .(col)), keyby= .(CP, variable)]

# Plot ----
pdf("pdf/draft/review_examples_add_mult_hkCP_dCP.pdf", 2.9, 2.7)
vl_par(mgp= c(1, .2, 0),
       mai= c(.9, .7, .9, .7),
       cex.axis= 7/12,
       cex.lab= 8/12,
       lwd= .75)
pl[, {
  bar <- vl_barplot(mean,
                    sd,
                    space= c(rep(0.1, 5), 2, rep(0.1, 4)),
                    individual.var = value,
                    ind.col = adjustcolor(unlist(col), .6),
                    ind.jitter = .4,
                    col= "lightgrey",
                    border= "lightgrey",
                    lend= 2,
                    ind.cex= .4,
                    ylab= "Activity",
                    density= c(NA, NA,50,50,NA))
  vl_tilt_xaxis(bar, labels = variable)
}]
dev.off()