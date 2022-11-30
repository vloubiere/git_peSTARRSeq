setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/linear_models/FC_vllib002_with_predictions.rds")

rsqAdd <- vl_model_eval(dat$log2FoldChange, dat$additive)$Rsquare
adj.rsqAdd <- 1-(((1-rsqAdd)*(nrow(dat)-1))/(nrow(dat)-2-1))

rsqMult <- vl_model_eval(dat$log2FoldChange, dat$multiplicative)$Rsquare
adj.rsqMult <- 1-(((1-rsqMult)*(nrow(dat)-1))/(nrow(dat)-2-1))

model <- readRDS("db/linear_models/lm_vllib002.rds")
adj.rsqlm <- summary(model)$adj.r.square

pdf("pdf/draft/Compare_add_mult_vllib002.pdf",
    height = 3,
    width = 1.5)
par(las= 1,
    mar= c(4.5,3,1,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    cex= 1)
bar <- barplot(c(adj.rsqAdd, adj.rsqMult, adj.rsqlm),
        ylab= "Adj. R2 (goodness of fit)")
vl_tilt_xaxis(bar, 
              labels = c("Additive", "Multiplicative", "Linear model"),)
dev.off()
