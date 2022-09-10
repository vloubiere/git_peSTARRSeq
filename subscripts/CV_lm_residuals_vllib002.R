setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(glmnet)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
cl <- readRDS("Rdata/clustering_lm_residuals_vllib002.rds")
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & L %in% cl$rows$name & R %in% cl$cols$name]
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")
lib[, predicted:= predict(model, newdata = .SD)]
lib[, residuals:= log2FoldChange-predicted]
set.seed(1)
dat <- lib[sample(nrow(lib), nrow(lib))]

# Define train and test sets
set.seed(1)
dat[dat[, .(set= sample(5)), L], setL:= i.set, on= "L"]
set.seed(1)
dat[dat[, .(set= sample(5)), R], setR:= i.set, on= "R"]
dat[, set:= paste0("test", .GRP), .(setL, setR)]

# Get motif counts
feat <- fread("Rdata/final_300bp_enhancer_features.txt")
motL <- feat[dat$L, names(feat) %in% vl_Dmel_motifs_DB_full$motif, with= F, on= "ID"]
setnames(motL, function(x) paste0(x, "_L"))
motR <- feat[dat$R, names(feat) %in% vl_Dmel_motifs_DB_full$motif, with= F, on= "ID"]
setnames(motR, function(x) paste0(x, "_R"))
dat <- cbind(dat, motL, motR)

# perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(as.matrix(dat[, c(names(motL), names(motR)), with= F]), 
                      dat$residuals, 
                      alpha = 1)

# find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda

[1] 5.616345



# Train linear model for each train set and compute predicted values
.f <- paste0("`", c(names(motL), names(motR)), "`")
.f <- paste0(c("log2FoldChange~median_L*median_R", .f), collapse= "+")

model <- lm(formula = as.formula(.f),
            data= dat)
model$CV_rsqs <- dat[, {
  print(set)
  cL <- L
  cR <- R
  train <- dat[!(L %in% cL) & !(R %in% cR)]
  model <- lm(formula = log2FoldChange~median_L*median_R,
              data= train)
  pred <- predict(model, newdata = .SD)
  .(rsq= vl_model_eval(observed = log2FoldChange, 
                       predicted = pred)$Rsquare)
}, set]
saveRDS(model, 
        "Rdata/CV_linear_model_vllib002.rds")

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat[, predicted:= predict(model)]
# Select examples
ex <- dat[list("dev_strong_B_00302",
               c("dev_medium_C_00480", "dev_medium_C_00466", "dev_medium_B_00585")), on= c("L", "R")]
ex[, Cc:=  c("tomato", "limegreen", "cornflowerblue")]

pdf("pdf/draft/CV_linear_model_vllib002.pdf",
    height = 3,
    width = 3)
par(las= 1,
    mar= c(3,3,1,1),
    mgp= c(1.5, 0.5, 0),
    tcl= -0.2,
    cex= 1,
    bty= "n")
# Additive
dat[, {
  smoothScatter(predicted,
                log2FoldChange,
                colramp = colorRampPalette(c("white", gray.colors(3, rev= T))),
                xlab= "Predicted activity (log2)",
                ylab= "Observed activity (log2)")
  abline(0,1,lty=2)
  text(par("usr")[1],
       par("usr")[4]-strheight("M"),
       paste0("PCC= ", round(cor.test(predicted, log2FoldChange)$estimate, 2)),
       pos= 4)
}]

# Examples
ex[, {
  points(predicted,
         log2FoldChange,
         col= Cc, 
         pch= 16)
}]
par(mfrow= c(4,1),
    mar= c(0,0,0,0))
plot.new()
segments(0.25, 
         seq(0.3, 0.7, length.out=3), 
         0.75,
         seq(0.3, 0.7, length.out=3))
rect(rep(c(0.3, 0.55), each= 3),
     rep(seq(0.3, 0.7, length.out= 3)-0.06, 2),
     rep(c(0.3, 0.55), each= 3)+0.15,
     rep(seq(0.3, 0.7, length.out=3)+0.06, 2),
     col= c(rep("gold", 3), ex$Cc))
par(mar= c(2.75,8,0.1,8),
    cex.axis= 0.6,
    mgp= c(1, 0.3, 0),
    cex.axis= 0.6,
    cex= 0.66,
    tcl= -0.1)
ex[, {
  bar <- barplot(c(log2FoldChange,
                   median_L, 
                   median_R,
                   additive,
                   multiplicative,
                   predicted), 
                 col= c("black", "gold", Cc, grey.colors(3, start= 0.9, end= 0.4)),
                 xaxt= "n",
                 ylab= "Activity",
                 space= 1,
                 border= NA,
                 width= 0.5)
  segments(x0 = bar[1]-0.25,
           y0 = log2FoldChange,
           x1 = bar[6]+0.25,
           y1 = log2FoldChange,
           lty= "11")
  text(bar,
       par("usr")[3]-strheight("M", cex= 0.5),
       c("Obs.", "5'", "3'", "Exp. Add.", "Exp. Mult.", "Exp. lm"),
       pos= 2,
       offset= -0.25,
       xpd= T,
       srt= 45,
       cex= 0.8)
}, (ex)]
dev.off()