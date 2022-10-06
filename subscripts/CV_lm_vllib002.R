setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & !grepl("^control", L) & !grepl("^control", R)]
set.seed(1)
dat <- lib[sample(nrow(lib), nrow(lib))]

#-----------------------------------------------#
# Train model
#-----------------------------------------------#
if(!file.exists("Rdata/CV_linear_model_vllib002.rds"))
{
  # Define train and test sets
  set.seed(1)
  dat[dat[, .(set= sample(5)), L], setL:= i.set, on= "L"]
  set.seed(1)
  dat[dat[, .(set= sample(5)), R], setR:= i.set, on= "R"]
  dat[, set:= paste0("test", .GRP), .(setL, setR)]
  
  # Train linear model for each train set and compute predicted values
  model <- lm(formula = log2FoldChange~median_L*median_R,
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
} else
  model <- readRDS("Rdata/CV_linear_model_vllib002.rds")


#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat[, predicted:= predict(model)]
# Select examples
ex <- dat[list("dev_medium_B_00524",
               c("dev_weak_C_00365", "dev_strong_B_00302", "dev_medium_C_00542")), on= c("L", "R")]
ex[, Cc:=  c("cornflowerblue", "limegreen", "tomato")]

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
dev.off()