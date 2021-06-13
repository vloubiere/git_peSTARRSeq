# Import motif counts
.feat <- readRDS("Rdata/master_lib_features.rds")
cols <- c("ID", grep("^motif", colnames(.feat), value= T))
feat <- .feat[, ..cols]

# Import vllib002 STARR-Seq data
dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")
dat <- dat[!(spike_in) & lib=="libvl002" & median_L>1 & median_R>1]
dat <- merge(dat, feat, by.x= "L", by.y= "ID")
dat <- merge(dat, feat, by.x= "R", by.y= "ID", suffixes= c("_L", "_R"))

# Compute formulas and corresponding models
pw <- CJ(grep("motif.*_L", colnames(dat), value = T), grep("motif.*_R", colnames(dat), value = T))
pw <- paste0(pw[, paste0(V1, "*", V2)], collapse = "+")
.f <- list('Additive_exp'= "log2FoldChange~add",
           'L'= "log2FoldChange~median_L",
           'R'= "log2FoldChange~median_R",
           'L+R'= "log2FoldChange~median_L+median_R",
           'L*R'= "log2FoldChange~median_L*median_R",
           'L*R+A'= "log2FoldChange~median_L*median_R+add",
           'L*R*A'= "log2FoldChange~median_L*median_R*add",
           'mL1*mR1+mL2*mR1...'= paste0("log2FoldChange~", pw),
           'L+R+mL1*mR1+mL2*mR1...'= paste0("log2FoldChange~median_L+median_R+", pw))
if(!file.exists("Rdata/linear_models_prediction.rds"))
{
  .mod <- lapply(.f, function(x) lm(formula = as.formula(x), data = dat))
  names(.mod) <- names(.f)
  saveRDS(.mod, "Rdata/linear_models_prediction.rds")
}else
{
  .mod <- readRDS("Rdata/linear_models_prediction.rds")
}

# Boxplot additive model
pdf("pdf/boxplot_add_observed.pdf", width = 2, height = 5)
par(las= 2)
vl_boxplot(dat[, .(Additive= add,
                   Observed= log2FoldChange)], 
           col= "black",
           ylab= "Activity (log2)",
           ylim= c(1,10))
segments(1,8.75, 2, 8.75)
text(1.5, 
     8.75,
     pos=3, 
     labels = cut(cor.test(dat$add, dat$log2FoldChange)$p.value, 
                  c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), 
                  labels = c("****","***","**","*","N.S")))
dev.off()

pdf("pdf/linear_models_peSTARRSeq.pdf", width = 11, height = 12)
par(mfrow= c(3, 3), las= 1)
# Additive score
smoothScatter(dat$add, 
              dat$log2FoldChange,
              main= "additive", 
              xlim= c(1, 10),
              ylim= c(1, 10), 
              xlab= "Additive (log2)",
              ylab= "Observed (log2)")
pcc <- paste("PCC=", round(cor.test(dat$add, dat$log2FoldChange)$estimate, 2))
rsq <- paste("R2=", round(summary(.mod$Additive_exp)$r.squared, 2))
abline(.mod$Additive_exp, lty= 2)
legend("topleft", 
       bty= "n", 
       legend = c(pcc, rsq))
abline(0,1)
inset <- c(grconvertX(c(7.25, 10), "user", "ndc"), 
           grconvertY(c(1, 3.75), "user", "ndc"))
lapply(seq(.mod)[-1], function(i)
{
  x <- predict(.mod[[i]])
  smoothScatter(x, 
                dat$log2FoldChange, 
                main= names(.f)[i],
                xlim= c(1, 10),
                ylim= c(1, 10), 
                xlab= "Predicted (log2)",
                ylab= "Observed (log2)")
  abline(0,1)
  pcc <- paste("PCC=", round(cor.test(x, dat$log2FoldChange)$estimate, 2))
  rsq <- paste("R2=", round(summary(.mod[[i]])$r.squared, 2))
  legend("topleft", 
         bty= "n", 
         legend = c(pcc, rsq))
})
par(fig = inset, 
    mar= c(0,0,0,0),
    xaxs= "i", 
    yaxs= "i", 
    new = T)
set.seed(1)
.shuf <- sample(dat$add, nrow(dat))
smoothScatter(.shuf, 
              dat$log2FoldChange,
              xlim= c(1, 10),
              ylim= c(1, 10), 
              xlab= NA,
              ylab= NA,
              xaxt= "n",
              yaxt= "n")
mtext(paste0("Shuffled PCC= ", 
             round(cor.test(.shuf, dat$log2FoldChange)$estimate, 3)), 
      xpd= T, 
      line= 0.25,
      cex= 0.5)
dev.off()

#-------------------------------------------------------#
# Identify motifs associated with higher/lower residuals
#-------------------------------------------------------#
dat[, residuals:= .mod[["L+R"]]$residuals]
.fish <- data.table(motif= grep("^motif", colnames(dat), value = T))
.fish <- .fish[, .(cdition= c("up", "down")), motif]
.fish[, c("estimate", "pval"):= {
  if(cdition=="up")
    .class <- dat$residuals>1.5
  else if(cdition=="down")
    .class <- dat$residuals<(-1.5)
  .c <- fisher.test(dat[[motif]]>0, .class)
  .(log2(.c$estimate), .c$p.value)
}, .(motif, cdition)]
.fish[, padj:= p.adjust(pval, method = "fdr"), cdition]
.fish[, motif_id:= gsub("motif__|_L$|_R$", "", motif), motif]
.fish[readRDS("Rdata/som_enriched_motifs.rds")$info, Dmel:= i.Dmel, on= "motif_id==motif"]
.fish <- .fish[padj<0.001 
               & estimate>=1 
               & Dmel!=""]
# Plot
pdf("pdf/fisher_residuals.pdf", width = 6.5, height = 4.5)
layout(matrix(c(1,1,2,3), ncol= 2), widths = c(1, 0.6))
plot(NA, 
     xlim= c(0, nrow(dat)),
     ylim= range(dat$residuals),
     ylab= "Best model residuals", 
     xlab= "index",
     las= 1)
rect(par("usr")[1], 
     par("usr")[3], 
     par("usr")[2], 
     -1.5,
     border= NA, 
     col= adjustcolor("cornflowerblue", 0.5))
rect(par("usr")[1], 
     1.5, 
     par("usr")[2], 
     par("usr")[4],
     border= NA, 
     col= adjustcolor("tomato", 0.5))
lines(sort(dat$residuals))
abline(h= c(-1.5, 1.5))
par(mar= c(4,4,2,2))
for(resid in c("up", "down"))
{
  if(resid=="up")
  {
    .c <- .fish[cdition=="up"][order(estimate)]
    Cc <- adjustcolor("tomato", 0.5)
  }
  else if(resid=="down")
  {
    .c <- .fish[cdition=="down"][order(estimate)]
    Cc <- adjustcolor("cornflowerblue", 0.5)
  }
  bar <- barplot(.c$estimate, 
                 ylab= "Odds ratio (log2)", 
                 names.arg = gsub("__", "\n", .c$Dmel), 
                 border= NA,
                 las= 2, 
                 cex.names= 0.6, 
                 col= Cc)
}
dev.off()

#-------------------------------------------------------#
# Identify Enhancers for which motifs are informative
#-------------------------------------------------------#
x <- .mod[["L+R"]]$residuals # Compare residuals of best model with and without motifs
y <- .mod[["L+R+mL1*mR1+mL2*mR1..."]]$residuals
diff <- .mod[["L+R+mL1*mR1+mL2*mR1..."]]$residuals-.mod[["L+R"]]$residuals

# Fraction of enhancer group between all and the ones with lowered residuals
all <- copy(dat)
all[.feat, group_L:= i.group, on= "L==ID"]
all[.feat, group_R:= i.group, on= "R==ID"]
sub <- all[diff>1]
res <- rbindlist(list(all= all[, .N/nrow(all)*100, .(group_L, group_R)], 
                      improved= sub[, .N/nrow(sub)*100, .(group_L, group_R)]), idcol = T)
res <- dcast(res, group_L+group_R~.id, value.var = "V1")
res <- as.matrix(res[, 3:4], rownames= res[, paste0(group_L, "+", group_R)])
Cc <- adjustcolor(c("blue", 
                    "darkorchid4",
                    "red3",
                    "red2", 
                    "indianred2", 
                    "pink"), 0.5)

# PLOT
pdf("pdf/enhancers_improved_prediction_with_motifs.pdf", width = 5.75, height = 2.75)
layout(matrix(1:3, nrow= 1), widths = c(1, 0.6, 0.8))
plot(x, 
     y,
     xlab= "L+R w/o motifs",
     ylab= "L+R+motifs",
     col= adjustcolor(ifelse(diff>1, "tomato", "grey"), 0.5),
     pch= 19,
     las= 1)
par(las= 2)
vl_boxplot(data.table("L+R"= x[diff>1], "L+R+mot."= y[diff>1]), 
           col= "tomato", 
           ylab= "residuals")
par(mar= c(5.1,4.1,4.1,6))
barplot(res, 
        col= Cc)
legend(2.5, 
       105, 
       bty= "n", 
       legend = rownames(res), 
       fill= Cc, 
       xpd= T,
       cex= 0.7)
dev.off()

